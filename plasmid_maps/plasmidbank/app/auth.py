"""Who is using the app, so each person gets their own workspace (outbox) on a shared server.

PLASMIDBANK_AUTH selects the mode:
  none    single user on a laptop (default); everyone is "local".
  local   accounts in PLASMIDBANK_USERS (JSON, scrypt hashes) + a signed session cookie.
          Manage accounts with:  python -m plasmidbank.app.auth add|remove|list|role <name> [role]
          or, once one admin exists, from the admin page (/admin).
  header  a login proxy in front of the app (e.g. oauth2-proxy with AU single sign-on) sets
          X-Forwarded-User / X-Forwarded-Email. ONLY safe when the app port is reachable from
          the proxy alone (the docker-compose network, never published).
PLASMIDBANK_SECRET signs cookies (required for local mode; any long random string).
Stdlib only.
"""

from __future__ import annotations

import base64
import getpass
import hashlib
import hmac
import json
import os
import re
import secrets
import sys
import time
from pathlib import Path

MODE = os.environ.get("PLASMIDBANK_AUTH", "none").lower()
USERS_FILE = Path(
    os.environ.get("PLASMIDBANK_USERS", Path.home() / "PlasmidBank" / "users.json")
)
SECRET = os.environ.get("PLASMIDBANK_SECRET", "")
COOKIE = "plasmidbank_session"
MAX_AGE = 12 * 3600
_NAME = re.compile(r"^[a-z0-9][a-z0-9._-]{1,40}$")
ROLES = ("admin", "user", "student")  # student limits are not enforced yet (they wait for the approval queue)


def safe_name(user: str) -> str:
    """Filesystem-safe workspace name for a user id or e-mail address."""
    s = re.sub(r"[^a-z0-9._-]+", "_", user.lower()).strip("._")
    return s[:60] or "unknown"


# ── password store ────────────────────────────────────────────────────────────


def _hash(password: str, salt: bytes) -> str:
    return hashlib.scrypt(
        password.encode(), salt=salt, n=2**14, r=8, p=1, dklen=32
    ).hex()


def load_users() -> dict:
    try:
        return json.loads(USERS_FILE.read_text())
    except FileNotFoundError:
        return {}


def save_users(users: dict) -> None:
    USERS_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp = USERS_FILE.with_suffix(".tmp")
    tmp.write_text(json.dumps(users, indent=1))
    os.chmod(tmp, 0o600)
    tmp.replace(USERS_FILE)


def set_password(name: str, password: str, must_change: bool = False) -> None:
    """Set (or create) an account's password; keeps its role and other fields."""
    if not _NAME.match(name):
        raise ValueError("user names: lowercase letters, digits, . _ - (2-41 chars)")
    if len(password) < 10:
        raise ValueError("password must be at least 10 characters")
    users = load_users()
    salt = secrets.token_bytes(16)
    u = users.get(name, {})
    u.update({"salt": salt.hex(), "hash": _hash(password, salt), "must_change": must_change})
    if not must_change:  # chosen by the user themselves; a one-time password does not count
        from datetime import datetime

        u["password_set"] = datetime.now().isoformat(timespec="minutes")
    users[name] = u
    save_users(users)


def check_password(name: str, password: str) -> bool:
    u = load_users().get(name)
    if not u or "hash" not in u:
        _hash(password, b"0" * 16)  # same cost whether or not the user exists
        return False
    return hmac.compare_digest(u["hash"], _hash(password, bytes.fromhex(u["salt"])))


# ── roles and account management (admin page and CLI) ─────────────────────────

AUDIT_FILE = USERS_FILE.with_name("audit.jsonl")


def role_of(name: str) -> str:
    """'admin' | 'user' | 'student'. Single-user mode is admin; accounts without a role are users."""
    if MODE == "none":
        return "admin"
    return load_users().get(name, {}).get("role", "user")


def must_change(name: str) -> bool:
    u = load_users().get(name, {})
    return MODE == "local" and bool(u.get("must_change"))


def audit(actor: str, action: str, target: str = "", detail: str = "") -> None:
    """One JSON line per account change, next to users.json. Never holds a password."""
    from datetime import datetime

    AUDIT_FILE.parent.mkdir(parents=True, exist_ok=True)
    row = {"time": datetime.now().isoformat(timespec="seconds"), "by": actor, "action": action,
           "target": target, "detail": detail}
    with AUDIT_FILE.open("a", encoding="utf-8") as f:
        f.write(json.dumps(row) + "\n")


def _admins(users: dict) -> list[str]:
    return [n for n, u in users.items() if u.get("role") == "admin"]


def one_time_password() -> str:
    return secrets.token_urlsafe(12)  # 16 characters


def create_user(name: str, role: str = "user") -> str:
    """A new account with a one-time password (returned once, never stored in clear)."""
    if role not in ROLES:
        raise ValueError(f"role must be one of {', '.join(ROLES)}")
    if not _NAME.match(name):
        raise ValueError("user names: lowercase letters, digits, . _ - (2-41 chars)")
    if name in load_users():
        raise ValueError(f"{name} already exists")
    if MODE != "local":  # header mode: the login proxy owns passwords; only the role lives here
        users = load_users()
        users[name] = {"role": role}
        save_users(users)
        return ""
    pw = one_time_password()
    set_password(name, pw, must_change=True)
    users = load_users()
    users[name]["role"] = role
    save_users(users)
    return pw


_DANISH = str.maketrans({"æ": "ae", "ø": "oe", "å": "aa", "ß": "ss"})


def username_for(full_name: str) -> str:
    """'Lærke Maria Bjerrum Nielsen' -> 'laerke.nielsen': first + last name, Danish letters spelled
    out, other accents dropped. Middle names and initials are ignored."""
    import unicodedata

    s = unicodedata.normalize("NFKD", full_name.lower().translate(_DANISH))
    s = "".join(c for c in s if not unicodedata.combining(c))
    parts = [re.sub(r"[^a-z0-9-]", "", p) for p in s.split()]
    parts = [p.strip("-") for p in parts if p.strip("-")]
    if not parts:
        raise ValueError(f"no usable name in {full_name!r}")
    return ".".join([parts[0], parts[-1]] if len(parts) > 1 else parts)[:41]


def bulk_create(full_names: list[str], role: str = "user") -> list[tuple[str, str, str]]:
    """One account per name, each with a one-time password the user must change at first sign-in.
    Returns (full name, username, password) rows; an existing account gets 'exists' and is left
    untouched. A username clash between two new people gets a number (alex.smith2)."""
    rows, taken = [], set(load_users())
    for full in (n.strip() for n in full_names):
        if not full:
            continue
        name = username_for(full)
        users = load_users()
        if name in users and users[name].get("full_name", full) == full:
            rows.append((full, name, "exists"))
            continue
        base, i = name, 2
        while name in taken:
            name, i = f"{base}{i}", i + 1
        pw = create_user(name, role)
        users = load_users()
        users[name]["full_name"] = full
        save_users(users)
        audit("cli", "bulk add", name, role)
        taken.add(name)
        rows.append((full, name, pw))
    return rows


def reset_password(name: str) -> str:
    if name not in load_users():
        raise ValueError(f"no user {name}")
    pw = one_time_password()
    set_password(name, pw, must_change=True)
    return pw


def change_own_password(name: str, current: str, new: str) -> None:
    if not check_password(name, current):
        raise ValueError("current password is wrong")
    if new == current:
        raise ValueError("choose a new password")
    set_password(name, new, must_change=False)


def set_role(name: str, role: str) -> None:
    if role not in ROLES:
        raise ValueError(f"role must be one of {', '.join(ROLES)}")
    users = load_users()
    if name not in users:
        raise ValueError(f"no user {name}")
    if role != "admin" and _admins(users) == [name]:
        raise ValueError("the last admin cannot be demoted; make someone else admin first")
    users[name]["role"] = role
    save_users(users)


def remove_user(name: str) -> None:
    users = load_users()
    if name not in users:
        raise ValueError(f"no user {name}")
    if _admins(users) == [name]:
        raise ValueError("the last admin cannot be removed")
    users.pop(name)
    save_users(users)


def note_sign_in(name: str) -> None:
    from datetime import datetime

    users = load_users()
    if name in users:
        users[name]["last_sign_in"] = datetime.now().isoformat(timespec="minutes")
        save_users(users)


def list_users() -> list[dict]:
    return [{"name": n, "role": u.get("role", "user"), "has_password": "hash" in u,
             "must_change": bool(u.get("must_change")), "last_sign_in": u.get("last_sign_in", ""),
             "password_set": u.get("password_set", "")}
            for n, u in sorted(load_users().items())]


# ── session cookie ────────────────────────────────────────────────────────────


def _sign(payload: str) -> str:
    return hmac.new(SECRET.encode(), payload.encode(), hashlib.sha256).hexdigest()


def make_cookie(name: str) -> str:
    payload = base64.urlsafe_b64encode(f"{name}|{int(time.time())}".encode()).decode()
    return f"{payload}.{_sign(payload)}"


def read_cookie(value: str | None) -> str | None:
    if not value or "." not in value or not SECRET:
        return None
    payload, sig = value.rsplit(".", 1)
    if not hmac.compare_digest(sig, _sign(payload)):
        return None
    try:
        name, ts = base64.urlsafe_b64decode(payload.encode()).decode().split("|")
    except Exception:
        return None
    if time.time() - int(ts) > MAX_AGE or name not in load_users():
        return None
    return name


def current_user(headers, cookies) -> str | None:
    """The signed-in user, or None. In 'none' mode everyone is 'local'."""
    if MODE == "none":
        return "local"
    if MODE == "header":
        u = headers.get("x-forwarded-email") or headers.get("x-forwarded-user")
        return safe_name(u) if u else None
    if MODE == "local":
        return read_cookie(cookies.get(COOKIE))
    return None


def check_config() -> None:
    if MODE not in ("none", "local", "header"):
        raise SystemExit(
            f"PLASMIDBANK_AUTH must be none, local or header (got {MODE!r})"
        )
    if MODE == "local" and len(SECRET) < 32:
        raise SystemExit(
            "PLASMIDBANK_SECRET must be set to a random string of at least 32 characters in local mode"
        )
    if MODE == "local" and "change-me" in SECRET:
        raise SystemExit("PLASMIDBANK_SECRET is still the placeholder from .env.example: set a fresh random one")


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if not argv or argv[0] not in ("add", "reset", "remove", "list", "role", "bulk"):
        print("usage: python -m plasmidbank.app.auth add|reset|remove|list [name] | role <name> admin|user|student"
              " | bulk [role] < names.txt")
        return 2
    cmd = argv[0]
    if cmd == "bulk":  # full names on stdin, one per line -> CSV of one-time passwords on stdout
        import csv

        try:
            rows = bulk_create(sys.stdin.read().splitlines(), argv[1] if len(argv) > 1 else "user")
        except ValueError as e:
            print(e, file=sys.stderr)
            return 1
        w = csv.writer(sys.stdout)
        w.writerow(["name", "username", "one_time_password"])
        w.writerows(rows)
        return 0
    if cmd == "list":
        rows = list_users()
        print("\n".join(f"{u['name']:<24} {u['role']}" for u in rows) or "(no users)")
        return 0
    if len(argv) < 2:
        print(f"usage: python -m plasmidbank.app.auth {cmd} <name>")
        return 2
    name = argv[1]
    try:
        if cmd == "role":
            if len(argv) < 3:
                print("usage: python -m plasmidbank.app.auth role <name> admin|user|student")
                return 2
            users = load_users()
            if name not in users and MODE == "header":  # header mode: an entry with a role only
                users[name] = {}
                save_users(users)
            set_role(name, argv[2])
            audit("cli", "role", name, argv[2])
            print(f"{name} is now {argv[2]}")
            return 0
        if cmd == "reset":  # a one-time password to copy and pass on; they choose their own at sign-in
            new = name not in load_users()
            pw = create_user(name) if new else reset_password(name)
            audit("cli", "add" if new else "reset password", name)
            print(f"{'new account' if new else 'one-time password for'} {name}: {pw}")
            print("(shown only now; they choose their own password at their next sign-in)")
            return 0
        if cmd == "remove":
            remove_user(name)
            audit("cli", "remove", name)
            print(f"removed {name}")
            return 0
    except ValueError as e:
        print(e)
        return 1
    pw = os.environ.get("PLASMIDBANK_NEW_PASSWORD")
    if not pw:  # typed: ask twice, since the prompt does not echo
        pw = getpass.getpass(f"password for {name} (10+ characters): ")
        if getpass.getpass("again: ") != pw:
            print("the two passwords differ; nothing changed")
            return 1
    try:
        set_password(name, pw)
    except ValueError as e:  # too short, or a bad user name: say so, no traceback
        print(f"{e}; nothing changed")
        return 1
    audit("cli", "set password", name)
    print(f"saved {name} in {USERS_FILE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
