# deploy/

One `docker compose` stack: `app` (plasmid bank), `genomes` (genome browser) and `caddy` (HTTPS,
one hostname). Only Caddy publishes ports, so the apps are reachable only through it.

**Demo:** `deploy/quickstart.sh up`. It creates `deploy/.env` with a random secret, indexes
`examples/bank`, starts the stack and prints a first admin password.

**Your own server:**

1. `cp deploy/.env.example deploy/.env`, then set:
   - `SITE_ADDRESS`: the server's DNS name;
   - `HTTP_PORT=80`, `HTTPS_PORT=443`;
   - `TLS_DIRECTIVE`: your certificate, or empty for Let's Encrypt;
   - a long random `PLASMIDBANK_SECRET`.
2. Point `BANK_MOUNT` at your map folder (mounted read-only) and `OLIGO_FILE` at your oligo list.
3. Index the bank, then start:
   ```bash
   dc() { docker compose -f deploy/docker-compose.yml --env-file deploy/.env "$@"; }
   dc build && dc --profile tools run --rm scan && dc up -d app genomes caddy
   dc exec app python -m plasmidbank.app.auth add <name>        # asks for a password
   dc exec app python -m plasmidbank.app.auth role <name> admin
   ```
4. Schedule the scan, and `sync` if the app should write to a read-write mount of the share
   (`BANK_MOUNT_RW`). `deploy/systemd/` has example units: the scan daily at 07:00, the sync hourly.
   Edit their `WorkingDirectory` and `.mount` names, then:
   ```bash
   sudo cp deploy/systemd/plasmidbank-{scan,sync}.{service,timer} /etc/systemd/system/
   sudo systemctl daemon-reload && sudo systemctl enable --now plasmidbank-scan.timer plasmidbank-sync.timer
   journalctl -u plasmidbank-scan.service -n 50
   ```
   The scan only rewrites the index. The app notices within about two minutes and swaps the new
   list in, so there is no restart. cron works as well.

**What `sync` writes to the share** (nothing, until `BANK_MOUNT_RW` is mounted; uploads then wait):
- **New uploads**, add-only, into the example ID layout (`plasmidbank/uploads.py`, `bin_for`). It never
  overwrites, renames or deletes, and refuses an ID that is already on the share.
- **My workspace**, both ways, with `<share>/Workspaces/<user>/<folder>/`. An edit on one side is
  copied to the other. An edit on both sides keeps both copies. A delete is honoured only if the other
  side did not change. New GenBank/FASTA files people drop in their folder are imported. Nothing is
  deleted: share copies go to `Workspaces/<user>/_deleted/`, app copies to `/state/workspaces_trash/`.
  The scan skips `Workspaces/`.
- **Overwrites** someone confirmed twice in the app. They are written only if the file on the share is
  unchanged since the user saw it. The old file goes to `_replaced/` beside it; to undo, move it back.
- Alerts (share unreachable, conflicts, refusals) go to `PLASMIDBANK_DISCORD_WEBHOOK` if set.

**Label images:** `GET /api/qr/<text>.png` and `GET /api/code128/<text>.png` (`?px=` sets the
size; `<text>` is `[A-Za-z0-9._-]{1,40}`). They return a 1-bit PNG for a lab database or label printer. Such
clients cannot sign in, so these two routes, and nothing else, also accept
`Authorization: Bearer <PLASMIDBANK_LABEL_TOKEN>`. With the token unset they refuse it (503), and a
signed-in browser needs no token. FileMaker example: `Insert from URL [ <container> ;
"https://<your server>/api/qr/" & Table::Number & ".png" ; cURL options: "-H \"Authorization: Bearer " & $token & "\"" ]`.

**Restart notice:** while a service is starting, Caddy serves `deploy/errors/starting.html`. Before a
planned restart, write `{"deploy_at": <unix seconds>}` to `deploy/errors/deploy-notice.json`; every
page then shows a countdown bar. Delete the file afterwards. The pages ignore a notice more than 15 min old.

**Usage log:** the app appends one line per request and one per minute of CPU/RAM/disk to
`/state/usage/`. `deploy/usage_report.py --host <ssh host>` (or `--container lab_server-app-1` for a
local stack) fetches it and builds an HTML report: who uses what, where people wait, and whether the
server needs more cores, RAM or disk. The report holds personal data, so it is written outside the repo.
`deploy/loadtest.py` measures how many people at once the server takes.

**Genome browser:**
- Genomes live in the `genomes` volume. They are uploaded on its page, and the largest single
  upload is `GENOMES_MAX_UPLOAD_GB`. Watch the disk: a gzipped FASTA expands about 4×, and each
  genome keeps its FASTA plus a database.
- Set `GENOMES_IMPORT_DIR` to a folder on the server (e.g. a read-only mount of the lab share) to
  let admins import big genomes without uploading them.
- Sign-in comes from the bank. The `Caddyfile` sends every `/genomes/` request through
  `forward_auth` to `app:8765/auth/check`, inside a `route` block so the order holds.

**Back up** the `state` volume (accounts, projects, workspaces, uploads, usage log) and the `genomes` volume.
