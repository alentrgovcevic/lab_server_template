"""Last active per user, for the admin page.

The auth middleware calls touch() on every signed-in request. Times are kept in memory and written
to <state>/last_active.json at most once a minute, so a request never waits on a disk write. A
restart loses at most the last minute.
"""

from __future__ import annotations

import json
import threading
import time
from datetime import datetime
from pathlib import Path

FLUSH_SECONDS = 60


class Activity:
    def __init__(self, state: Path):
        self.file = Path(state) / "last_active.json"
        try:
            self.seen: dict[str, str] = json.loads(self.file.read_text())
        except (OSError, ValueError):
            self.seen = {}
        self._dirty = False
        self._flushed = time.monotonic()
        self._lock = threading.Lock()

    def touch(self, user: str) -> None:
        self.seen[user] = datetime.now().isoformat(timespec="minutes")
        self._dirty = True
        if time.monotonic() - self._flushed >= FLUSH_SECONDS:
            self.flush()

    def flush(self) -> None:
        with self._lock:
            if not self._dirty:
                return
            self._dirty, self._flushed = False, time.monotonic()
            try:
                tmp = self.file.with_suffix(".tmp")
                tmp.write_text(json.dumps(self.seen, indent=1))
                tmp.replace(self.file)
            except OSError:
                self._dirty = True

    def get(self, user: str) -> str:
        return self.seen.get(user, "")
