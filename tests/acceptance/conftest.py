"""Acceptance-test integrity check (LOCKED). Fails the run if any locked file changed."""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.check_lock import check  # noqa: E402


def pytest_sessionstart(session):
    problems = check(ROOT)
    if problems:
        raise SystemExit("ACCEPTANCE LOCK BROKEN:\n  " + "\n  ".join(problems))
