#!/usr/bin/env python3
"""PreToolUse guard: blocks edits to locked acceptance tests, lock tooling and .claude/.
Exit code 2 blocks the tool call and shows the message to Claude."""
import json
import os
import pathlib
import re
import sys

ROOT = pathlib.Path(os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd())).resolve()
ALWAYS = ["tests/acceptance/_ref.py", "tests/acceptance/conftest.py",
          "tests/acceptance/LOCK.json", "tools/check_lock.py", "tools/lock.py", ".claude/"]


def protected_prefixes():
    out = list(ALWAYS)
    try:
        lock = json.loads((ROOT / "tests/acceptance/LOCK.json").read_text())
        out += [f"tests/acceptance/{m}/" for m in lock.get("milestones", [])]
    except Exception:
        pass
    return out


def rel(path_str):
    p = pathlib.Path(path_str)
    if not p.is_absolute():
        p = ROOT / p
    try:
        return p.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return None


def block(msg):
    print(f"BLOCKED by guard: {msg}. Locked files must not change; if a locked test looks "
          f"wrong, record evidence in notes/STATUS.md and stop.", file=sys.stderr)
    sys.exit(2)


def main():
    try:
        data = json.load(sys.stdin)
    except Exception:
        sys.exit(0)
    tool = data.get("tool_name", "")
    inp = data.get("tool_input", {}) or {}
    prots = protected_prefixes()
    if tool in ("Edit", "Write", "MultiEdit", "NotebookEdit"):
        r = rel(inp.get("file_path") or inp.get("notebook_path") or "")
        if r and any(r == p.rstrip("/") or r.startswith(p) for p in prots):
            block(f"{tool} on protected file {r}")
    elif tool == "Bash":
        cmd = inp.get("command", "")
        if re.search(r"git\s+push\s+.*(--force|-f\b)", cmd):
            block("force push")
        # 1) redirection targets (">", ">>"), ignoring fd redirects like 2>&1
        for target in re.findall(r"(?<![0-9&])>{1,2}\s*([^\s&|;<>]+)", cmd):
            r = rel(target)
            if r and any(r == p.rstrip("/") or r.startswith(p) for p in prots):
                block(f"redirect into protected file {r}")
        # 2) commands that can modify files, if a protected path is mentioned
        writes = re.search(r"(\btee\b|\bsed\s+-i|\bmv\b|\bcp\b|\brm\b|\bchmod\b|\btruncate\b|"
                           r"\bln\b|git\s+(checkout|restore|rm|mv)\b|\bpython3?\s+-c\b|\bperl\b)", cmd)
        if writes:
            for p in prots:
                if p.rstrip("/") in cmd:
                    block(f"shell command touching protected path {p}")
    sys.exit(0)


if __name__ == "__main__":
    main()
