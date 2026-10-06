#!/usr/bin/env bash
# SessionStart: print a compact context block (keeps the coordinator's reload cheap).
P="${CLAUDE_PROJECT_DIR:-$(pwd)}"
cd "$P" || exit 0
export PATH="$HOME/.local/bin:$PATH"
[ -n "${CLAUDE_ENV_FILE:-}" ] && echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$CLAUDE_ENV_FILE"
echo "=== su2qc-qcadvantage session start ==="
sed -n '1,60p' notes/STATUS.md 2>/dev/null
echo "--- lock: $(python3 tools/check_lock.py 2>&1 | head -3 | tr '\n' ' ')"
echo "--- last commits:"; git log --oneline -3 2>/dev/null
if ! python3 -c "import numpy, scipy, qiskit, pytest" 2>/dev/null; then
  echo "--- WARNING: core packages missing -> python3 -m pip install --user -r requirements.txt"
fi
echo "--- Plan: plan/PLAN.md | Spec: docs/METHODS.md | Rules: CLAUDE.md"
exit 0
