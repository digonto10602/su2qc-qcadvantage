# CLAUDE.md — su2qc-qcadvantage (read fully at session start)

You are building the code for the "SU(2) torus frontier run". The specification is `docs/METHODS.md`, the milestones are `plan/PLAN.md` (M1–M6), and the progress record is `notes/STATUS.md`. Physics correctness comes first. Credit is limited, so work efficiently (section 4).

## 0. First session only: setup check
- If `.claude/settings.json` does not exist (hidden folders can be lost in a browser upload), copy `claude-setup-backup/` to `.claude/` (`settings.json`, `agents/`, `hooks/`), copy `claude-setup-backup/gitignore.txt` to `.gitignore`, commit, push, then STOP and tell Digonto: "Setup restored — please start a new session."
- Run `python3 -c "import numpy, scipy, qiskit, pytest"`. If that fails, run `python3 -m pip install --user -r requirements.txt` (add `--break-system-packages` if pip refuses).
- Run `python3 tools/check_lock.py`. It must print `LOCK OK`.

## 1. Roles (the team)
| Role | Who | Model | When |
|---|---|---|---|
| Builder + coordinator | you (main session) | Sonnet 5.5 | all routine work: code, unit tests, runs, commits, STATUS |
| `reviewer` | subagent | Opus 5.5 | (a) start of M4, M5, M6: writes acceptance tests from PLAN criteria; (b) every milestone gate: independent physics/numerics audit |
| `fixer` | subagent | Opus 5.5 | only on escalation (section 3) |

No other subagents. Never use the built-in general-purpose or Explore agents.

## 2. Build loop (per milestone)
1. Read only that milestone's section of `plan/PLAN.md` and the relevant METHODS sections.
2. M4–M6 only: ask `reviewer` to write `tests/acceptance/mN/` from the criteria. Run `python3 tools/lock.py mN`, then commit "lock mN".
3. Implement in small steps, writing unit tests in `tests/unit/` as needed. After each step run only the relevant tests: `python3 -m pytest tests/acceptance/mN -x -q`.
4. Commit after each passing step: `git add -A && git commit -m "mN: <what>" && git push`.
5. Milestone gate:
   - run the full suite `python3 -m pytest -q`, which must be green
   - ask `reviewer` for a gate audit, passing file paths and a 5-line summary (not file contents)
   - fix blocking findings
   - `git tag mN-done && git push --tags`
   - update STATUS
6. Go straight on to the next milestone. Do not wait for Digonto unless a stop rule fires.

## 3. Escalation and stop rules
- Count attempts per failing acceptance test in STATUS ("attempt log").
- After 3 failed attempts on the same test, call `fixer` once with the test name, the error (last 30 lines) and file paths.
- If it still fails after `fixer`, STOP.
- Also STOP if:
  - a locked test looks wrong (never edit it; write the evidence in STATUS)
  - a package cannot be installed
  - a job would exceed 12 GB RAM or 2 h CPU
  - G2, G3, G5 or G8 fails after `fixer` (this means the model or cost claim in `docs/idea.md` is wrong)
- When stopping, write the reason and the evidence file under "Decisions waiting for Digonto" in STATUS, commit, push, then end your turn with a 5-line summary.
- STOP after M6 is complete (tag `v0.1-pre-cluster`).
- Locked files (tests/acceptance m1..m3 and any locked milestone, `_ref.py`, `conftest.py`, `LOCK.json`, `tools/check_lock.py`, `tools/lock.py`, `.claude/`) are protected by a hook. Never try to get around it.

## 4. Credit-saving rules (important)
- Do routine work yourself. Use subagents only as in section 1. When you brief a subagent, give file paths, not contents, and ask for a verdict of at most 300 words. It writes long output to `reports/`.
- Read files once per session. Use `grep -n` / `sed -n 'a,bp'` to look at parts. Don't print large arrays, JSON or logs: use `| tail -20` and `head -c 2000`. Don't read `docs/idea.md` unless the spec is unclear.
- Run targeted tests with `-x -q`. Run the full suite only at gates.
- Long computations: one batch script, `nohup python3 scripts/X.py > logs/X.log 2>&1 &`, then poll with `sleep 590; tail -5 logs/X.log` (each poll is one cheap turn). Never poll more often than every 5 minutes.
- Keep `notes/STATUS.md` under 60 lines. Use compact commit messages.
- No web access (blocked). Everything needed is in `docs/`.

## 5. Conventions
- Python 3, numpy/scipy/qiskit; package `su2torus/`, scripts in `scripts/`, results in `results/mN/`, reports in `reports/`.
- Follow METHODS exactly (bit order, Trotter ordering, gate-count convention). The gate count of ~9.5 per plaquette per step and the Pauli form are claims to verify, not facts.
- Reports for Digonto: plain English, LaTeX math, define every technical term on first use, order = what was asked → what was done → what it means → problems/assumptions.
- Never force-push or rewrite history. Push to the branch the session is on.
