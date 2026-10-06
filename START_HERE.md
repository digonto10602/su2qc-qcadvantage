# START HERE — upload to GitHub and run in a Claude Code cloud session

## 1. Create the repo (browser)
github.com → **New repository** → name `su2qc-qcadvantage` → **Public** → tick nothing else (no README, no .gitignore, no license) → **Create repository**.

## 2. Upload the files (browser)
1. Unzip `su2qc-qcadvantage.zip` on your computer and open the `su2qc-qcadvantage` folder.
2. **Show hidden files** so that `.claude` and `.gitignore` are visible:
   - Linux (GNOME Files/Nautilus): `Ctrl+H`
   - macOS Finder: `Cmd+Shift+.`
   - Windows Explorer: View → Show → Hidden items
3. On the empty repo page, click **uploading an existing file**. Select **everything inside** the folder (`Ctrl+A`, including `.claude` and `.gitignore`) and drag it into the browser. Commit to `main`.
4. Check on GitHub that `.claude/settings.json`, `.claude/agents/reviewer.md` and `.claude/hooks/guard.py` exist. If `.claude` is missing, that's OK: the first session restores it from `claude-setup-backup/` and asks you to restart.

## 3. Cloud environment (claude.ai/code)
1. Claim the cloud credit if you haven't yet (deadline 7 Oct 2026).
2. Connect GitHub (authorize the Claude GitHub App for `su2qc-qcadvantage`).
3. Create an environment: **Network access: Trusted**. As the **Setup script**, paste:

```bash
#!/bin/bash
pipi() { python3 -m pip install --user --quiet "$@" || python3 -m pip install --user --quiet --break-system-packages "$@"; }
pipi numpy scipy qiskit qiskit-aer matplotlib pytest quimb
exit 0
```

## 4. Start the session
Select repo `su2qc-qcadvantage` and the environment, then send:

```
Follow CLAUDE.md. Do the section-0 setup check, then build milestones M1 to M6 from
plan/PLAN.md in order, using the build loop and credit-saving rules. Continue without
asking me until a stop rule in CLAUDE.md section 3 fires or M6 is tagged v0.1-pre-cluster.
```

## 5. When it stops
Open `notes/STATUS.md` on GitHub and read "Decisions waiting for Digonto". Reply in the same session and say "continue".
