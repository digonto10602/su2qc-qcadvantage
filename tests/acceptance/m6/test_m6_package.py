"""M6: collaborator package -- reports/collaborator_summary.md, figures, scripts/check_summary.py,
results/m6/roadmap.json and the README reproduction section.

PLAN M6 acceptance: "every figure is referenced and its data file exists"; "the numbers quoted in
the summary are reproduced from the JSON files by scripts/check_summary.py"; README explains how
to reproduce everything.

1. Figures (fixed names, PNG, written by scripts/make_figures.py), each embedded in the summary as
   a Markdown image with a path relative to reports/, e.g. ![...](figures/fig1_exact_dynamics.png).
   The CAPTION of a figure = the first paragraph after the image line (skip blank lines, then the
   non-blank lines up to the next blank line). It must name every required data file below by its
   repo-relative path, e.g. "Data: `results/m3/summary.json`".
     (i)   reports/figures/fig1_exact_dynamics.png   <- results/m3/summary.json
           (exact N = 12 and 24 dynamics with thermal and decohered lines)
     (ii)  reports/figures/fig2_noise_mitigation.png <- results/m4/noise_runs.json
     (iii) reports/figures/fig3_resources.png        <- results/m5/resources.json
     (iv)  reports/figures/fig4_mps_convergence.png  <- results/m6/convergence.json,
                                                        results/m6/n48_stripe_g1.25.json
     (v)   reports/figures/fig5_cluster_roadmap.png  <- results/m6/roadmap.json
   Any other image referenced in the summary must also exist.

2. results/m6/roadmap.json -- what remains for the GPU cluster and hardware (figure v):
   {"items": [{"key": str, "title": str, "N": [int, ...], "purpose": str, "needs": str}, ...]}
   keys required: "mps_large_chi", "peps_bp", "pauli_propagation", "qmc_thermal_anchor"
   (all its N in 48..64), "hardware_runs".

3. Quoted numbers. Every number the summary quotes from a result file is tagged INLINE as
       <!--num:FILE#PATH-->VALUE
   FILE  repo-relative path of a .json file;
   PATH  '/'-separated segments, no whitespace; each segment is a dict key, a list index (int),
         or a selector "k1=v1,k2=v2" picking the UNIQUE dict in a list whose fields match
         (numbers compared numerically, true/false as booleans, otherwise strings);
   VALUE the visible number, immediately after "-->": [-+]?digits[.digits][e[-+]digits].
   A VALUE matches the JSON number v (int/float, not bool) when it equals v rounded to the
   precision shown: |VALUE - v| <= 0.5 * 10^(e - d) (d = digits after the point, e = exponent).
   Example:  the chi=256 run is converged up to step <!--num:results/m6/convergence.json#last_converged_step/2-->5.
   Required: >= 1 marker into each of results/m3/summary.json, results/m4/noise_runs.json,
   results/m5/resources.json, results/m6/convergence.json, and the headline marker
   results/m6/convergence.json#last_converged_step/2.

   scripts/check_summary.py [--summary PATH]     (default reports/collaborator_summary.md;
       FILE paths are relative to the repo root, not to PATH)
   Exit 0 iff the summary has >= 1 marker and every marker is well-formed, resolves, and matches;
   otherwise non-zero, printing one line per failing marker.

4. Summary text: LaTeX math ($...$), and a section whose heading contains "Glossary" or
   "Definitions" that defines (at least) the terms in GLOSSARY below.

5. README.md: a section whose heading contains "Reproduc" (until the next heading of the same or
   higher level) containing the strings in README_MUST; every scripts/... and configs/... path
   named anywhere in README.md exists.
"""
import json
import pathlib
import re
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[3]
SUMMARY = ROOT / "reports" / "collaborator_summary.md"
FIGURES = {
    "i": ("reports/figures/fig1_exact_dynamics.png", ["results/m3/summary.json"]),
    "ii": ("reports/figures/fig2_noise_mitigation.png", ["results/m4/noise_runs.json"]),
    "iii": ("reports/figures/fig3_resources.png", ["results/m5/resources.json"]),
    "iv": ("reports/figures/fig4_mps_convergence.png",
           ["results/m6/convergence.json", "results/m6/n48_stripe_g1.25.json"]),
    "v": ("reports/figures/fig5_cluster_roadmap.png", ["results/m6/roadmap.json"]),
}
ROADMAP_KEYS = {"mps_large_chi", "peps_bp", "pauli_propagation", "qmc_thermal_anchor", "hardware_runs"}
GLOSSARY = ["plaquette", "Trotter", "matrix product state", "bond dimension", "fidelity",
            "decohered", "thermal", "Loschmidt echo", "zero-noise extrapolation", "PEPS",
            "belief propagation", "Pauli propagation", "quantum Monte Carlo"]
README_MUST = ["requirements.txt", "run_m3", "scripts/run_m4.py", "scripts/run_m5_check.py",
               "scripts/resource_table.py",
               "scripts/run_classical.py --config configs/n48_stripe_g1.25.json",
               "scripts/analyze_convergence.py", "scripts/make_figures.py",
               "scripts/check_summary.py", "pytest", "results/m3", "results/m4", "results/m5",
               "results/m6", "reports/collaborator_summary.md"]
MARKER = re.compile(r"<!--num:([^#\s>]+)#(\S*?)-->([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")
PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


@pytest.fixture(scope="module")
def summary():
    assert SUMMARY.exists(), f"missing {SUMMARY}"
    return SUMMARY.read_text()


def _image_refs(text):
    return [(i, m.group(1)) for i, line in enumerate(text.splitlines())
            for m in re.finditer(r"!\[[^\]]*\]\(([^)\s]+)\)", line)]


def _caption(lines, i):
    j = i + 1
    while j < len(lines) and not lines[j].strip():
        j += 1
    cap = []
    while j < len(lines) and lines[j].strip():
        cap.append(lines[j])
        j += 1
    return lines[i] + "\n" + "\n".join(cap)


# ------------------------------------------------------------------ figures
@pytest.mark.parametrize("fig", list(FIGURES))
def test_figure_referenced_with_data(summary, fig):
    """PLAN M6: every figure is referenced and its data file exists (named in its caption)."""
    png, data = FIGURES[fig]
    path = ROOT / png
    blob = path.read_bytes() if path.exists() else b""
    # > 5 kB: a blank matplotlib canvas is ~2-3 kB, so this rejects empty placeholders
    assert blob[:8] == PNG_MAGIC and len(blob) > 5000, f"missing or empty figure {png}"
    lines = summary.splitlines()
    hits = [i for i, ref in _image_refs(summary)
            if (SUMMARY.parent / ref).resolve() == path.resolve()]
    assert hits, f"figure ({fig}) {png} not embedded in {SUMMARY.name}"
    cap = _caption(lines, hits[0])
    for d in data:
        assert (ROOT / d).exists(), f"data file {d} of figure ({fig}) missing"
        assert d in cap, f"caption of figure ({fig}) does not name its data file {d}"


def test_no_dangling_images(summary):
    for _, ref in _image_refs(summary):
        assert (SUMMARY.parent / ref).exists(), f"image {ref} referenced but missing"


def test_figure_i_data_has_reference_lines():
    """Figure (i) needs exact N=12 and N=24 dynamics with thermal and decohered values."""
    d = json.loads((ROOT / "results" / "m3" / "summary.json").read_text())
    for N in (12, 24):
        runs = [r for r in d["runs"] if r["N"] == N]
        assert runs, N
        assert all("thermal" in r and "infinite_temperature" in r for r in runs), N


def test_roadmap():
    """Figure (v): large-chi MPS, PEPS-BP, Pauli propagation, QMC anchor at N=48-64, hardware."""
    d = json.loads((ROOT / "results" / "m6" / "roadmap.json").read_text())
    items = {it["key"]: it for it in d["items"]}
    assert ROADMAP_KEYS <= set(items), ROADMAP_KEYS - set(items)
    for k in ROADMAP_KEYS:
        it = items[k]
        assert it["title"] and it["purpose"] and it["needs"], k
        assert it["N"] and all(isinstance(n, int) and n > 0 for n in it["N"]), k
    assert all(48 <= n <= 64 for n in items["qmc_thermal_anchor"]["N"])


# ------------------------------------------------------------------ quoted numbers
def _resolve(obj, path):
    for seg in path.split("/") if path else []:
        if isinstance(obj, dict):
            assert seg in obj, f"no key {seg!r}"
            obj = obj[seg]
        elif isinstance(obj, list) and re.fullmatch(r"\d+", seg):
            obj = obj[int(seg)]
        elif isinstance(obj, list):
            conds = [c.split("=", 1) for c in seg.split(",")]

            def ok(el):
                if not isinstance(el, dict):
                    return False
                for k, v in conds:
                    if k not in el:
                        return False
                    f = el[k]
                    if isinstance(f, bool):
                        if v.lower() != str(f).lower():
                            return False
                    elif isinstance(f, (int, float)):
                        try:
                            if abs(float(v) - f) > 1e-12 * max(1.0, abs(f)):
                                return False
                        except ValueError:
                            return False
                    elif str(f) != v:
                        return False
                return True

            hits = [el for el in obj if ok(el)]
            assert len(hits) == 1, f"selector {seg!r} matches {len(hits)} elements"
            obj = hits[0]
        else:
            raise AssertionError(f"cannot descend into {type(obj).__name__} with {seg!r}")
    return obj


def _matches(token, v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return False
    m = re.fullmatch(r"[-+]?(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?", token)
    d = len(m.group(2))
    e = int(m.group(3)) if m.group(3) else 0
    half = 0.5 * 10.0 ** (e - d)
    return abs(float(token) - v) <= half * (1 + 1e-9) + 1e-12 * max(1.0, abs(v))


def test_quoted_numbers_independent(summary):
    """Reviewer's own parser: every marker is well-formed, resolves and matches its JSON value."""
    marks = MARKER.findall(summary)
    assert summary.count("<!--num:") == len(marks), "malformed <!--num:...--> marker(s)"
    files = {f for f, _, _ in marks}
    for need in ("results/m3/summary.json", "results/m4/noise_runs.json",
                 "results/m5/resources.json", "results/m6/convergence.json"):
        assert need in files, f"no quoted number from {need}"
    assert ("results/m6/convergence.json", "last_converged_step/2") in {(f, p) for f, p, _ in marks}
    cache = {}
    for f, p, val in marks:
        if f not in cache:
            assert (ROOT / f).exists(), f
            cache[f] = json.loads((ROOT / f).read_text())
        v = _resolve(cache[f], p)
        assert _matches(val, v), f"{f}#{p}: quoted {val}, JSON {v}"


def _check(*args):
    return subprocess.run([sys.executable, "scripts/check_summary.py", *args], cwd=ROOT,
                          capture_output=True, text=True, timeout=120)


def test_check_summary_passes():
    """PLAN M6: scripts/check_summary.py reproduces the quoted numbers (exit 0)."""
    p = _check()
    assert p.returncode == 0, (p.stdout + p.stderr)[-2000:]


def test_check_summary_detects_errors(summary, tmp_path):
    """The checker must fail on a wrong value, an unresolvable path, and a marker-free file."""
    m = MARKER.search(summary)
    assert m
    bad_val = f"{abs(float(m.group(3))) * 3 + 7:.6g}"
    wrong = summary[:m.start(3)] + bad_val + summary[m.end(3):]
    unres = summary[:m.start(2)] + "no_such_key_xyz" + summary[m.end(2):]
    for name, text in (("wrong", wrong), ("unresolved", unres), ("empty", "# nothing quoted\n")):
        f = tmp_path / f"{name}.md"
        f.write_text(text)
        p = _check("--summary", str(f))
        assert p.returncode != 0, f"check_summary.py accepted the {name} summary"


# ------------------------------------------------------------------ text requirements
def test_summary_math_and_glossary(summary):
    assert re.search(r"\$[^$\n]+\$", summary), "no LaTeX math"
    heads = [(i, line) for i, line in enumerate(summary.splitlines()) if re.match(r"#+\s", line)]
    gl = [i for i, line in heads if re.search(r"glossary|definitions", line, re.I)]
    assert gl, "no Glossary/Definitions section"
    lines = summary.splitlines()
    start = gl[0]
    level = len(re.match(r"#+", lines[start]).group(0))
    end = next((i for i, line in heads if i > start and len(re.match(r"#+", line).group(0)) <= level),
               len(lines))
    body = "\n".join(lines[start:end]).lower()
    missing = [t for t in GLOSSARY if t.lower() not in body]
    assert not missing, f"glossary does not define: {missing}"


def test_readme_reproduction():
    text = (ROOT / "README.md").read_text()
    lines = text.splitlines()
    heads = [(i, len(re.match(r"#+", l).group(0))) for i, l in enumerate(lines) if re.match(r"#+\s", l)]
    rep = [(i, lv) for i, lv in heads if re.search(r"reproduc", lines[i], re.I)]
    assert rep, "README has no reproduction section"
    start, level = rep[0]
    end = next((i for i, lv in heads if i > start and lv <= level), len(lines))
    section = "\n".join(lines[start:end])
    missing = [s for s in README_MUST if s not in section]
    assert not missing, f"README reproduction section lacks: {missing}"
    for ref in set(re.findall(r"\b((?:scripts|configs)/[\w./-]+\.(?:py|sh|sbatch|json))", text)):
        assert (ROOT / ref).exists(), f"README names missing file {ref}"
