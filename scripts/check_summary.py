"""Verify every number quoted in the collaborator summary against the JSON result files.
Markers: <!--num:FILE#PATH-->VALUE (see tests/acceptance/m6/test_m6_package.py).
Usage: python3 scripts/check_summary.py [--summary PATH]; exit 0 iff all markers match."""
import argparse
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
MARKER = re.compile(r"<!--num:([^#\s>]+)#(\S*?)-->([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")
ANY = re.compile(r"<!--num:")


def _match_field(v, s):
    if isinstance(v, bool):
        return s in ("true", "false") and v == (s == "true")
    if isinstance(v, (int, float)):
        try:
            return float(s) == float(v)
        except ValueError:
            return False
    return str(v) == s


def resolve(obj, path):
    for seg in path.split("/") if path else []:
        if isinstance(obj, dict) and seg in obj:
            obj = obj[seg]
        elif isinstance(obj, list) and re.fullmatch(r"-?\d+", seg):
            obj = obj[int(seg)]
        elif isinstance(obj, list) and "=" in seg:
            conds = [c.split("=", 1) for c in seg.split(",")]
            hits = [o for o in obj if isinstance(o, dict)
                    and all(k in o and _match_field(o[k], v) for k, v in conds)]
            if len(hits) != 1:
                raise KeyError(f"selector {seg!r} matched {len(hits)} items")
            obj = hits[0]
        else:
            raise KeyError(f"segment {seg!r} does not resolve")
    return obj


def tolerance(text):
    m = re.fullmatch(r"[-+]?(\d*)\.?(\d*)(?:[eE]([-+]?\d+))?", text)
    d = len(m.group(2)) if "." in text.split("e")[0].split("E")[0] else 0
    e = int(m.group(3)) if m.group(3) else 0
    return 0.5 * 10 ** (e - d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=str(ROOT / "reports" / "collaborator_summary.md"))
    text = pathlib.Path(ap.parse_args().summary).read_text()
    marks = MARKER.findall(text)
    bad = []
    if len(ANY.findall(text)) != len(marks):
        bad.append("malformed marker(s) present")
    cache = {}
    for f, path, val in marks:
        try:
            if f not in cache:
                cache[f] = json.loads((ROOT / f).read_text())
            v = resolve(cache[f], path)
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise TypeError(f"not a number: {v!r}")
            if abs(float(val) - v) > tolerance(val) * (1 + 1e-9):
                bad.append(f"{f}#{path}: summary says {val}, JSON has {v}")
        except Exception as exc:  # noqa: BLE001
            bad.append(f"{f}#{path}: {exc}")
    if not marks:
        bad.append("no <!--num:...--> markers found")
    for b in bad:
        print("FAIL", b)
    print(f"{len(marks)} markers, {len(bad)} failures")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
