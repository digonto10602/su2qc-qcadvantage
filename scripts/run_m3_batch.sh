#!/bin/bash
# M3 batch: N=12 (if missing), open string, the six N=24 runs (one process each), assemble.
cd "$(dirname "$0")/.."
[ -f results/m3/parts/n12.json ] || python3 scripts/run_m3.py n12
[ -f results/m3/parts/string_open_3x4_g1.4.json ] || python3 scripts/run_m3.py string
for m in su2 abelian; do for s in vacuum neel_half stripe; do
  [ -f results/m3/parts/n24_${s}_${m}.json ] || python3 scripts/run_m3.py n24 $s $m
done; done
python3 scripts/run_m3.py assemble
echo BATCH DONE
