#!/bin/bash
# Recreate the milestone tags at the exact commits used in the cloud session and push them.
# Run from a clone that can push tags (e.g. your laptop):  bash scripts/push_tags.sh
# The cloud session's git proxy only allows branch pushes, so these tags never reached GitHub.
set -euo pipefail
git fetch origin claude/modest-mayer-kh4s6b
while read -r tag sha; do
  git cat-file -e "$sha^{commit}"            # fails loudly if the commit is missing
  git tag -f "$tag" "$sha"
done <<'TAGS'
m1-done 252712a4369f
m2-done 42e03a7730f0
m3-done 42e03a7730f0
m4-done 42e03a7730f0
m5-done 42e03a7730f0
m6-done a049a011cb4b
v0.1-pre-cluster a049a011cb4b
TAGS
git push origin m1-done m2-done m3-done m4-done m5-done m6-done v0.1-pre-cluster
