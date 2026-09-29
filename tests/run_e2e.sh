#!/usr/bin/env bash
#
# End-to-end test against the real somalier binary.
#
#   bash tests/run_e2e.sh [workdir]
#
# Builds a synthetic cohort with a known identity structure, runs the whole
# pipeline on it, and asserts that the verdicts come back as designed. Needs
# Docker (or swap -profile for singularity/conda).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
OUT="${1:-${TMPDIR:-/tmp}/concordance-e2e}"
PROFILE="${PROFILE:-docker}"

rm -rf "$OUT"
mkdir -p "$OUT"

echo "==> running pipeline (-profile test,$PROFILE)"
# uses the committed fixture in assets/test_data; regenerate it with
#   python3 tests/make_fixture_vcfs.py assets/test_data
NXF_WORK="$OUT/work" nextflow run "$ROOT" -profile "test,$PROFILE" \
    --outdir "$OUT/results" \
    -ansi-log false

echo "==> asserting verdicts"
python3 - "$OUT/results/concordance/concordance.pairs.tsv" <<'PY'
import sys

# ground truth built into make_fixture_vcfs.py
EXPECTED = {
    ("s_a", "s_b"): "CONFIRMED_MATCH",       # same subject, same genotypes
    ("s_c", "s_e"): "UNEXPECTED_DUPLICATE",  # different subjects, same genotypes
    ("s_f", "s_g"): "SAMPLE_SWAP",           # same subject, different genotypes
}

rows = {}
with open(sys.argv[1]) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(header, line.rstrip("\n").split("\t")))
        rows[(r["sample_a"], r["sample_b"])] = r

failures = []
for key, want in EXPECTED.items():
    got = rows.get(key, {}).get("verdict")
    print("%s %-11s expected %-22s got %s"
          % ("ok " if got == want else "FAIL", "/".join(key), want, got))
    if got != want:
        failures.append("%s: expected %s, got %s" % ("/".join(key), want, got))

# every other pair is a genuinely independent draw and must come back distinct
for key, r in rows.items():
    if key in EXPECTED:
        continue
    if r["verdict"] != "CONFIRMED_DISTINCT":
        failures.append("%s: expected CONFIRMED_DISTINCT, got %s (relatedness %s)"
                        % ("/".join(key), r["verdict"], r["relatedness"]))
print("ok  %d independent pairs all CONFIRMED_DISTINCT"
      % sum(1 for k in rows if k not in EXPECTED))

if failures:
    print("\nFAILED:")
    for f in failures:
        print("  - " + f)
    sys.exit(1)
print("\ne2e passed")
PY
