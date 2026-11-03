#!/bin/bash
# merge_parts.sh -- merge the part outputs of a SPLIT run of run_playlists.sh.
#
# Usage (GPVM, after `source setup_mat.sh`, once every part has finished):
#   ./merge_parts.sh logs/<PL>_<timestamp>[_<TAG>][_TEST]_parts.txt
#
# Reads the part logs listed in the parts file, checks that each ended with
# "=== root exit code 0", and merges their output files with MAT's `madd`
# (adds MnvH1D/MnvH2D universe by universe; the 1-bin mc_pot / data_pot hists
# therefore sum to the total POT; other objects such as git_describe are copied
# from the first part).  The merged file drops the _p<i>of<N> part of the name:
#   XSecInputs_1110_ME1D_<date>_wk45_p1of2.root + ..._p2of2.root
#     -> XSecInputs_1110_ME1D_<date>_wk45.root
# The part files are left in place.  Refuses to overwrite an existing target
# unless FORCE=1.

set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARTS="${1:-}"
if [ -z "$PARTS" ]; then sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 1; fi
[ -f "$PARTS" ] || PARTS="$REPO/$PARTS"
[ -f "$PARTS" ] || { echo "ERROR: parts file not found: $1" >&2; exit 1; }
cd "$REPO"

MADD="$(command -v madd 2>/dev/null || true)"
if [ -z "$MADD" ] && [ -x "${PLOTUTILSROOT:-/nonexistent}/../bin/madd" ]; then
  MADD="$PLOTUTILSROOT/../bin/madd"
fi
[ -n "$MADD" ] || { echo "ERROR: madd not found (installed with MAT in opt/bin) -- source setup_mat.sh first." >&2; exit 1; }

outs=()
while read -r log; do
  [ -z "$log" ] && continue
  if ! grep -q '^=== root exit code 0' "$log" 2>/dev/null; then
    echo "ERROR: $log has not finished cleanly (no '=== root exit code 0') -- not merging." >&2
    exit 1
  fi
  out=$(sed -n 's/^Saving output to \(.*\.root\).*/\1/p' "$log" | tail -1)
  if [ -z "$out" ] || [ ! -f "$out" ]; then
    echo "ERROR: output file for $log not found ('$out')." >&2
    exit 1
  fi
  outs+=("$out")
done < "$PARTS"
[ "${#outs[@]}" -ge 2 ] || { echo "ERROR: fewer than 2 part outputs listed in $PARTS." >&2; exit 1; }

target=$(echo "${outs[0]}" | sed 's/_p1of[0-9][0-9]*//')
if [ "$target" = "${outs[0]}" ]; then
  echo "ERROR: could not derive the merged name from ${outs[0]} (no _p1of<N>)." >&2
  exit 1
fi
if [ -e "$target" ] && [ "${FORCE:-0}" != "1" ]; then
  echo "ERROR: $target already exists (FORCE=1 to overwrite)." >&2
  exit 1
fi

echo "Merging ${#outs[@]} parts -> $target"
printf '  %s\n' "${outs[@]}"
"$MADD" "$target" "${outs[@]}" || { echo "ERROR: madd failed." >&2; exit 1; }

# Best-effort POT check: the merged POT should equal the sum of the parts.
python - "$target" "${outs[@]}" <<'PY' || echo "(POT check skipped)"
import sys
import ROOT
ROOT.gROOT.SetBatch(True)
def pot(path):
    f = ROOT.TFile.Open(path)
    vals = tuple(f.Get(n).GetBinContent(1) for n in ("mc_pot", "data_pot"))
    f.Close()
    return vals
target, parts = sys.argv[1], sys.argv[2:]
mc = sum(pot(p)[0] for p in parts); data = sum(pot(p)[1] for p in parts)
tmc, tdata = pot(target)
ok = abs(tmc - mc) <= 1e-6 * mc and abs(tdata - data) <= 1e-6 * max(data, 1.0)
print("POT  parts sum: mc=%.4e data=%.4e   merged: mc=%.4e data=%.4e   %s"
      % (mc, data, tmc, tdata, "OK" if ok else "MISMATCH"))
sys.exit(0 if ok else 1)
PY
echo "Done: $target"
