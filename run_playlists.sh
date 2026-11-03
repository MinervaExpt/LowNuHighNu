#!/bin/bash
# run_playlists.sh -- run makeCrossSectionMCInputs on several Medium Energy (FHC)
# playlists in parallel, one detached screen session per job.
#
# Usage (interactive shell on the GPVM):
#   source setup_mat.sh                       # once; the jobs inherit this environment
#   ./run_playlists.sh 1B 1C 1D               # full runs (systematics + truth)
#   TEST=1 ./run_playlists.sh 1A              # truncated smoke test (test_run = true)
#   TAG=wk45 ./run_playlists.sh 1D            # tagged run: coexists with untagged/other-tag runs
#   TAG=wk45 SPLIT=2 ./run_playlists.sh 1D    # split each playlist's file lists in 2 parallel jobs
#
# Playlists may be given as 1B or ME1B.  For each job (P = ME1B, plus _<TAG>,
# plus _p<i>of<N> when split, plus _TEST in test mode):
#   screen session  lnhn_<P>     (screen -r lnhn_<P> to watch; Ctrl-a d to detach)
#   log             logs/<PL>_<UTC timestamp>[_<TAG>][_p<i>of<N>][_TEST].log
#   output          XSecInputs_1110_<PL>_<UTC date>[_<TAG>][_p<i>of<N>][_TEST].root
#                   in this directory.  TEST runs write their _TEST name from the
#                   start, so a smoke test can never touch a full run's file.
#
# SPLIT=N (default 1) deals each playlist's MC and data file lists round-robin
# into N part lists (written to logs/) and runs one job per part.  POT is
# counted per part, so the merged file's POT is the sum.  When every part has
# finished, merge them with the command printed at the end:
#   ./merge_parts.sh logs/<PL>_<timestamp>[_<TAG>]_parts.txt
# Keep SPLIT small on the shared GPVMs (2 = two cores); the grid is the proper
# route for wide splits.
#
# Jobs are started one at a time: loadLibs.C force-recompiles the include
# libraries ("kf") on every launch, so two jobs must never compile at once.
# The next job starts only after the previous one's log shows its data loop
# has begun (all compilation done).  Already-running jobs keep the libraries
# they loaded.  If a job dies or stalls before that, no further jobs start.
#
# For several playlists, run this inside a screen (e.g. screen -S launcher):
# it waits between starts, and an ssh drop would kill the waiting loop.  For a
# single job it returns as soon as that job reaches its data loop.
#
# Overnight runs need the token renewer alive in its own screen:
#   screen -S renew_token ;  nohup python renew_token.py ;  Ctrl-a d
#
# Env overrides: WAIT_MAX (s, default 2700) = longest wait for a job to reach
# its data loop; OPENDATA_DIR = playlist-file root (default below).

set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OD="${OPENDATA_DIR:-/pnfs/minerva/persistent/OpenData/MediumEnergy_FHC}"
TEST="${TEST:-0}"
TAG="${TAG:-}"
SPLIT="${SPLIT:-1}"
WAIT_MAX="${WAIT_MAX:-2700}"
POLL="${POLL:-15}"

if [ $# -eq 0 ]; then
  sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi
if [ -n "$TAG" ] && ! [[ "$TAG" =~ ^[A-Za-z0-9-]+$ ]]; then
  echo "ERROR: TAG may contain only letters, digits and '-' (got '$TAG')." >&2
  exit 1
fi
if ! [[ "$SPLIT" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: SPLIT must be a positive integer (got '$SPLIT')." >&2
  exit 1
fi
if [ -z "${TOPDIR:-}" ] || ! command -v root >/dev/null 2>&1; then
  echo "ERROR: environment not set up -- run 'source setup_mat.sh' in this shell first." >&2
  exit 1
fi
if [ "$(cd "${TOPDIR%/}/LowNuHighNu" 2>/dev/null && pwd)" != "$REPO" ]; then
  echo "ERROR: TOPDIR=$TOPDIR does not contain this checkout ($REPO); loadLibs.C would" >&2
  echo "       compile another checkout's includes.  Source this checkout's setup_mat.sh." >&2
  exit 1
fi
if ! pgrep -u "$USER" -f renew_token.py >/dev/null 2>&1; then
  echo "WARNING: renew_token.py is not running; long jobs may lose their token." >&2
fi
if [ "$TEST" = "1" ]; then TESTFLAG=true; SUFFIX=_TEST; else TESTFLAG=false; SUFFIX=; fi
TAGSFX="${TAG:+_$TAG}"
mkdir -p "$REPO/logs"

session_alive() { screen -ls | grep -q "[.]$1[[:space:]]"; }

# launch_job <PL> <mc list> <data list> <out_tag> <name suffix> <stamp>
# Starts one job in a detached screen and waits until it reaches its data loop.
# Returns 2 if its screen session already exists (not started), 1 if it died or
# stalled before the data loop.
launch_job() {
  local PL="$1" MCL="$2" DATAL="$3" OUT_TAG="$4" NSFX="$5" STAMP="$6"
  local SESSION="lnhn_${PL}${NSFX}${SUFFIX}"
  [ "$TESTFLAG" = "true" ] && OUT_TAG="${OUT_TAG:+${OUT_TAG}_}TEST"
  local LOG="$REPO/logs/${PL}_${STAMP}${NSFX}${SUFFIX}.log"
  local JOB="$REPO/logs/${PL}_${STAMP}${NSFX}${SUFFIX}.job.sh"
  if session_alive "$SESSION"; then
    echo "ERROR: screen session $SESSION is already running -- not starting it (use a different TAG)" >&2
    return 2
  fi
  cat > "$JOB" <<JOBEOF
#!/bin/bash
cd "$REPO"
echo "=== start \$(date -u)  (TOPDIR=\$TOPDIR)" > "$LOG"
root -l -q -b loadLibs.C+ 'xsec/makeCrossSectionMCInputs.C+(1,"$PL",true,true,false,"$MCL","$DATAL",$TESTFLAG,"$OUT_TAG")' >> "$LOG" 2>&1
rc=\$?
echo "=== root exit code \$rc, \$(date -u)" >> "$LOG"
JOBEOF
  screen -dmS "$SESSION" bash "$JOB"
  echo "$PL$NSFX: started (screen $SESSION, log logs/$(basename "$LOG"))"
  [ -n "${PARTS_FILE:-}" ] && echo "$LOG" >> "$PARTS_FILE"

  # Wait until this job is past compilation before starting the next one.
  local waited=0
  while ! grep -q "Starting Data Loop" "$LOG" 2>/dev/null; do
    if ! session_alive "$SESSION"; then
      echo "ERROR: $PL$NSFX ended before reaching its data loop -- see $LOG. Not starting any more jobs." >&2
      return 1
    fi
    if [ "$waited" -ge "$WAIT_MAX" ]; then
      echo "ERROR: $PL$NSFX has not reached its data loop after ${WAIT_MAX}s -- check $LOG. Not starting any more jobs." >&2
      return 1
    fi
    sleep "$POLL"
    waited=$((waited + POLL))
  done
  echo "$PL$NSFX: compiled; data loop started after ~${waited}s"
}

MERGE_HINTS=()
for arg in "$@"; do
  PL="ME${arg#ME}"
  NUM="${PL#ME}"
  MC="$OD/MC/StandardMC/MediumEnergy_FHC_StandardMC_Playlist${NUM}.txt"
  DATA="$OD/Data/MediumEnergy_FHC_Data_Playlist${NUM}.txt"
  for f in "$MC" "$DATA"; do
    if [ ! -f "$f" ]; then
      echo "ERROR: $PL: missing $f -- skipping" >&2
      continue 2
    fi
  done
  STAMP=$(date -u +%Y-%m-%dT%H%M%SZ)
  if [ "$SPLIT" -eq 1 ]; then
    PARTS_FILE= launch_job "$PL" "$MC" "$DATA" "$TAG" "$TAGSFX" "$STAMP"
    rc=$?
    [ "$rc" -eq 2 ] && continue   # already running: skip this playlist
    [ "$rc" -ne 0 ] && exit 1
  else
    export PARTS_FILE="$REPO/logs/${PL}_${STAMP}${TAGSFX}${SUFFIX}_parts.txt"
    : > "$PARTS_FILE"
    for i in $(seq 1 "$SPLIT"); do
      PSFX="${TAGSFX}_p${i}of${SPLIT}"
      MCL="$REPO/logs/${PL}_${STAMP}${PSFX}_mc.txt"
      DATAL="$REPO/logs/${PL}_${STAMP}${PSFX}_data.txt"
      # Round-robin over non-empty lines, so each part gets a similar mix of runs.
      awk -v n="$SPLIT" -v i="$i" 'NF { k++; if ((k - 1) % n == i - 1) print }' "$MC" > "$MCL"
      awk -v n="$SPLIT" -v i="$i" 'NF { k++; if ((k - 1) % n == i - 1) print }' "$DATA" > "$DATAL"
      echo "$PL${PSFX}: $(wc -l < "$MCL") MC files, $(wc -l < "$DATAL") data files"
      launch_job "$PL" "$MCL" "$DATAL" "${TAG:+${TAG}_}p${i}of${SPLIT}" "$PSFX" "$STAMP" || exit 1
    done
    MERGE_HINTS+=("./merge_parts.sh logs/$(basename "$PARTS_FILE")")
    unset PARTS_FILE
  fi
done
echo "Done launching.  Watch with:  screen -ls ;  tail -f logs/<playlist>_*.log"
for h in ${MERGE_HINTS[@]+"${MERGE_HINTS[@]}"}; do echo "When all parts have finished:  $h"; done
