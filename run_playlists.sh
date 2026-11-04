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
#   TUPLES=p9 TAG=wk45-p9 ./run_playlists.sh 1D   # same, on the P9 data-preservation tuples
#   TUNE=legacy TAG=legacy ./run_playlists.sh 1D  # CV tune: v4.4.1 (default) or legacy
#
# Playlists may be given as 1B or ME1B.
#
# RUN DIRECTORY.  Each launch writes everything into one new directory on the
# data area, $RUNS_ROOT/<ID>/, with
#   ID = <UTC date>_g<commit>[-dirty][_<TAG>][_TEST]   e.g. 2026-11-04_g8e1c6f2_wk45
# (-2, -3, ... appended if that directory already exists).  Inside it:
#   RUN.txt         what was launched: code version, tune, tuples, host, command;
#                   each job appends its start and its exit code
#   XSecInputs_1110_<PL>_<UTC date>[_<TAG>][_TEST].root   one output per playlist
#   parts/          split runs only: the part outputs (..._p<i>of<N>...root)
#   logs/           <PL>[_p<i>of<N>].log, the .job.sh with the exact call, and for
#                   split runs the part lists and <PL>_parts.txt
# Screen sessions are lnhn_<PL>[_<TAG>][_p<i>of<N>][_TEST] (screen -r <name> to
# watch; Ctrl-a d to detach).  TEST runs get their own run directory and _TEST
# names, so a smoke test can never touch a full run's files.  Jobs still run in
# this checkout (loadLibs.C compiles its includes); only the outputs go to the
# run directory.
#
# SPLIT=N (default 1) deals each playlist's MC and data file lists round-robin
# into N part lists and runs one job per part.  POT is counted per part, so the
# merged file's POT is the sum.  When every part has finished, merge them with
# the command printed at the end:
#   ./merge_parts.sh <run dir>/logs/<PL>_parts.txt
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
# Input tuples (TUPLES, default opendata) -- the playlist files are, for <N> = 1A ... 1P:
#   opendata  $OPENDATA_DIR/MC/StandardMC/MediumEnergy_FHC_StandardMC_Playlist<N>.txt
#             $OPENDATA_DIR/Data/MediumEnergy_FHC_Data_Playlist<N>.txt
#             (OPENDATA_DIR default /pnfs/minerva/persistent/OpenData/MediumEnergy_FHC)
#   p9        $P9_DIR/FullDetector/P9_Playlist<N>_MC.txt
#             $P9_DIR/Data/P9_Playlist<N>_Data.txt
#             (P9_DIR default /pnfs/minerva/persistent/DataPreservation/p9)
# The playlist name passed to the event loop is ME<N> either way (MAT uses it
# for playlist-dependent weights); use TAG to keep the outputs apart.
#
# Env overrides: RUNS_ROOT (default /exp/minerva/data/users/$USER/MATAna/runs)
# = where run directories are made; WAIT_MAX (s, default 2700) = longest wait
# for a job to reach its data loop; OPENDATA_DIR, P9_DIR = tuple-list roots
# (defaults above).

set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OD="${OPENDATA_DIR:-/pnfs/minerva/persistent/OpenData/MediumEnergy_FHC}"
P9D="${P9_DIR:-/pnfs/minerva/persistent/DataPreservation/p9}"
TUPLES="${TUPLES:-opendata}"
TUNE="${TUNE:-v4.4.1}"
TEST="${TEST:-0}"
TAG="${TAG:-}"
SPLIT="${SPLIT:-1}"
WAIT_MAX="${WAIT_MAX:-2700}"
RUNS_ROOT="${RUNS_ROOT:-/exp/minerva/data/users/$USER/MATAna/runs}"
POLL="${POLL:-15}"

if [ $# -eq 0 ]; then
  sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi
case "$TUNE" in
  v4.4.1|legacy) ;;
  *) echo "ERROR: TUNE must be 'v4.4.1' or 'legacy' (got '$TUNE')." >&2; exit 1 ;;
esac
case "$TUPLES" in
  opendata|p9) ;;
  *) echo "ERROR: TUPLES must be 'opendata' or 'p9' (got '$TUPLES')." >&2; exit 1 ;;
esac
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

# The run directory: one per launch, named for the date, commit and tag.
GIT="git -C $REPO --no-optional-locks"
COMMIT=$($GIT rev-parse --short=7 HEAD) || { echo "ERROR: $REPO is not a git checkout." >&2; exit 1; }
$GIT diff --quiet HEAD || COMMIT="${COMMIT}-dirty"
RUN_ID="$(date -u +%Y-%m-%d)_g${COMMIT}${TAGSFX}${SUFFIX}"
if ! mkdir -p "$RUNS_ROOT" 2>/dev/null || [ ! -w "$RUNS_ROOT" ]; then
  echo "ERROR: cannot write to RUNS_ROOT=$RUNS_ROOT" >&2
  exit 1
fi
RUN_DIR="$RUNS_ROOT/$RUN_ID"
n=1
while [ -e "$RUN_DIR" ]; do n=$((n + 1)); RUN_DIR="$RUNS_ROOT/${RUN_ID}-$n"; done
RUN_ID="$(basename "$RUN_DIR")"
mkdir -p "$RUN_DIR/logs" || { echo "ERROR: cannot create $RUN_DIR" >&2; exit 1; }
{
  echo "run_id      $RUN_ID"
  echo "created     $(date -u +%Y-%m-%dT%H:%M:%SZ)  host $(hostname -s)  user $USER"
  echo "checkout    $REPO"
  echo "git         $($GIT describe --always --dirty)  (diff md5 $($GIT diff HEAD | md5sum | cut -c1-12))"
  echo "tune        $TUNE"
  echo "tuples      $TUPLES"
  echo "tag         ${TAG:-<none>}"
  echo "split       $SPLIT"
  echo "test        $TESTFLAG"
  echo "playlists   $*"
  echo "command     TUNE=$TUNE TUPLES=$TUPLES TAG=$TAG SPLIT=$SPLIT TEST=$TEST $0 $*"
  echo "--- jobs (start / exit lines are appended by each job)"
} > "$RUN_DIR/RUN.txt"
echo "Run directory: $RUN_DIR"
df -h "$RUNS_ROOT" | tail -n 1 | awk '{print "Free space there: " $4 " of " $2}'

session_alive() { screen -ls | grep -q "[.]$1[[:space:]]"; }

# launch_job <PL> <mc list> <data list> <out_tag> <part suffix> <output dir>
# Starts one job in a detached screen and waits until it reaches its data loop.
# Returns 2 if its screen session already exists (not started), 1 if it died or
# stalled before the data loop.
launch_job() {
  local PL="$1" MCL="$2" DATAL="$3" OUT_TAG="$4" PSFX="$5" OUTDIR="$6"
  local SESSION="lnhn_${PL}${TAGSFX}${PSFX}${SUFFIX}"
  [ "$TESTFLAG" = "true" ] && OUT_TAG="${OUT_TAG:+${OUT_TAG}_}TEST"
  local LOG="$RUN_DIR/logs/${PL}${PSFX}.log"
  local JOB="$RUN_DIR/logs/${PL}${PSFX}.job.sh"
  if session_alive "$SESSION"; then
    echo "ERROR: screen session $SESSION is already running -- not starting it (use a different TAG)" >&2
    return 2
  fi
  cat > "$JOB" <<JOBEOF
#!/bin/bash
cd "$REPO"
echo "=== start \$(date -u)  (TOPDIR=\$TOPDIR, TUPLES=$TUPLES, TUNE=$TUNE, RUN_DIR=$RUN_DIR)" > "$LOG"
echo "=== MC list: $MCL" >> "$LOG"
echo "=== data list: $DATAL" >> "$LOG"
echo "job \$(printf '%-16s' "$PL$PSFX") start \$(date -u +%Y-%m-%dT%H:%M:%SZ)  host \$(hostname -s)  screen $SESSION" >> "$RUN_DIR/RUN.txt"
root -l -q -b loadLibs.C+ 'xsec/makeCrossSectionMCInputs.C+(1,"$PL",true,true,false,"$MCL","$DATAL",$TESTFLAG,"$OUT_TAG","$TUNE","$OUTDIR")' >> "$LOG" 2>&1
rc=\$?
echo "=== root exit code \$rc, \$(date -u)" >> "$LOG"
echo "job \$(printf '%-16s' "$PL$PSFX") exit  \$(date -u +%Y-%m-%dT%H:%M:%SZ)  code \$rc" >> "$RUN_DIR/RUN.txt"
JOBEOF
  screen -dmS "$SESSION" bash "$JOB"
  echo "$PL$PSFX: started (screen $SESSION, log $LOG)"
  [ -n "${PARTS_FILE:-}" ] && echo "$LOG" >> "$PARTS_FILE"

  # Wait until this job is past compilation before starting the next one.
  local waited=0
  while ! grep -q "Starting Data Loop" "$LOG" 2>/dev/null; do
    if ! session_alive "$SESSION"; then
      echo "ERROR: $PL$PSFX ended before reaching its data loop -- see $LOG. Not starting any more jobs." >&2
      return 1
    fi
    if [ "$waited" -ge "$WAIT_MAX" ]; then
      echo "ERROR: $PL$PSFX has not reached its data loop after ${WAIT_MAX}s -- check $LOG. Not starting any more jobs." >&2
      return 1
    fi
    sleep "$POLL"
    waited=$((waited + POLL))
  done
  echo "$PL$PSFX: compiled; data loop started after ~${waited}s"
}

MERGE_HINTS=()
for arg in "$@"; do
  PL="ME${arg#ME}"
  NUM="${PL#ME}"
  if [ "$TUPLES" = "p9" ]; then
    MC="$P9D/FullDetector/P9_Playlist${NUM}_MC.txt"
    DATA="$P9D/Data/P9_Playlist${NUM}_Data.txt"
  else
    MC="$OD/MC/StandardMC/MediumEnergy_FHC_StandardMC_Playlist${NUM}.txt"
    DATA="$OD/Data/MediumEnergy_FHC_Data_Playlist${NUM}.txt"
  fi
  for f in "$MC" "$DATA"; do
    if [ ! -f "$f" ]; then
      echo "ERROR: $PL: missing $f -- skipping" >&2
      continue 2
    fi
  done
  if [ "$SPLIT" -eq 1 ]; then
    PARTS_FILE= launch_job "$PL" "$MC" "$DATA" "$TAG" "" "$RUN_DIR"
    rc=$?
    [ "$rc" -eq 2 ] && continue   # already running: skip this playlist
    [ "$rc" -ne 0 ] && exit 1
  else
    mkdir -p "$RUN_DIR/parts"
    export PARTS_FILE="$RUN_DIR/logs/${PL}_parts.txt"
    : > "$PARTS_FILE"
    for i in $(seq 1 "$SPLIT"); do
      PSFX="_p${i}of${SPLIT}"
      MCL="$RUN_DIR/logs/${PL}${PSFX}_mc.txt"
      DATAL="$RUN_DIR/logs/${PL}${PSFX}_data.txt"
      # Round-robin over non-empty lines, so each part gets a similar mix of runs.
      awk -v n="$SPLIT" -v i="$i" 'NF { k++; if ((k - 1) % n == i - 1) print }' "$MC" > "$MCL"
      awk -v n="$SPLIT" -v i="$i" 'NF { k++; if ((k - 1) % n == i - 1) print }' "$DATA" > "$DATAL"
      echo "$PL${PSFX}: $(wc -l < "$MCL") MC files, $(wc -l < "$DATAL") data files"
      launch_job "$PL" "$MCL" "$DATAL" "${TAG:+${TAG}_}p${i}of${SPLIT}" "$PSFX" "$RUN_DIR/parts" || exit 1
    done
    MERGE_HINTS+=("./merge_parts.sh $PARTS_FILE")
    unset PARTS_FILE
  fi
done
echo "Done launching.  Run directory: $RUN_DIR"
echo "Watch with:  screen -ls ;  tail -f $RUN_DIR/logs/<playlist>*.log ;  cat $RUN_DIR/RUN.txt"
for h in ${MERGE_HINTS[@]+"${MERGE_HINTS[@]}"}; do echo "When all parts have finished:  $h"; done
echo "Register the run in RUNS.md as $RUN_ID."
