#!/bin/bash
# run_playlists.sh -- run makeCrossSectionMCInputs on several Medium Energy (FHC)
# playlists in parallel, one detached screen session per playlist.
#
# Usage (interactive shell on the GPVM):
#   source setup_mat.sh                # once; the jobs inherit this environment
#   ./run_playlists.sh 1B 1C 1D        # full runs (systematics + truth)
#   TEST=1 ./run_playlists.sh 1A       # truncated smoke test (test_run = true)
#
# Playlists may be given as 1B or ME1B.  For each job:
#   screen session  lnhn_ME1B   (screen -r lnhn_ME1B to watch; Ctrl-a d to detach)
#   log             logs/ME1B_<UTC timestamp>[_TEST].log
#   output          XSecInputs_1110_ME1B_<UTC date>.root in this directory.
#                   A TEST run's file is renamed ..._TEST.root when it finishes,
#                   so a smoke test never overwrites a full run's file.
#
# Jobs are started one at a time: loadLibs.C force-recompiles the include
# libraries ("kf") on every launch, so two jobs must never compile at once.
# The next job starts only after the previous one's log shows its data loop
# has begun (all compilation done).  Already-running jobs keep the libraries
# they loaded.  If a job dies or stalls before that, no further jobs start.
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
WAIT_MAX="${WAIT_MAX:-2700}"
POLL="${POLL:-15}"

if [ $# -eq 0 ]; then
  sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi
if [ -z "${TOPDIR:-}" ] || ! command -v root >/dev/null 2>&1; then
  echo "ERROR: environment not set up -- run 'source setup_mat.sh' in this shell first." >&2
  exit 1
fi
if ! pgrep -u "$USER" -f renew_token.py >/dev/null 2>&1; then
  echo "WARNING: renew_token.py is not running; long jobs may lose their token." >&2
fi
if [ "$TEST" = "1" ]; then TESTFLAG=true; SUFFIX=_TEST; else TESTFLAG=false; SUFFIX=; fi
mkdir -p "$REPO/logs"

session_alive() { screen -ls | grep -q "[.]$1[[:space:]]"; }

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
  SESSION="lnhn_${PL}"
  if session_alive "$SESSION"; then
    echo "ERROR: $PL: screen session $SESSION is already running -- skipping" >&2
    continue
  fi
  STAMP=$(date -u +%Y-%m-%dT%H%M%SZ)
  LOG="$REPO/logs/${PL}_${STAMP}${SUFFIX}.log"
  JOB="$REPO/logs/${PL}_${STAMP}${SUFFIX}.job.sh"
  cat > "$JOB" <<JOBEOF
#!/bin/bash
cd "$REPO"
root -l -q -b loadLibs.C+ 'xsec/makeCrossSectionMCInputs.C+(1,"$PL",true,true,false,"$MC","$DATA",$TESTFLAG)' > "$LOG" 2>&1
rc=\$?
echo "=== root exit code \$rc, \$(date -u)" >> "$LOG"
if [ "$TESTFLAG" = "true" ]; then
  out=\$(sed -n 's/^Saving output to \(.*\.root\).*/\1/p' "$LOG" | tail -1)
  if [ -n "\$out" ] && [ -f "\$out" ]; then
    mv "\$out" "\${out%.root}_TEST.root" && echo "=== renamed \$out -> \${out%.root}_TEST.root" >> "$LOG"
  fi
fi
JOBEOF
  screen -dmS "$SESSION" bash "$JOB"
  echo "$PL: started (screen $SESSION, log logs/$(basename "$LOG"))"

  # Wait until this job is past compilation before starting the next one.
  waited=0
  while ! grep -q "Starting Data Loop" "$LOG" 2>/dev/null; do
    if ! session_alive "$SESSION"; then
      echo "ERROR: $PL ended before reaching its data loop -- see $LOG. Not starting any more jobs." >&2
      exit 1
    fi
    if [ "$waited" -ge "$WAIT_MAX" ]; then
      echo "ERROR: $PL has not reached its data loop after ${WAIT_MAX}s -- check $LOG. Not starting any more jobs." >&2
      exit 1
    fi
    sleep "$POLL"
    waited=$((waited + POLL))
  done
  echo "$PL: compiled; data loop started after ~${waited}s"
done
echo "Done launching.  Watch with:  screen -ls ;  tail -f logs/<playlist>_*.log"
