# LowNuHighNu

MINERvA Medium Energy (FHC) charged-current inclusive analysis using the
**low-ν method** for an in-situ flux measurement, plus the high-ν method
(R. Fine, PhD thesis, University of Rochester, 2020). Built on the MINERvA
Analysis Toolkit: [MAT](https://github.com/MinervaExpt/MAT) and
[MAT-MINERvA](https://github.com/MinervaExpt/MAT-MINERvA). The input is the
MINERvA open-data ME FHC tuples. Working branch: `rob-dev`.

## Pipeline

```
open-data tuples ──► xsec/makeCrossSectionMCInputs.C ──► XSecInputs_1110_<PL>_<date>.root (one per playlist)
                     (compiled by loadLibs.C; run via run_playlists.sh)
XSecInputs ──► plotting/processHists.py ──► processedHists_<date>.root ──► plotting/new_makeAllPlots.py
XSecInputs ──► plotting/boundary_diagnostic.py   (low-ν staircase boundary diagnostics)
processedHists + thesis-era files ──► plotting/study_thesisComparison.py
```

## Running (FNAL GPVM)

```bash
cd /exp/minerva/app/users/$USER/MATAna/LowNuHighNu
source setup_mat.sh                      # ROOT, MAT/MAT-MINERvA (../opt), token
screen -S renew_token                    # long runs: keep the token alive
nohup python renew_token.py              #   then Ctrl-a d
TEST=1 ./run_playlists.sh 1D             # smoke test: 1,000 entries per loop, writes ..._TEST.root
./run_playlists.sh 1D                    # full run, one detached screen per job (lnhn_ME1D)
TAG=wk45 SPLIT=2 ./run_playlists.sh 1D   # tagged run split in two parallel jobs ...
./merge_parts.sh <run dir>/logs/ME1D_parts.txt   # ... merged with MAT's madd when both finish
python plotting/boundary_diagnostic.py <run dir>/XSecInputs_1110_ME1D_<date>.root [outdir]
```

Each launch writes into its own run directory on the data area,
`/exp/minerva/data/users/$USER/MATAna/runs/<date>_g<commit>[_<TAG>][_TEST]/`
(outputs, `logs/`, and a `RUN.txt` recording what was launched); `RUNS_ROOT=`
puts run directories elsewhere. `run_playlists.sh` with no arguments prints its full usage. Every output file
records the code version it was built from (`git_describe`, `git_diff_md5`).

## Layout

| Path | Contents |
|---|---|
| `includes/` | `CVUniverse` (variables, CV weight), `SignalDefinition.h` (signal, `IsLowNu`), `LowNuBoundary.h` (the ν₀(E_ν) staircase), `Cuts`, `Systematics.h` (error bands), `Binning.h`, the `Variable` / `Variable2D` hist containers, `LowNuHighNuEvent` (fills) |
| `xsec/makeCrossSectionMCInputs.C` | the event loop (data, MC reco, MC truth) |
| `plotting/`, `py_classes/` | post-processing and plots (Python + PyROOT) |
| `run_playlists.sh`, `merge_parts.sh`, `renew_token.py`, `setup_mat.sh`, `loadLibs.C` | operations |
| `studies/`, the rest of `xsec/`, `tuple-utils/`, `runTransWarp.sh`, CCPi grid scripts | inherited from the CC-CH-pip-ana fork; **held as templates** for later studies, not maintained and not expected to compile |

The code was forked from the CC-CH-pip-ana (single charged pion) analysis;
pion-specific leftovers are being removed as the analysis is cleaned up.
