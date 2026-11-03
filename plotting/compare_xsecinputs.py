#!/usr/bin/env python
"""compare_xsecinputs.py -- regression check between two output files.

Usage:
    python compare_xsecinputs.py <reference.root> <candidate.root> [--rtol 1e-9] [--verbose]

For every histogram present in both files, compares the CV contents
(including under/overflow) and, for MnvH1D / MnvH2D, every universe of every
vertical and lateral error band (band names and universe counts must match).
Lists histograms present in only one file.  Non-histogram objects (e.g. the
git_describe / git_diff_md5 provenance stamps) are reported, not compared.

Subdirectories are compared too, so the script also checks two processedHists
files (processHists.py output).  Objects are named by their path, "ME1D/eff_lowNu_ME1D";
a name written more than once in one directory is compared cycle by cycle,
"ME1D/migration_enu;2".

Use it on two runs over the same events (e.g. two TEST=1 runs of the same
playlist built from different code) to show that a change is
behavior-preserving: every common histogram should agree to rounding.

Exit status: 0 if all common histograms agree within rtol, 1 otherwise.
"""
import argparse
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
try:
    from ROOT import PlotUtils  # noqa: F401  (MnvH1D/MnvH2D dictionaries)
except Exception as exc:  # pragma: no cover
    print("WARNING: could not import PlotUtils (%s); error bands will not be compared" % exc)

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None


def contents(h):
    """All cell contents (incl. under/overflow) as a list or numpy array."""
    n = h.GetNcells()
    if np is not None:
        try:
            buf = h.GetArray()
            buf.reshape((n,))
            return np.array(buf, dtype=np.float64, copy=True)
        except Exception:
            pass
    return [h.GetBinContent(i) for i in range(n)]


def max_rel_diff(a, b):
    if np is not None and isinstance(a, np.ndarray):
        if a.shape != b.shape:
            return float("inf")
        scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-300)
        d = np.abs(a - b) / scale
        d[(a == b)] = 0.0  # exact matches, incl. both zero (and both NaN-free)
        both_nan = np.isnan(a) & np.isnan(b)
        d[both_nan] = 0.0
        return float(np.nanmax(d)) if d.size else 0.0
    if len(a) != len(b):
        return float("inf")
    worst = 0.0
    for x, y in zip(a, b):
        if x == y:
            continue
        worst = max(worst, abs(x - y) / max(abs(x), abs(y), 1e-300))
    return worst


def bands(h, kind):
    getnames = getattr(h, "Get%sErrorBandNames" % kind, None)
    if getnames is None:
        return {}
    out = {}
    for name in getnames():
        band = getattr(h, "Get%sErrorBand" % kind)(str(name))
        out[str(name)] = [band.GetHist(i) for i in range(band.GetNHists())]
    return out


def compare(href, hcan, rtol):
    """Return a list of problem strings (empty if they agree)."""
    problems = []
    if href.GetNcells() != hcan.GetNcells():
        return ["binning differs (%d vs %d cells)" % (href.GetNcells(), hcan.GetNcells())]
    d = max_rel_diff(contents(href), contents(hcan))
    if d > rtol:
        problems.append("CV max rel diff %.3g" % d)
    for kind in ("Vert", "Lat"):
        bref, bcan = bands(href, kind), bands(hcan, kind)
        if set(bref) != set(bcan):
            only_r = sorted(set(bref) - set(bcan))
            only_c = sorted(set(bcan) - set(bref))
            problems.append("%s band names differ (ref only %s, cand only %s)" % (kind, only_r[:5], only_c[:5]))
        for name in sorted(set(bref) & set(bcan)):
            ur, uc = bref[name], bcan[name]
            if len(ur) != len(uc):
                problems.append("%s band %s: %d vs %d universes" % (kind, name, len(ur), len(uc)))
                continue
            worst = max([max_rel_diff(contents(a), contents(b)) for a, b in zip(ur, uc)] or [0.0])
            if worst > rtol:
                problems.append("%s band %s max rel diff %.3g" % (kind, name, worst))
    return problems


def keys(directory, prefix=""):
    """{path: class name} for every object under directory, descending into subdirectories."""
    listed = list(directory.GetListOfKeys())
    counts = {}
    for k in listed:
        counts[k.GetName()] = counts.get(k.GetName(), 0) + 1
    out = {}
    for k in listed:
        name = prefix + k.GetName()
        cls = ROOT.TClass.GetClass(k.GetClassName())
        if cls and cls.InheritsFrom("TDirectory"):
            out.update(keys(k.ReadObj(), name + "/"))
            continue
        if counts[k.GetName()] > 1:
            name += ";%d" % k.GetCycle()
        out[name] = k.GetClassName()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("candidate")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--verbose", action="store_true", help="also list histograms that agree")
    args = ap.parse_args()

    fr, fc = ROOT.TFile.Open(args.reference), ROOT.TFile.Open(args.candidate)
    if not fr or fr.IsZombie() or not fc or fc.IsZombie():
        print("ERROR: could not open both files")
        return 1
    kr, kc = keys(fr), keys(fc)
    common = sorted(set(kr) & set(kc))
    only_r, only_c = sorted(set(kr) - set(kc)), sorted(set(kc) - set(kr))

    n_ok = n_bad = 0
    for name in common:
        or_, oc = fr.Get(name), fc.Get(name)
        if not (isinstance(or_, ROOT.TH1) and isinstance(oc, ROOT.TH1)):
            title_r = or_.GetTitle() if hasattr(or_, "GetTitle") else ""
            title_c = oc.GetTitle() if hasattr(oc, "GetTitle") else ""
            print("  (not compared) %-40s %s | %s" % (name, title_r, title_c))
            continue
        problems = compare(or_, oc, args.rtol)
        if problems:
            n_bad += 1
            print("  DIFF  %-48s %s" % (name, "; ".join(problems[:4])))
        else:
            n_ok += 1
            if args.verbose:
                print("  same  %s" % name)
    for name in only_r:
        print("  only in reference: %s (%s)" % (name, kr[name]))
    for name in only_c:
        print("  only in candidate: %s (%s)" % (name, kc[name]))
    print("SUMMARY: %d common histograms agree, %d differ (rtol %g); %d only in reference, %d only in candidate"
          % (n_ok, n_bad, args.rtol, len(only_r), len(only_c)))
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
