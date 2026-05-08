"""
boundary_diagnostic.py

Plot the low-nu/high-nu staircase boundary diagnostic histograms produced by
the event loop:
  - delta_nu          (1D, reco)   nu - nu_cut(Enu)
  - delta_nu_true     (1D, truth)  same in true kinematics
  - enu_ehad_boundary (2D)         (Enu, nu) on a fine uniform grid

Reads the raw event-loop output (XSecInputs_*.root) directly, NOT the
post-processed processedHists_*.root, because the diagnostic is concerned
with all events near the cut, not the lowNu/highNu split that processHists
imposes.  The Var2D already exposes an `_inclusive` variant which is what
we want.

Usage:
    python boundary_diagnostic.py <path-to-XSecInputs.root> [<output-dir>]

Histogram name conventions (verified against includes/Histograms.cxx and
includes/Variable2D.h):
  1D: selection_{data,mc}_<varname>  ; effnum_<varname> ; effden_<varname>
  2D: selection_{data,mc}_<varname>_{inclusive,lowNu,highNu}
"""

import os
import sys

import ROOT

# Match the import style of other plotting/ scripts; setup_mat.sh puts
# py_classes on PYTHONPATH.
from functions import *
from plottingClasses import *

ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
ROOT.myPlotStyle()
ROOT.TH1.AddDirectory(False)

# Staircase corners in (Enu_GeV, nu_GeV).  Must mirror LowNuBoundary.h.
STAIRCASE_CORNERS = [
    (0.0,  0.3),
    (3.0,  0.3),
    (3.0,  0.5),
    (7.0,  0.5),
    (7.0,  1.0),
    (12.0, 1.0),
    (12.0, 2.0),
    (15.0, 2.0),
]


def parse_args():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    in_path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "./boundary_diagnostic_plots"
    return in_path, out_dir


def staircase_polyline():
    """TPolyLine tracing the staircase, for overlay on the 2D map."""
    n = len(STAIRCASE_CORNERS)
    line = ROOT.TPolyLine(n)
    for i, (x, y) in enumerate(STAIRCASE_CORNERS):
        line.SetPoint(i, x, y)
    line.SetLineColor(ROOT.kRed + 2)
    line.SetLineWidth(3)
    line.SetLineStyle(2)
    return line


def get_cv(hist):
    """Return the CV TH1/TH2 from an MnvH1D/MnvH2D, with stat errors."""
    return hist.GetCVHistoWithStatError()


def cut_marker_line(y_min, y_max, x=0.0):
    line = ROOT.TLine(x, y_min, x, y_max)
    line.SetLineColor(ROOT.kRed + 2)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    return line


def plot_delta_nu_reco(h_data, h_mc, out_path):
    """1D Δν (reco) with data + MC overlay and a vertical line at Δν = 0."""
    cv_data = get_cv(h_data)
    cv_mc = get_cv(h_mc)

    cv_data.SetMarkerStyle(20)
    cv_data.SetMarkerColor(ROOT.kBlack)
    cv_data.SetLineColor(ROOT.kBlack)
    cv_mc.SetLineColor(ROOT.kBlue + 2)
    cv_mc.SetLineWidth(2)

    canvas = ROOT.TCanvas("c_delta_nu_reco", "delta_nu (reco)", 800, 600)
    canvas.SetLogy()

    cv_mc.SetTitle("#nu - #nu_{cut}(E_{#nu}) (reco)")
    cv_mc.GetXaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    cv_mc.GetYaxis().SetTitle("Events / bin")
    cv_mc.Draw("HIST")
    cv_data.Draw("PE SAME")

    y_min = max(cv_mc.GetMinimum(), 0.1)
    y_max = max(cv_mc.GetMaximum(), cv_data.GetMaximum()) * 2.0
    line = cut_marker_line(y_min, y_max, 0.0)
    line.Draw("SAME")

    legend = ROOT.TLegend(0.65, 0.75, 0.9, 0.9)
    legend.AddEntry(cv_data, "Data", "PE")
    legend.AddEntry(cv_mc, "MC reco", "L")
    legend.AddEntry(line, "#Delta#nu = 0 (cut)", "L")
    legend.Draw()

    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_delta_nu_truth(h_mc_true, out_path):
    """1D Δν (truth) -- MC only; data has no truth."""
    cv = get_cv(h_mc_true)
    cv.SetLineColor(ROOT.kBlue + 2)
    cv.SetLineWidth(2)

    canvas = ROOT.TCanvas("c_delta_nu_truth", "delta_nu (truth)", 800, 600)
    canvas.SetLogy()

    cv.SetTitle("#nu - #nu_{cut}(E_{#nu}) (truth, MC)")
    cv.GetXaxis().SetTitle("#Delta#nu_{true}  [GeV]")
    cv.GetYaxis().SetTitle("Events / bin")
    cv.Draw("HIST")

    y_min = max(cv.GetMinimum(), 0.1)
    y_max = cv.GetMaximum() * 2.0
    line = cut_marker_line(y_min, y_max, 0.0)
    line.Draw("SAME")

    legend = ROOT.TLegend(0.65, 0.78, 0.9, 0.88)
    legend.AddEntry(cv, "MC truth", "L")
    legend.AddEntry(line, "#Delta#nu = 0 (cut)", "L")
    legend.Draw()

    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_enu_ehad_boundary(h_2d, out_path, title_extra=""):
    """2D (E_nu, nu) heatmap with the staircase overlaid as a red dashed
    polyline."""
    cv = get_cv(h_2d)

    canvas = ROOT.TCanvas("c_enu_ehad_boundary", "enu_ehad_boundary", 900, 700)
    canvas.SetLogz()
    canvas.SetRightMargin(0.13)

    cv.SetTitle("E_{{#nu}} vs #nu near the staircase cut{0}".format(title_extra))
    cv.GetXaxis().SetTitle("E_{#nu}  [GeV]")
    cv.GetYaxis().SetTitle("#nu = E_{had}  [GeV]")
    cv.Draw("COLZ")

    line = staircase_polyline()
    line.Draw("SAME")

    legend = ROOT.TLegend(0.55, 0.80, 0.84, 0.88)
    legend.AddEntry(line, "Low-#nu cut boundary", "L")
    legend.Draw()

    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def safe_get(f, name, fallback_substring=None):
    """TFile.Get, returning None on miss.

    On miss, if `fallback_substring` is provided, list all top-level keys in
    the file whose names contain that substring -- helps diagnose naming
    drift (e.g. Var2D's auto-joined names) without requiring a separate
    file inspection.
    """
    h = f.Get(name)
    if h:
        return h
    print("  WARN: '{0}' not found.".format(name))
    if fallback_substring is not None:
        matches = [k.GetName() for k in f.GetListOfKeys()
                   if fallback_substring in k.GetName()]
        if matches:
            print("        Keys in file containing '{0}':".format(fallback_substring))
            for m in sorted(matches):
                print("          {0}".format(m))
        else:
            print("        (no keys containing '{0}' found either)".format(fallback_substring))
    return None


def main():
    in_path, out_dir = parse_args()

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        print("Created output directory: {0}".format(out_dir))

    f = ROOT.TFile.Open(in_path, "READ")
    if not f or f.IsZombie():
        sys.stderr.write("ERROR: Could not open input file: {0}\n".format(in_path))
        sys.exit(1)
    print("Reading: {0}".format(in_path))

    # ---- 1D delta_nu (reco) ----
    h_data = safe_get(f, "selection_data_delta_nu")
    h_mc = safe_get(f, "selection_mc_delta_nu")
    if h_data and h_mc:
        plot_delta_nu_reco(h_data, h_mc, "{0}/delta_nu_reco.png".format(out_dir))

    # ---- 1D delta_nu (truth) ----
    h_mc_true = safe_get(f, "selection_mc_delta_nu_true")
    if h_mc_true:
        plot_delta_nu_truth(h_mc_true, "{0}/delta_nu_truth.png".format(out_dir))

    # ---- 2D enu_ehad_boundary (inclusive: all events passing cuts) ----
    # Var2D auto-joins its two VarMAT names. With VarMATs named
    # "enu_boundary" and "ehad_boundary" this produces the (ugly) joined
    # name "enu_boundary_ehad_boundary".  Look there first, then fall back
    # to the cleaner "enu_ehad_boundary" in case the event-loop code is
    # later cleaned up to produce the better name.
    candidate_names_data = [
        "selection_data_enu_boundary_ehad_boundary_inclusive",
        "selection_data_enu_ehad_boundary_inclusive",
    ]
    candidate_names_mc = [
        "selection_mc_enu_boundary_ehad_boundary_inclusive",
        "selection_mc_enu_ehad_boundary_inclusive",
    ]
    h_2d_data = None
    for name in candidate_names_data:
        h_2d_data = f.Get(name)
        if h_2d_data:
            print("  using {0}".format(name))
            break
    if not h_2d_data:
        safe_get(f, candidate_names_data[0], fallback_substring="boundary")
    else:
        plot_enu_ehad_boundary(
            h_2d_data,
            "{0}/enu_ehad_boundary_data.png".format(out_dir),
            title_extra=" -- data",
        )

    h_2d_mc = None
    for name in candidate_names_mc:
        h_2d_mc = f.Get(name)
        if h_2d_mc:
            print("  using {0}".format(name))
            break
    if not h_2d_mc:
        safe_get(f, candidate_names_mc[0], fallback_substring="boundary")
    else:
        plot_enu_ehad_boundary(
            h_2d_mc,
            "{0}/enu_ehad_boundary_mc.png".format(out_dir),
            title_extra=" -- MC reco",
        )

    print("Done.  Plots in {0}".format(out_dir))


if __name__ == "__main__":
    main()
