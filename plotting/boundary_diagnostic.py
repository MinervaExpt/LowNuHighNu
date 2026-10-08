"""
boundary_diagnostic.py   (v3.0 -- 2026-08-14)

Low-nu/high-nu staircase boundary diagnostics from the raw event-loop output
(XSecInputs_*.root, NOT processedHists_*.root -- the diagnostic spans the
lowNu/highNu boundary rather than splitting on it).

Produces, per the advisor action items:

  1D delta_nu (reco):
    - POT-normalized data/MC overlay with the MC systematic error band
    - data/MC ratio, with the MC systematic band drawn around 1
    - a grouped systematic error-summary breakdown (uses errorMaps.error_bands)
  2D:
    - the original (E_nu, nu) heatmaps with the staircase overlay
    - a reconstructed  Delta-nu vs E_nu  map: fine nu binning, nominal-analysis
      E_nu binning, built by shifting the (E_nu, E_had) boundary histogram by
      -nu_cut(E_nu) column-by-column.  The MC reconstruction now propagates the
      full set of vertical systematic universes (see below), so it is a real
      MnvH2D with error bands, not a CV-only TH2D.
    - a COLZ map of the MC reconstruction's total fractional systematic
      uncertainty over (E_nu, Delta-nu)
    - Delta-nu distributions in the four staircase E_nu regions (slices), each
      with the MC systematic band drawn on the overlay and a companion grouped
      systematic error-summary

  ==================================================================
  STATE (2026-08-14): both blockers noted in v2.1 are resolved in the event
  loop and this script is updated to consume the fixed output.

  (1) Native 1D reco DATA delta_nu is now filled.  The data loop in
      makeCrossSectionMCInputs.C was repointed onto the full `variables`
      vector (the dead truth-free `variables_lessTruth` copy that was looped
      but never written is gone), so selection_data_delta_nu now carries real
      data.  This script uses that native histogram whenever it is non-empty;
      the old 2D-projection reconstruction of the pooled 1D data (and the
      matched 2D-projection MC) is retained ONLY as a legacy fallback for
      pre-re-run files and auto-disables the moment native data is present.

  (2) The 2D boundary map's systematic universes are now synced.  A
      SyncAllHists(Variable2D&) overload was added and is called in the write
      loop, so the MC (E_nu, E_had) boundary MnvH2D now carries its vertical
      error bands.  Everything derived from the 2D MC -- the reconstructed
      Delta-nu vs E_nu map, its fractional-uncertainty map, and all four E_nu
      slices -- therefore now carries MC systematics.

  DATA IS STILL STAT-ONLY ON THE 2D, BY DESIGN.  SyncAllHists(Variable2D&)
  syncs only the MC wrappers; the data 2D wrappers are filled straight onto
  their CV hist and carry statistical uncertainty only, exactly as the 1D data
  does.  So on the E_nu slices the DATA points have statistical error bars and
  the SYSTEMATIC band shown is the MC band.  This is intentional, not a gap.

  VALIDATION PENDING: written against the MnvH2D band structure ahead of the
  full ME1A re-run.  Confirm against the first real synced XSecInputs file that
  (a) the 2D MC boundary map reports GetVertErrorBandNames() non-empty and
  (b) the slice error summaries populate.  Until then treat the banded outputs
  as structurally-correct-but-unvalidated.
  ==================================================================

Usage:
    python boundary_diagnostic.py <path-to-XSecInputs.root> [<output-dir>]

Histogram name conventions (verified against Histograms.cxx / Variable2D.h):
  1D: selection_{data,mc}_<varname>   (delta_nu, delta_nu_true)
  2D: selection_{data,mc}_enu_boundary_ehad_boundary_inclusive
"""

import os
import sys
from array import array

import ROOT
from functions import *          # brings in PlotUtils, array, math, os
from plottingClasses import *    # makeEnv_TCanvas, setPlotSpecs_*, etc.
from errorMaps import *          # brings in error_bands (OrderedDict)

ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
try:
    ROOT.myPlotStyle()
except Exception:
    pass
ROOT.TH1.AddDirectory(False)

# ---- MnvPlotter with the standard grouped error summary (as new_makeAllPlots) ----
plotter = PlotUtils.MnvPlotter()
plotter.error_summary_group_map.clear()
for group in error_bands:
    for error in error_bands[group]:
        plotter.error_summary_group_map[group].push_back(error)
plotter.SetLegendNColumns(2)

# Instantiating MnvPlotter overrides the global colour palette.  Restore the
# ROOT default (kBird) -- Rob's default 2D colour scale, matching the original
# boundary plots -- and re-assert it before each 2D draw below.
DEFAULT_PALETTE = ROOT.kBird
ROOT.gStyle.SetPalette(DEFAULT_PALETTE)

# Keep Python-side references to ROOT drawables so they aren't garbage-collected
# out from under the canvas before SaveAs.
_keep = []

# Separate keepalive for the remapped systematic-universe TH2Ds that back the
# reconstructed MnvH2D's error bands.  Some PlotUtils builds store (rather than
# deep-copy) the universe pointers handed to AddVertErrorBand, so these must
# outlive the MnvH2D that references them -- never let them be GC'd.
_univ_keep = []

# Nominal analysis E_nu binning (GeV), clipped to the boundary hist's 15 GeV top.
NOMINAL_ENU_EDGES = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 14., 15.]

# Staircase corners (E_nu, nu) mirroring LowNuBoundary.h
STAIRCASE_CORNERS = [
    (0.0, 0.3), (3.0, 0.3), (3.0, 0.5), (7.0, 0.5),
    (7.0, 1.0), (12.0, 1.0), (12.0, 2.0), (15.0, 2.0),
]
ENU_CORNERS = [3.0, 7.0, 12.0]   # E_nu values where nu_cut steps

# (name, E_nu_lo, E_nu_hi, nu_cut) for the four slices
STAIRCASE_REGIONS = [
    ("Enu_lt3",    0.0,  3.0,  0.3),
    ("Enu_3to7",   3.0,  7.0,  0.5),
    ("Enu_7to12",  7.0,  12.0, 1.0),
    ("Enu_12to15", 12.0, 15.0, 2.0),
]


def nu_cut(enu):
    """Low-nu staircase cut value nu_cut(E_nu) in GeV.  Mirrors LowNuBoundary.h."""
    if enu < 3.0:
        return 0.3
    elif enu < 7.0:
        return 0.5
    elif enu < 12.0:
        return 1.0
    else:
        return 2.0


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def parse_args():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    in_path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "./boundary_diagnostic_plots"
    return in_path, out_dir


def read_pot(f):
    """(mc_pot, data_pot) from the WritePOT single-bin histos; safe fallbacks."""
    mc, dat = 1.0, 1.0
    h = f.Get("mc_pot")
    if h:
        mc = h.GetBinContent(1)
    h = f.Get("data_pot")
    if h:
        dat = h.GetBinContent(1)
    if mc <= 0:
        print("  WARN: mc_pot <= 0; using mcScale = 1.0")
        mc = 1.0
    print("  POT: mc={0:.4g}  data={1:.4g}  ->  mcScale = data/mc = {2:.4g}"
          .format(mc, dat, dat / mc))
    return mc, dat


def get_cv(hist):
    """CV TH1/TH2 (with stat errors) from an MnvH1D/MnvH2D."""
    return hist.GetCVHistoWithStatError()


def _apply_frac_band(band_hist, cv_hist, frac_hist):
    """Set band_hist bin errors = cv content * fractional error, per bin."""
    if frac_hist is None:
        return
    for b in range(0, band_hist.GetNbinsX() + 2):
        band_hist.SetBinError(b, cv_hist.GetBinContent(b) * frac_hist.GetBinContent(b))


def safe_get(f, name, fallback_substring=None):
    """TFile.Get returning None on miss; optionally list keys matching a substring."""
    h = f.Get(name)
    if h:
        return h
    print("  WARN: '{0}' not found.".format(name))
    if fallback_substring is not None:
        matches = [k.GetName() for k in f.GetListOfKeys()
                   if fallback_substring in k.GetName()]
        if matches:
            print("        Keys containing '{0}':".format(fallback_substring))
            for m in sorted(matches):
                print("          {0}".format(m))
    return None


def get_2d(f, which):
    """Fetch the 2D boundary MnvH2D for 'data' or 'mc', trying both name styles."""
    candidates = [
        "selection_{0}_enu_boundary_ehad_boundary_inclusive".format(which),
        "selection_{0}_enu_ehad_boundary_inclusive".format(which),
    ]
    for name in candidates:
        h = f.Get(name)
        if h:
            print("  using {0}".format(name))
            return h
    safe_get(f, candidates[0], fallback_substring="boundary")
    return None


def staircase_polyline():
    n = len(STAIRCASE_CORNERS)
    line = ROOT.TPolyLine(n)
    for i, (x, y) in enumerate(STAIRCASE_CORNERS):
        line.SetPoint(i, x, y)
    line.SetLineColor(ROOT.kRed + 2)
    line.SetLineWidth(3)
    line.SetLineStyle(2)
    return line


def draw_vline(x=0.0, color=None):
    """Vertical dashed line spanning the current pad (handles log-y)."""
    if color is None:
        color = ROOT.kRed + 2
    ROOT.gPad.Update()
    uy1, uy2 = ROOT.gPad.GetUymin(), ROOT.gPad.GetUymax()
    if ROOT.gPad.GetLogy():
        y1, y2 = 10 ** uy1, 10 ** uy2
    else:
        y1, y2 = uy1, uy2
    ln = ROOT.TLine(x, y1, x, y2)
    ln.SetLineColor(color)
    ln.SetLineStyle(2)
    ln.SetLineWidth(2)
    ln.Draw("SAME")
    _keep.append(ln)
    return ln


def data_delta_nu_from_2d(cv2d_data, template_cv):
    """STOPGAP: pooled 1D data delta_nu reconstructed from the 2D data map.

    Projects the (E_nu, E_had) data boundary histogram over all E_nu, shifting
    each column by -nu_cut(E_nu).  Filled onto `template_cv`'s binning (the MC
    delta_nu binning) so a bin-by-bin data/MC ratio is well defined.  CV only.
    """
    h = template_cv.Clone("data_delta_nu_from2d")
    h.Reset()
    nbins = h.GetNbinsX()
    err2 = {}
    for ix in range(1, cv2d_data.GetNbinsX() + 1):
        ex = cv2d_data.GetXaxis().GetBinCenter(ix)
        nc = nu_cut(ex)
        for iy in range(1, cv2d_data.GetNbinsY() + 1):
            val = cv2d_data.GetBinContent(ix, iy)
            if val == 0.0:
                continue
            err = cv2d_data.GetBinError(ix, iy)
            dnu = cv2d_data.GetYaxis().GetBinCenter(iy) - nc
            b = h.FindBin(dnu)
            if b < 1 or b > nbins:
                continue
            h.SetBinContent(b, h.GetBinContent(b) + val)
            err2[b] = err2.get(b, 0.0) + err * err
    for b, e2 in err2.items():
        h.SetBinError(b, e2 ** 0.5)
    return h


# --------------------------------------------------------------------------- #
# 1D delta_nu (reco) -- the headline: data/MC + ratio + systematics breakdown
# --------------------------------------------------------------------------- #
def plot_delta_nu_datamc(data_cv, mc_cv, mc_frac, mc_scale, out_path):
    """POT-normalized data/MC overlay with the MC systematic error band.

    All inputs are plain TH1 CVs on a common binning (see main).  In stopgap
    mode BOTH data and mc_cv come from the 2D projection, so they share the same
    nu < 2.5 GeV truncation -- this is what keeps the overlay honest (otherwise
    a capped data vs an uncapped native MC produces a spurious high-Delta-nu
    deficit).  `mc_frac` is the native (synced) fractional systematic, drawn as
    the band.
    """
    mc = mc_cv.Clone("dnu_ov_mc")
    mc.Scale(mc_scale)
    band = mc.Clone("dnu_ov_band")
    _apply_frac_band(band, mc, mc_frac)
    data = data_cv.Clone("dnu_ov_data")
    data.SetMarkerStyle(20)
    data.SetMarkerColor(ROOT.kBlack)
    data.SetLineColor(ROOT.kBlack)

    canvas = ROOT.TCanvas("c_dnu_datamc", "", 1000, 750)
    if mc.GetMaximum() > 0 or data.GetMaximum() > 0:   # guard: don't log an empty frame
        canvas.SetLogy()
    band.SetFillColorAlpha(ROOT.kRed - 9, 0.5)
    band.SetLineColor(ROOT.kRed + 1)
    band.SetMarkerSize(0)
    band.SetTitle("#Delta#nu (reco)")
    band.GetXaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    band.GetYaxis().SetTitle("Events / bin")
    band.Draw("E2")
    mc.SetLineColor(ROOT.kRed + 1)
    mc.SetLineWidth(2)
    mc.Draw("HIST SAME")
    data.Draw("PE SAME")
    leg = ROOT.TLegend(0.62, 0.74, 0.88, 0.88)
    leg.AddEntry(data, "Data", "PE")
    leg.AddEntry(band, "Simulation", "LF")
    leg.Draw()
    draw_vline(0.0)
    _keep.extend([mc, band, data, leg])
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_delta_nu_ratio(data_cv, mc_cv, mc_frac, mc_scale, out_path):
    """data/MC ratio (POT-normalized) with the MC systematic band around 1.

    data_cv and mc_cv must be on a common binning with a shared truncation
    (see plot_delta_nu_datamc / main).
    """
    m = mc_cv.Clone("dnu_ratio_mc")
    m.Scale(mc_scale)
    ratio = data_cv.Clone("dnu_ratio")
    ratio.Divide(m)  # data/(mc*scale), stat errors propagated

    band = None
    if mc_frac is not None:
        band = ratio.Clone("dnu_ratio_band")
        for b in range(0, band.GetNbinsX() + 2):
            band.SetBinContent(b, 1.0)
            band.SetBinError(b, mc_frac.GetBinContent(b))

    canvas = ROOT.TCanvas("c_dnu_ratio", "", 1000, 750)
    ratio.SetTitle("#Delta#nu  data/MC (POT-normalized)")
    ratio.GetXaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    ratio.GetYaxis().SetTitle("data / MC")
    ratio.GetYaxis().SetRangeUser(0.0, 2.0)
    ratio.SetMarkerStyle(20)
    ratio.SetLineColor(ROOT.kBlack)
    if band is not None:
        band.SetFillColorAlpha(ROOT.kBlue - 9, 0.6)
        band.SetLineColor(ROOT.kBlue - 9)
        band.SetMarkerSize(0)
        band.SetTitle(ratio.GetTitle())
        band.GetXaxis().SetTitle(ratio.GetXaxis().GetTitle())
        band.GetYaxis().SetTitle("data / MC")
        band.GetYaxis().SetRangeUser(0.0, 2.0)
        band.Draw("E2")
        ratio.Draw("PE SAME")
    else:
        ratio.Draw("PE")
    one = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                     ratio.GetXaxis().GetXmax(), 1.0)
    one.SetLineStyle(2)
    one.Draw("SAME")
    draw_vline(0.0)
    leg = ROOT.TLegend(0.60, 0.78, 0.88, 0.88)
    leg.AddEntry(ratio, "data / MC", "PE")
    if band is not None:
        leg.AddEntry(band, "MC systematic band", "F")
    leg.Draw()
    _keep.extend([ratio, band, one, leg])
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_delta_nu_errorsummary(h_mc, out_path):
    """Grouped fractional systematic uncertainty for reco delta_nu."""
    try:
        h_mc.GetXaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    except Exception:
        pass
    canvas = ROOT.TCanvas("c_dnu_errsum", "", 1000, 750)
    try:
        # Same signature as plottingClasses.localDrawErrorSummary, minus the box.
        plotter.DrawErrorSummary(h_mc, "TR", True, True, 0.00001, False,
                                 "", True, "", True)
    except Exception as e:
        print("  WARN DrawErrorSummary failed: {0}".format(e))
        return
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_delta_nu_truth(h_mc_true, out_path):
    """1D delta_nu (truth) -- MC only; data has no truth."""
    cv = get_cv(h_mc_true)
    cv.SetLineColor(ROOT.kBlue + 2)
    cv.SetLineWidth(2)
    canvas = ROOT.TCanvas("c_dnu_truth", "", 1000, 750)
    if cv.GetMaximum() > 0:
        canvas.SetLogy()
    cv.SetTitle("#nu - #nu_{cut}(E_{#nu}) (truth, MC)")
    cv.GetXaxis().SetTitle("#Delta#nu_{true}  [GeV]")
    cv.GetYaxis().SetTitle("Events / bin")
    cv.Draw("HIST")
    draw_vline(0.0)
    _keep.append(cv)
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


# --------------------------------------------------------------------------- #
# 2D: original maps, reconstructed Delta-nu vs E_nu (MC with systematic bands),
#     its fractional-uncertainty map, and E_nu slices (MC band + error summary)
# --------------------------------------------------------------------------- #
def plot_enu_ehad_boundary(h_2d, out_path, title_extra=""):
    """Original (E_nu, nu) heatmap with the staircase overlay (CV)."""
    ROOT.gStyle.SetPalette(DEFAULT_PALETTE)   # Rob's default 2D colour scale
    cv = get_cv(h_2d)
    canvas = ROOT.TCanvas("c_enu_ehad_boundary", "", 900, 700)
    canvas.SetLogz()
    canvas.SetRightMargin(0.13)
    cv.SetTitle("E_{{#nu}} vs #nu near the staircase cut{0}".format(title_extra))
    cv.GetXaxis().SetTitle("E_{#nu}  [GeV]")
    cv.GetYaxis().SetTitle("#nu = E_{had}  [GeV]")
    cv.Draw("COLZ")
    line = staircase_polyline()
    line.Draw("SAME")
    _keep.extend([cv, line])
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def _new_dnu_enu_th2(name):
    """Empty (E_nu nominal, Delta-nu fine) TH2D with the standard axes/titles."""
    xedges = array('d', NOMINAL_ENU_EDGES)
    nx = len(NOMINAL_ENU_EDGES) - 1
    dnu_lo, dnu_hi, dnu_w = -2.0, 2.5, 0.05
    ny = int(round((dnu_hi - dnu_lo) / dnu_w))
    out = ROOT.TH2D(name, "", nx, xedges, ny, dnu_lo, dnu_hi)
    out.SetDirectory(0)
    out.GetXaxis().SetTitle("E_{#nu}  [GeV]")
    out.GetYaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    return out


def _remap_dnu_vs_enu(src2d, name, track_errors=True):
    """Column-shift remap of one (E_nu, E_had) TH2 into (E_nu nominal, Delta-nu).

    Each 0.25-GeV E_nu column is remapped to the nominal E_nu bin it falls in,
    and its nu (=E_had) axis is shifted by -nu_cut(E_nu) to give Delta-nu.
    nu_cut is constant within any nominal bin (corners are nominal edges).

    Works on any TH2 -- a CV hist OR a single systematic universe.  Statistical
    errors are summed in quadrature only when track_errors is set (the CV path).
    For universe hists the error band is built from the spread of bin CONTENTS,
    so per-universe stat errors are irrelevant and skipped for speed.
    """
    out = _new_dnu_enu_th2(name)
    err2 = {} if track_errors else None
    for ix in range(1, src2d.GetNbinsX() + 1):
        ex = src2d.GetXaxis().GetBinCenter(ix)
        if ex >= NOMINAL_ENU_EDGES[-1]:      # outside [0, 15]
            continue
        nc = nu_cut(ex)
        bx = out.GetXaxis().FindBin(ex)
        for iy in range(1, src2d.GetNbinsY() + 1):
            val = src2d.GetBinContent(ix, iy)
            err = src2d.GetBinError(ix, iy)
            if val == 0.0 and err == 0.0:
                continue
            ey = src2d.GetYaxis().GetBinCenter(iy)   # nu = E_had
            by = out.GetYaxis().FindBin(ey - nc)
            out.SetBinContent(bx, by, out.GetBinContent(bx, by) + val)
            if err2 is not None:
                key = (bx, by)
                err2[key] = err2.get(key, 0.0) + err * err
    if err2 is not None:
        for (bx, by), e2 in err2.items():
            out.SetBinError(bx, by, e2 ** 0.5)
    return out


def build_dnu_vs_enu(cv2d, name):
    """CV TH2D of (E_nu nominal, Delta-nu fine) from the (E_nu, E_had) CV map.

    Thin wrapper over _remap_dnu_vs_enu with statistical-error tracking on.
    Used for the DATA map (stat-only) and as the CV base of the MC MnvH2D.
    """
    return _remap_dnu_vs_enu(cv2d, name, track_errors=True)


def build_dnu_vs_enu_mnv(h2_mnv, name):
    """MnvH2D of (E_nu nominal, Delta-nu) carrying the MC vertical systematics.

    Builds the CV by remapping GetCVHistoWithStatError(), then transports EACH
    vertical universe of EACH error band through the IDENTICAL column-shift and
    reassembles them into vertical error bands on the output MnvH2D.  Because
    every universe is moved by the same deterministic remap as the CV, the band
    spread (hence the systematic covariance) is preserved bin-for-bin.

    If the input carries no vertical bands (e.g. a pre-re-run, unsynced 2D map)
    the result is a CV-only MnvH2D and every downstream systematic draw degrades
    gracefully to nothing.
    """
    cv = build_dnu_vs_enu(get_cv(h2_mnv), name)
    out = PlotUtils.MnvH2D(cv)
    out.SetDirectory(0)
    try:
        band_names = list(h2_mnv.GetVertErrorBandNames())
    except Exception as e:
        print("  WARN GetVertErrorBandNames failed ({0}); CV-only 2D map.".format(e))
        band_names = []
    if not band_names:
        print("  NOTE: 2D map '{0}' has no vertical systematic bands "
              "(unsynced/old file?) -> reconstruction is CV-only.".format(name))
    for bname in band_names:
        band = h2_mnv.GetVertErrorBand(bname)
        nuniv = band.GetNHists()
        vec = ROOT.std.vector("TH2D*")()
        for iu in range(nuniv):
            uh = _remap_dnu_vs_enu(band.GetHist(iu),
                                   "{0}_{1}_u{2}".format(name, bname, iu),
                                   track_errors=False)
            _univ_keep.append(uh)     # must outlive `out` (see _univ_keep note)
            vec.push_back(uh)
        out.AddVertErrorBand(bname, vec)
    return out


def plot_dnu_vs_enu(th2, out_path, title_extra=""):
    """COLZ Delta-nu vs E_nu (CV) with the cut at Delta-nu = 0 and corner guides.

    Accepts a plain TH2 (data) or an MnvH2D (MC) -- either way the heatmap shows
    the central value; the MC's systematics are visualized separately by
    plot_dnu_vs_enu_fracunc and drawn as bands on the E_nu slices.
    """
    ROOT.gStyle.SetPalette(DEFAULT_PALETTE)   # Rob's default 2D colour scale
    canvas = ROOT.TCanvas("c_dnu_enu", "", 1000, 750)
    canvas.SetLogz()
    canvas.SetRightMargin(0.13)
    th2.SetTitle("#Delta#nu vs E_{{#nu}}{0}".format(title_extra))
    th2.Draw("COLZ")
    xax, yax = th2.GetXaxis(), th2.GetYaxis()
    ln0 = ROOT.TLine(xax.GetXmin(), 0.0, xax.GetXmax(), 0.0)
    ln0.SetLineColor(ROOT.kRed + 2)
    ln0.SetLineStyle(2)
    ln0.SetLineWidth(3)
    ln0.Draw("SAME")
    _keep.append(ln0)
    for ec in ENU_CORNERS:
        lv = ROOT.TLine(ec, yax.GetXmin(), ec, yax.GetXmax())
        lv.SetLineColor(ROOT.kGray + 2)
        lv.SetLineStyle(3)
        lv.Draw("SAME")
        _keep.append(lv)
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def _has_vert_bands(mnv):
    """True iff `mnv` is an MnvHxD carrying at least one vertical error band."""
    try:
        return len(mnv.GetVertErrorBandNames()) > 0
    except Exception:
        return False


def plot_dnu_vs_enu_fracunc(mnv2d, out_path, title_extra=""):
    """COLZ total fractional systematic uncertainty of the MC reconstruction.

    Directly visualizes the newly-synced 2D bands: for each (E_nu, Delta-nu)
    bin, the total (all-band) fractional systematic uncertainty.  No-op on a
    CV-only map (unsynced/old file).
    """
    if not _has_vert_bands(mnv2d):
        print("  (no 2D systematic bands -> skipping fractional-uncertainty map)")
        return
    frac = mnv2d.GetTotalError(False, True, False)   # sys-only, fractional TH2D
    frac.SetName("dnu_vs_enu_mc_fracunc")
    ROOT.gStyle.SetPalette(DEFAULT_PALETTE)
    canvas = ROOT.TCanvas("c_dnu_enu_fracunc", "", 1000, 750)
    canvas.SetRightMargin(0.13)
    frac.SetTitle("#Delta#nu vs E_{{#nu}} -- MC fractional systematic unc.{0}"
                  .format(title_extra))
    frac.GetXaxis().SetTitle("E_{#nu}  [GeV]")
    frac.GetYaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    frac.GetZaxis().SetTitle("fractional systematic uncertainty")
    frac.Draw("COLZ")
    xax, yax = frac.GetXaxis(), frac.GetYaxis()
    ln0 = ROOT.TLine(xax.GetXmin(), 0.0, xax.GetXmax(), 0.0)
    ln0.SetLineColor(ROOT.kRed + 2)
    ln0.SetLineStyle(2)
    ln0.SetLineWidth(3)
    ln0.Draw("SAME")
    _keep.extend([frac, ln0])
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_dnu_slice(dnu2d_data, dnu2d_mc, mc_scale, region, out_path):
    """Delta-nu data/MC overlay for one staircase E_nu region.

    MC (from the banded MnvH2D) is drawn with its systematic band; DATA points
    carry statistical error bars only -- the 2D data map is stat-only by design
    (SyncAllHists syncs only the MC 2D wrappers), matching the 1D data.
    """
    name, elo, ehi, nc = region
    xax = dnu2d_mc.GetXaxis()
    blo = xax.FindBin(elo + 1e-6)
    bhi = xax.FindBin(ehi - 1e-6)

    # MC: project the banded MnvH2D -> MnvH1D (bands preserved), then POT-scale.
    hm_mnv = dnu2d_mc.ProjectionY("dnu_" + name + "_mc", blo, bhi)
    hm_mnv.Scale(mc_scale)
    hm = hm_mnv.GetCVHistoWithStatError()          # scaled CV (TH1D)
    hd = dnu2d_data.ProjectionY("dnu_" + name + "_data", blo, bhi)  # data (stat)

    # MC systematic band = CV content * per-bin fractional systematic.
    band = None
    if _has_vert_bands(hm_mnv):
        mc_frac = hm_mnv.GetTotalError(False, True, False)  # sys-only, fractional
        band = hm.Clone("dnu_" + name + "_band")
        _apply_frac_band(band, hm, mc_frac)

    title = ("#Delta#nu,  {0:g} < E_{{#nu}} < {1:g} GeV  (#nu_{{cut}} = {2:g})"
             .format(elo, ehi, nc))
    canvas = ROOT.TCanvas("c_slice_" + name, "", 1000, 750)
    if hm.GetMaximum() > 0:
        canvas.SetLogy()
    hm.SetLineColor(ROOT.kBlue + 2)
    hm.SetLineWidth(2)
    if band is not None:
        band.SetFillColorAlpha(ROOT.kBlue - 9, 0.5)
        band.SetLineColor(ROOT.kBlue + 2)
        band.SetMarkerSize(0)
        band.SetTitle(title)
        band.GetXaxis().SetTitle("#Delta#nu  [GeV]")
        band.GetYaxis().SetTitle("Events / bin")
        band.Draw("E2")
        hm.Draw("HIST SAME")
    else:
        hm.SetTitle(title)
        hm.GetXaxis().SetTitle("#Delta#nu  [GeV]")
        hm.GetYaxis().SetTitle("Events / bin")
        hm.Draw("HIST")
    hd.SetMarkerStyle(20)
    hd.SetLineColor(ROOT.kBlack)
    hd.Draw("PE SAME")
    draw_vline(0.0)

    di, mi = hd.Integral(), hm.Integral()
    r = di / mi if mi > 0 else 0.0
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.035)
    tex.DrawLatex(0.16, 0.20, "data/MC (integral) = {0:.3f}".format(r))
    leg = ROOT.TLegend(0.63, 0.74, 0.88, 0.88)
    leg.AddEntry(hd, "Data (stat)", "PE")
    if band is not None:
        leg.AddEntry(band, "MC #pm syst", "LF")
    else:
        leg.AddEntry(hm, "MC (CV)", "L")
    leg.Draw()
    # hm_mnv kept alive too: the band's universes back the projected MnvH1D.
    _keep.extend([hd, hm, hm_mnv, tex, leg])
    if band is not None:
        _keep.append(band)
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_dnu_slice_errorsummary(dnu2d_mc, mc_scale, region, out_path):
    """Grouped fractional systematic breakdown for one staircase E_nu slice.

    The 2D analogue payoff of syncing the boundary map: the same grouped
    error-summary shown for the 1D delta_nu, now per E_nu region.  No-op on a
    CV-only map.
    """
    name, elo, ehi, nc = region
    if not _has_vert_bands(dnu2d_mc):
        print("  (no 2D bands -> skipping slice error summary for {0})".format(name))
        return
    xax = dnu2d_mc.GetXaxis()
    blo = xax.FindBin(elo + 1e-6)
    bhi = xax.FindBin(ehi - 1e-6)
    proj = dnu2d_mc.ProjectionY("dnu_" + name + "_mc_errsum", blo, bhi)
    proj.Scale(mc_scale)
    try:
        proj.GetXaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    except Exception:
        pass
    canvas = ROOT.TCanvas("c_slice_errsum_" + name, "", 1000, 750)
    try:
        plotter.DrawErrorSummary(proj, "TR", True, True, 0.00001, False,
                                 "", True, "", True)
    except Exception as e:
        print("  WARN slice DrawErrorSummary failed ({0})".format(e))
        return
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.03)
    tex.DrawLatex(0.16, 0.92,
                  "{0:g} < E_{{#nu}} < {1:g} GeV".format(elo, ehi))
    _keep.extend([proj, tex])
    canvas.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    in_path, out_dir = parse_args()
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        print("Created output directory: {0}".format(out_dir))

    f = ROOT.TFile.Open(in_path, "READ")
    if not f or f.IsZombie():
        sys.stderr.write("ERROR: could not open {0}\n".format(in_path))
        sys.exit(1)
    print("Reading: {0}".format(in_path))

    mc_pot, data_pot = read_pot(f)
    mc_scale = data_pot / mc_pot

    # ---- load histograms (2D needed before 1D for the data stopgap) ----
    h_mc = safe_get(f, "selection_mc_delta_nu")
    h_data = safe_get(f, "selection_data_delta_nu")
    h_mc_true = safe_get(f, "selection_mc_delta_nu_true")
    h2_data = get_2d(f, "data")
    h2_mc = get_2d(f, "mc")

    # ---- 1D delta_nu (reco): headline data/MC + ratio + error summary ----
    if h_mc:
        # Determine the data CV: native if actually filled, else the stopgap.
        # In stopgap mode the MC must be built the SAME way (2D projection) so
        # data and MC share the nu<2.5 truncation -- otherwise capped data vs
        # uncapped native MC fakes a high-Delta-nu deficit.
        data_cv = None
        if h_data:
            nat = get_cv(h_data)
            if nat.Integral() > 0:
                data_cv = nat
        mc_cv_pooled = get_cv(h_mc)          # native CV (correct once data is native)
        if data_cv is None and h2_data is not None:
            print("  FALLBACK (pre-re-run file): native 1D data delta_nu is empty "
                  "-> reconstructing data AND MC pooled delta_nu from the 2D "
                  "projection so both share the nu < 2.5 GeV truncation "
                  "(apples-to-apples). On synced re-run output native data is "
                  "present and this path is skipped automatically.")
            data_cv = data_delta_nu_from_2d(get_cv(h2_data), get_cv(h_mc))
            if h2_mc is not None:
                mc_cv_pooled = data_delta_nu_from_2d(get_cv(h2_mc), get_cv(h_mc))

        try:
            mc_frac = h_mc.GetTotalError(False, True, False)  # sys-only, fractional
        except Exception as e:
            print("  WARN MC fractional systematic unavailable ({0})".format(e))
            mc_frac = None

        if data_cv is not None:
            plot_delta_nu_datamc(data_cv, mc_cv_pooled, mc_frac, mc_scale,
                                 "{0}/delta_nu_reco_datamc.png".format(out_dir))
            plot_delta_nu_ratio(data_cv, mc_cv_pooled, mc_frac, mc_scale,
                                "{0}/delta_nu_reco_ratio.png".format(out_dir))
        else:
            print("  WARN: no data available for delta_nu; skipping overlay/ratio.")
        plot_delta_nu_errorsummary(h_mc,
                                   "{0}/delta_nu_reco_errorsummary.png".format(out_dir))

    # ---- 1D delta_nu (truth) ----
    if h_mc_true:
        plot_delta_nu_truth(h_mc_true,
                            "{0}/delta_nu_truth.png".format(out_dir))

    # ---- 2D original (E_nu, nu) maps ----
    if h2_data:
        plot_enu_ehad_boundary(h2_data,
                               "{0}/enu_ehad_boundary_data.png".format(out_dir),
                               title_extra=" -- data")
    if h2_mc:
        plot_enu_ehad_boundary(h2_mc,
                               "{0}/enu_ehad_boundary_mc.png".format(out_dir),
                               title_extra=" -- MC reco")

    # ---- reconstructed Delta-nu vs E_nu + slices ----
    # DATA: stat-only CV TH2D.  MC: full MnvH2D whose vertical systematic bands
    # are transported through the same column-shift remap as the CV.
    dnu2d_data = build_dnu_vs_enu(get_cv(h2_data), "dnu_vs_enu_data") if h2_data else None
    dnu2d_mc = build_dnu_vs_enu_mnv(h2_mc, "dnu_vs_enu_mc") if h2_mc else None
    if dnu2d_data:
        plot_dnu_vs_enu(dnu2d_data,
                        "{0}/dnu_vs_enu_data.png".format(out_dir),
                        title_extra=" -- data")
    if dnu2d_mc:
        plot_dnu_vs_enu(dnu2d_mc,
                        "{0}/dnu_vs_enu_mc.png".format(out_dir),
                        title_extra=" -- MC reco")
        plot_dnu_vs_enu_fracunc(dnu2d_mc,
                                "{0}/dnu_vs_enu_mc_fracunc.png".format(out_dir),
                                title_extra=" -- MC reco")
    if dnu2d_data and dnu2d_mc:
        for region in STAIRCASE_REGIONS:
            plot_dnu_slice(dnu2d_data, dnu2d_mc, mc_scale, region,
                           "{0}/dnu_slice_{1}.png".format(out_dir, region[0]))
            plot_dnu_slice_errorsummary(dnu2d_mc, mc_scale, region,
                                        "{0}/dnu_slice_{1}_errorsummary.png"
                                        .format(out_dir, region[0]))

    print("Done.  Plots in {0}".format(out_dir))
    print("NOTE: with synced re-run output the 1D data/MC uses NATIVE data and "
          "the 2D-derived MC plots (reconstruction, fractional-uncertainty map, "
          "and E_nu slices + error summaries) carry systematic bands; the 2D "
          "DATA stays stat-only by design. On a pre-re-run file the 1D falls "
          "back to the 2D-projection data and the 2D bands are absent (CV-only).")


if __name__ == "__main__":
    main()
