"""
boundary_diagnostic.py   (v3.1 -- 2026-10-07)

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
    - a Delta-nu vs E_nu map (data; MC with its vertical systematic bands)
    - a COLZ map of the MC's total fractional systematic uncertainty over
      (E_nu, Delta-nu)
    - Delta-nu distributions in the four staircase E_nu regions (slices), each
      with the MC systematic band drawn on the overlay and a companion grouped
      systematic error-summary
  Low-nu misclassification (needs the 2026-10 event loop):
    - purity LL/(LL+LH) and completeness LL/(LL+HL) of the reco low-nu call vs
      E_nu, and the net reco/truth low-nu ratio (LL+LH)/(LL+HL), each with the
      MC systematic band; once binned in reco E_nu, once in true E_nu
    - a printed cross-check: summed over E_nu the four categories must equal
      the sign quadrants of migration_delta_nu in the same file (that pooled
      table is printed for older files too)

  ==================================================================
  STATE (2026-10-07)

  v3.0 VALIDATED on XSecInputs_1110_ME1A_2026-08-14.root: the 2D MC boundary
  map carries 55 synced vertical bands (band CV == hist CV), native 1D data
  delta_nu is filled (263,949 events), and the fractional-uncertainty map and
  the slice error summaries populate.

  v3.1 reads two new inputs from the 2026-10 event-loop changes, each only
  when present, so older files run exactly as under v3.0:

  (1) NATIVE (E_nu, Delta-nu) map  selection_{data,mc}_enu_delta_nu_inclusive:
      Delta-nu filled per universe, nominal E_nu bins x 0.05 GeV on [-2, 3].
      When present it replaces the column-shift remap of the boundary map for
      the Delta-nu vs E_nu maps, the fractional-uncertainty map and the
      slices.  The top slice becomes 12-20 GeV (the range the NOMAD
      normalization integrates) and is drawn in 0.2 GeV Delta-nu bins for MC
      statistics.  Without it the remap path runs as before: E_nu < 15 GeV,
      nu < 2.5 GeV, top slice 12-15 GeV.
  (2) misclass_lownu_enu_{reco,true}: see "Low-nu misclassification" above.

  The boundary map is now written as selection_{data,mc}_enu_ehad_boundary_inclusive;
  the old auto-joined name (..._enu_boundary_ehad_boundary_...) is still read.

  DATA IS STILL STAT-ONLY ON THE 2D, BY DESIGN.  SyncAllHists(Variable2D&)
  syncs only the MC wrappers; on the slices the DATA points carry statistical
  errors and the band shown is the MC systematic band.
  ==================================================================

Usage:
    python boundary_diagnostic.py <path-to-XSecInputs.root> [<output-dir>]

Histogram name conventions (Histograms.cxx / Variable2D.h; 2D suffix = sample):
  1D: selection_{data,mc}_<varname>   (delta_nu, delta_nu_true); migration_delta_nu
  2D: selection_{data,mc}_enu_ehad_boundary_inclusive
        (before 2026-10: selection_{data,mc}_enu_boundary_ehad_boundary_inclusive)
      selection_{data,mc}_enu_delta_nu_inclusive              (2026-10 on)
  misclass_lownu_enu_reco, misclass_lownu_enu_true            (2026-10 on)
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

# (name, E_nu_lo, E_nu_hi, nu_cut, Delta-nu rebin factor) for the four slices.
# STAIRCASE_REGIONS: remap path (the boundary map stops at E_nu = 15 GeV).
STAIRCASE_REGIONS = [
    ("Enu_lt3",    0.0,  3.0,  0.3, 1),
    ("Enu_3to7",   3.0,  7.0,  0.5, 1),
    ("Enu_7to12",  7.0,  12.0, 1.0, 1),
    ("Enu_12to15", 12.0, 15.0, 2.0, 1),
]
# NATIVE_REGIONS: native enu_delta_nu path.  The top slice runs to the 20 GeV
# nominal edge (the NOMAD normalization integral is 12-20 GeV) and is rebinned
# x4 (0.2 GeV in Delta-nu) because MC statistics are thin there.
NATIVE_REGIONS = [
    ("Enu_lt3",    0.0,  3.0,  0.3, 1),
    ("Enu_3to7",   3.0,  7.0,  0.5, 1),
    ("Enu_7to12",  7.0,  12.0, 1.0, 1),
    ("Enu_12to20", 12.0, 20.0, 2.0, 4),
]
NATIVE_ENU_MAX = 20.0      # displayed E_nu range (GeV) of the native 2D maps
MISCLASS_ENU_MAX = 20.0    # displayed E_nu range (GeV) of the misclassification plots
# y bins of misclass_lownu_enu_*: first letter = reco call, second = truth call
MISCLASS_BINS = [("LL", 1), ("LH", 2), ("HL", 3), ("HH", 4)]


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
        "selection_{0}_enu_ehad_boundary_inclusive".format(which),           # 2026-10 on
        "selection_{0}_enu_boundary_ehad_boundary_inclusive".format(which),  # older files
    ]
    for name in candidates:
        h = f.Get(name)
        if h:
            print("  using {0}".format(name))
            return h
    safe_get(f, candidates[0], fallback_substring="boundary")
    return None


def get_native_dnu(f, which):
    """Native (E_nu, Delta-nu) MnvH2D for 'data' or 'mc' (2026-10 event loop on), else None."""
    name = "selection_{0}_enu_delta_nu_inclusive".format(which)
    h = f.Get(name)
    if h:
        print("  using {0}".format(name))
        return h
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


def _cv_copy(h, name):
    """Plain TH1/TH2 copy for drawing: the CV of an MnvHxD, else a clone."""
    try:
        c = h.GetCVHistoWithStatError()
    except AttributeError:
        c = h.Clone()
    c.SetName(name)
    c.SetDirectory(0)
    return c


def _enu_display_max(xax, enu_max):
    """Upper E_nu shown: enu_max if given (native axis runs to 120 GeV), else the axis top."""
    return xax.GetXmax() if enu_max is None else min(enu_max, xax.GetXmax())


def plot_dnu_vs_enu(th2, out_path, title_extra="", enu_max=None):
    """COLZ Delta-nu vs E_nu (CV) with the cut at Delta-nu = 0 and corner guides.

    Accepts a plain TH2 (data) or an MnvH2D (MC) -- either way the heatmap shows
    the central value; the MC's systematics are visualized separately by
    plot_dnu_vs_enu_fracunc and drawn as bands on the E_nu slices.  enu_max
    limits the displayed E_nu range.
    """
    ROOT.gStyle.SetPalette(DEFAULT_PALETTE)   # Rob's default 2D colour scale
    h = _cv_copy(th2, "draw_" + os.path.splitext(os.path.basename(out_path))[0])
    canvas = ROOT.TCanvas("c_dnu_enu", "", 1000, 750)
    canvas.SetLogz()
    canvas.SetRightMargin(0.13)
    h.SetTitle("#Delta#nu vs E_{{#nu}}{0}".format(title_extra))
    h.GetXaxis().SetTitle("E_{#nu}  [GeV]")
    h.GetYaxis().SetTitle("#Delta#nu = #nu - #nu_{cut}(E_{#nu})  [GeV]")
    xax, yax = h.GetXaxis(), h.GetYaxis()
    x_hi = _enu_display_max(xax, enu_max)
    if enu_max is not None:
        xax.SetRangeUser(xax.GetXmin(), x_hi - 1e-6)
    h.Draw("COLZ")
    ln0 = ROOT.TLine(xax.GetXmin(), 0.0, x_hi, 0.0)
    ln0.SetLineColor(ROOT.kRed + 2)
    ln0.SetLineStyle(2)
    ln0.SetLineWidth(3)
    ln0.Draw("SAME")
    _keep.extend([h, ln0])
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


def plot_dnu_vs_enu_fracunc(mnv2d, out_path, title_extra="", enu_max=None):
    """COLZ total fractional systematic uncertainty of the MC Delta-nu vs E_nu map.

    For each (E_nu, Delta-nu) bin, the total (all-band) fractional systematic
    uncertainty.  No-op on a CV-only map (unsynced/old file).
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
    xax, yax = frac.GetXaxis(), frac.GetYaxis()
    x_hi = _enu_display_max(xax, enu_max)
    if enu_max is not None:
        xax.SetRangeUser(xax.GetXmin(), x_hi - 1e-6)
    frac.Draw("COLZ")
    ln0 = ROOT.TLine(xax.GetXmin(), 0.0, x_hi, 0.0)
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
    name, elo, ehi, nc, rebin = region
    xax = dnu2d_mc.GetXaxis()
    blo = xax.FindBin(elo + 1e-6)
    bhi = xax.FindBin(ehi - 1e-6)

    # MC: project the banded MnvH2D -> MnvH1D (bands preserved), then POT-scale.
    hm_mnv = dnu2d_mc.ProjectionY("dnu_" + name + "_mc", blo, bhi)
    hm_mnv.Scale(mc_scale)
    hd = dnu2d_data.ProjectionY("dnu_" + name + "_data", blo, bhi)  # data (stat)
    if rebin > 1:
        hm_mnv.Rebin(rebin)   # in place; MnvH1D::Rebin carries the error bands
        hd.Rebin(rebin)
    hm = hm_mnv.GetCVHistoWithStatError()          # scaled CV (TH1D)
    ytitle = ("Events / bin" if rebin == 1 else
              "Events / {0:g} GeV".format(hm.GetXaxis().GetBinWidth(1)))

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
        band.GetYaxis().SetTitle(ytitle)
        band.Draw("E2")
        hm.Draw("HIST SAME")
    else:
        hm.SetTitle(title)
        hm.GetXaxis().SetTitle("#Delta#nu  [GeV]")
        hm.GetYaxis().SetTitle(ytitle)
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
    name, elo, ehi, nc, rebin = region
    if not _has_vert_bands(dnu2d_mc):
        print("  (no 2D bands -> skipping slice error summary for {0})".format(name))
        return
    xax = dnu2d_mc.GetXaxis()
    blo = xax.FindBin(elo + 1e-6)
    bhi = xax.FindBin(ehi - 1e-6)
    proj = dnu2d_mc.ProjectionY("dnu_" + name + "_mc_errsum", blo, bhi)
    proj.Scale(mc_scale)
    if rebin > 1:
        proj.Rebin(rebin)     # in place; bands follow
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
# Low-nu misclassification (reco call vs truth call), 2026-10 event loop on
# --------------------------------------------------------------------------- #
def migration_quadrants(m):
    """(LL, LH, HL, HH) CV sums of migration_delta_nu (x = reco, y = true Delta-nu).

    The sign of Delta-nu is the low-nu call (IsLowNu: nu < nu_cut, so Delta-nu = 0
    counts as high-nu, matching the bin that starts at 0), hence the four sign
    quadrants -- under/overflow included -- are the misclassification table
    pooled over E_nu.  Selected signal MC, CV weights, not POT-scaled.
    """
    b0x = m.GetXaxis().FindBin(0.0)
    b0y = m.GetYaxis().FindBin(0.0)
    nx, ny = m.GetNbinsX(), m.GetNbinsY()
    return (m.Integral(0, b0x - 1, 0, b0y - 1), m.Integral(0, b0x - 1, b0y, ny + 1),
            m.Integral(b0x, nx + 1, 0, b0y - 1), m.Integral(b0x, nx + 1, b0y, ny + 1))


def misclass_totals(h2):
    """(LL, LH, HL, HH) CV totals of a misclass_lownu_enu_* map over all E_nu."""
    nx = h2.GetNbinsX()
    return tuple(h2.Integral(0, nx + 1, b, b) for _, b in MISCLASS_BINS)


def _print_misclass_row(label, q):
    ll, lh, hl, hh = q
    pur = ll / (ll + lh) if (ll + lh) > 0 else 0.0
    comp = ll / (ll + hl) if (ll + hl) > 0 else 0.0
    net = (ll + lh) / (ll + hl) if (ll + hl) > 0 else 0.0
    print("    {0:<26s}{1:>13.1f}{2:>13.1f}{3:>13.1f}{4:>14.1f}   {5:.4f}  {6:.4f}  {7:.4f}"
          .format(label, ll, lh, hl, hh, pur, comp, net))


def misclass_check(f, h_mis):
    """Print the pooled table from migration_delta_nu and, if present, the
    misclass_lownu_enu_* totals, which must agree with it exactly."""
    m = f.Get("migration_delta_nu")
    print("  Low-nu misclassification, selected signal MC (CV, not POT-scaled):")
    print("    {0:<26s}{1:>13s}{2:>13s}{3:>13s}{4:>14s}   {5:>6s}  {6:>6s}  {7:>6s}"
          .format("", "LL", "LH", "HL", "HH", "purity", "compl.", "net"))
    q_ref = None
    if m:
        q_ref = migration_quadrants(m)
        _print_misclass_row("migration_delta_nu", q_ref)
    else:
        print("    WARN: migration_delta_nu not found")
    for tag, h2 in h_mis:
        q = misclass_totals(h2)
        _print_misclass_row("misclass_lownu_enu_" + tag, q)
        if q_ref is not None:
            worst = max(abs(a - b) / max(abs(b), 1.0) for a, b in zip(q, q_ref))
            print("      -> vs migration_delta_nu: max relative difference {0:.2e}  {1}"
                  .format(worst, "OK" if worst < 1e-6 else "MISMATCH -- investigate"))


def _visible_extent(hists, x_max):
    """(min, max) of content -/+ error over bins with centre < x_max and content > 0."""
    lo, hi = None, None
    for h in hists:
        for b in range(1, h.GetNbinsX() + 1):
            if h.GetXaxis().GetBinCenter(b) >= x_max or h.GetBinContent(b) <= 0:
                continue
            c, e = h.GetBinContent(b), h.GetBinError(b)
            lo = c - e if lo is None else min(lo, c - e)
            hi = c + e if hi is None else max(hi, c + e)
    return lo, hi


def plot_misclass(h2, tag, out_dir):
    """Purity, completeness and net reco/truth ratio of the low-nu call vs E_nu.

    h2 = misclass_lownu_enu_<tag> (tag 'reco' or 'true': which E_nu bins it).
    Each category is projected to an MnvH1D (bands kept), so the ratios carry
    the MC systematic universes; the bands drawn are stat + sys.
    """
    xlabel = ("E_{#nu}^{reco}  [GeV]" if tag == "reco" else "E_{#nu}^{true}  [GeV]")
    c = {}
    for k, b in MISCLASS_BINS:
        c[k] = h2.ProjectionX("misclass_{0}_{1}".format(tag, k), b, b)
    reco_low = c["LL"].Clone("misclass_{0}_recolow".format(tag))
    reco_low.Add(c["LH"])
    true_low = c["LL"].Clone("misclass_{0}_truelow".format(tag))
    true_low.Add(c["HL"])
    pur = c["LL"].Clone("misclass_{0}_purity".format(tag))
    pur.Divide(c["LL"], reco_low, 1.0, 1.0, "B")
    comp = c["LL"].Clone("misclass_{0}_completeness".format(tag))
    comp.Divide(c["LL"], true_low, 1.0, 1.0, "B")
    net = reco_low.Clone("misclass_{0}_net".format(tag))
    net.Divide(reco_low, true_low)
    _keep.extend(list(c.values()) + [reco_low, true_low, pur, comp, net])

    def band_and_line(h, color, name):
        band = h.GetCVHistoWithError()          # stat + sys
        band.SetName(name + "_band")
        band.SetDirectory(0)
        band.SetFillColorAlpha(color, 0.35)
        band.SetLineColor(color)
        band.SetMarkerSize(0)
        line = h.GetCVHistoWithStatError()
        line.SetName(name + "_line")
        line.SetDirectory(0)
        line.SetLineColor(color)
        line.SetLineWidth(2)
        for hh in (band, line):
            hh.GetXaxis().SetRangeUser(hh.GetXaxis().GetXmin(), MISCLASS_ENU_MAX - 1e-6)
        _keep.extend([band, line])
        return band, line

    # -- purity + completeness
    canvas = ROOT.TCanvas("c_misclass_pc_" + tag, "", 1000, 750)
    pb, pl = band_and_line(pur, ROOT.kBlue + 1, "pur_" + tag)
    cb, cl = band_and_line(comp, ROOT.kRed + 1, "comp_" + tag)
    lo, hi = _visible_extent([pb, cb], MISCLASS_ENU_MAX)
    lo = 0.0 if lo is None else max(0.0, int(lo * 10.0) / 10.0 - 0.05)
    pb.SetTitle("Low-#nu call: purity and completeness ({0} E_{{#nu}})".format(tag))
    pb.GetXaxis().SetTitle(xlabel)
    pb.GetYaxis().SetTitle("fraction")
    pb.GetYaxis().SetRangeUser(lo, 1.05)
    pb.Draw("E2")
    cb.Draw("E2 SAME")
    pl.Draw("HIST SAME")
    cl.Draw("HIST SAME")
    for ec in ENU_CORNERS:
        draw_vline(ec, color=ROOT.kGray + 2)
    leg = ROOT.TLegend(0.45, 0.16, 0.88, 0.32)
    leg.AddEntry(pb, "purity  LL/(LL+LH)", "LF")
    leg.AddEntry(cb, "completeness  LL/(LL+HL)", "LF")
    leg.Draw()
    _keep.append(leg)
    out = "{0}/misclass_lownu_enu_{1}_purity_completeness.png".format(out_dir, tag)
    canvas.SaveAs(out)
    print("  wrote {0}".format(out))

    # -- net reco/truth low-nu ratio
    canvas = ROOT.TCanvas("c_misclass_net_" + tag, "", 1000, 750)
    nb, nl = band_and_line(net, ROOT.kGreen + 2, "net_" + tag)
    lo, hi = _visible_extent([nb], MISCLASS_ENU_MAX)
    lo = 0.9 if lo is None else min(0.9, lo - 0.02)
    hi = 1.1 if hi is None else max(1.1, hi + 0.02)
    nb.SetTitle("Low-#nu call: reco low-#nu / truth low-#nu ({0} E_{{#nu}})".format(tag))
    nb.GetXaxis().SetTitle(xlabel)
    nb.GetYaxis().SetTitle("(LL+LH) / (LL+HL)")
    nb.GetYaxis().SetRangeUser(lo, hi)
    nb.Draw("E2")
    nl.Draw("HIST SAME")
    one = ROOT.TLine(nb.GetXaxis().GetXmin(), 1.0, MISCLASS_ENU_MAX, 1.0)
    one.SetLineStyle(2)
    one.Draw("SAME")
    for ec in ENU_CORNERS:
        draw_vline(ec, color=ROOT.kGray + 2)
    _keep.append(one)
    out = "{0}/misclass_lownu_enu_{1}_net.png".format(out_dir, tag)
    canvas.SaveAs(out)
    print("  wrote {0}".format(out))


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

    # ---- Delta-nu vs E_nu + slices ----
    # The native (E_nu, Delta-nu) map when the file has one (2026-10 event loop
    # on); otherwise the column-shift remap of the (E_nu, E_had) boundary map.
    # DATA: stat-only CV TH2D.  MC: MnvH2D carrying its vertical systematic bands.
    h2n_data = get_native_dnu(f, "data")
    h2n_mc = get_native_dnu(f, "mc")
    if h2n_data is not None and h2n_mc is not None:
        print("  NATIVE (E_nu, Delta-nu) map present -> maps and slices use it "
              "(top slice 12-20 GeV, 0.2 GeV Delta-nu bins).")
        dnu2d_data = get_cv(h2n_data)
        dnu2d_mc = h2n_mc
        regions = NATIVE_REGIONS
        enu_max = NATIVE_ENU_MAX
    else:
        print("  native (E_nu, Delta-nu) map absent (pre-2026-10 file) -> column-shift "
              "remap of the boundary map (E_nu < 15 GeV, nu < 2.5 GeV).")
        dnu2d_data = build_dnu_vs_enu(get_cv(h2_data), "dnu_vs_enu_data") if h2_data else None
        dnu2d_mc = build_dnu_vs_enu_mnv(h2_mc, "dnu_vs_enu_mc") if h2_mc else None
        regions = STAIRCASE_REGIONS
        enu_max = None
    if dnu2d_data:
        plot_dnu_vs_enu(dnu2d_data,
                        "{0}/dnu_vs_enu_data.png".format(out_dir),
                        title_extra=" -- data", enu_max=enu_max)
    if dnu2d_mc:
        plot_dnu_vs_enu(dnu2d_mc,
                        "{0}/dnu_vs_enu_mc.png".format(out_dir),
                        title_extra=" -- MC reco", enu_max=enu_max)
        plot_dnu_vs_enu_fracunc(dnu2d_mc,
                                "{0}/dnu_vs_enu_mc_fracunc.png".format(out_dir),
                                title_extra=" -- MC reco", enu_max=enu_max)
    if dnu2d_data and dnu2d_mc:
        for region in regions:
            plot_dnu_slice(dnu2d_data, dnu2d_mc, mc_scale, region,
                           "{0}/dnu_slice_{1}.png".format(out_dir, region[0]))
            plot_dnu_slice_errorsummary(dnu2d_mc, mc_scale, region,
                                        "{0}/dnu_slice_{1}_errorsummary.png"
                                        .format(out_dir, region[0]))

    # ---- low-nu misclassification ----
    h_mis = []
    for tag in ("reco", "true"):
        h = f.Get("misclass_lownu_enu_" + tag)
        if h:
            h_mis.append((tag, h))
    misclass_check(f, h_mis)
    if h_mis:
        for tag, h2 in h_mis:
            plot_misclass(h2, tag, out_dir)
    else:
        print("  misclass_lownu_enu_* absent (pre-2026-10 file) -> no per-E_nu "
              "misclassification plots; pooled table above is from migration_delta_nu.")

    print("Done.  Plots in {0}".format(out_dir))
    print("NOTE: 2026-10 event-loop output uses the native (E_nu, Delta-nu) map and "
          "adds the misclassification plots; older files fall back to the boundary-map "
          "remap (and, before 2026-08, to the 2D-projection 1D data and CV-only 2D). "
          "The 2D DATA stays stat-only by design.")


if __name__ == "__main__":
    main()
