"""
nucut_efficiency.py   (v1.0 -- 2026-11-20)

The efficiency of a nu cut in slices of muon energy (Kevin, 2026-11-13).  In
each E_mu slice, the fraction of selected events below a recoil threshold,

    F(nu_c) = N(nu < nu_c) / N(all),

at many thresholds nu_c, with every event in the total (nu > 5 GeV included).
Data, MC and data/MC, with the MC systematic uncertainty on F taken from the
universes: each universe's own fraction (its numerator over its own total), so
numerator and denominator move together and a pure normalization (the flux
normalization, for example) drops out.  Does the band cover the data?

Inputs, from the raw event-loop output (XSecInputs_*.root), merged:
  selection_{data,mc}_emu_ehad_fine_<sample>   nu in 0.1 GeV bins to 5 GeV (every
                                               file since 2026-10)
  selection_{data,mc}_emu_ehad_nucut_<sample>  nu in 0.01 GeV bins to 0.6 GeV and
                                               0.05 GeV bins to 2.4 GeV (patch 0036;
                                               the runs of 2026-11-13 on)
Thresholds: from emu_ehad_nucut up to 2.4 GeV (every 0.01 GeV to 0.6, every
0.05 GeV to 2.4), then from emu_ehad_fine every 0.1 GeV to 5 GeV.  Without
emu_ehad_nucut (older files, or --no-nucut) every threshold comes from the
fine map.  Where both maps have a threshold the log prints how well they agree:
they hold the same events, so data, MC and every universe must agree exactly.
The MC is all selected MC (signal + background), as in recoil_emu_slices.py.

Statistical errors on F are binomial (weighted for the MC):
    var F = [(1-F)^2 S2(below) + F^2 S2(above)] / W^2,
with W the total weight and S2 the sum of squared weights.  Neighbouring
thresholds share events, so the points of a curve are correlated: read the
curve, not the scatter between points.

Per E_mu slice it writes (Emu_<lo>to<hi>):
  nucut_eff_<slice>.png          F against nu_c (0-5 GeV): data, MC +- syst;
                                 data/MC with the band around 1
  nucut_eff_<slice>_zoom.png     the same for nu_c < 1.2 GeV (nu_0 = 0.3, 0.5, 1.0)
  nucut_eff_<slice>_errsum.png   the fractional uncertainty on F by group
  breakdowns/nucut_eff_<slice>_<group>.png   the knobs inside GENIE, Model_Fits,
                                 Response and GEANT (those present in the file;
                                 knobs under 0.1 % everywhere are not drawn)
and per staircase value nu_0 = 0.3, 0.5, 1.0, 2.0 (tag nu0_0p3 ...):
  nucut_eff_vs_Emu_<nu0>.png     F(nu < nu_0) against E_mu, data/MC below; the
                                 shaded E_mu range is where nu_0 is the staircase
                                 cut (E_mu + nu_0 inside its E_nu range)
  nucut_eff_vs_Emu_<nu0>_errsum.png, breakdowns/nucut_eff_vs_Emu_<nu0>_<group>.png
and over all slices:
  nucut_eff_ratio_grid.png, nucut_eff_ratio_grid_zoom.png   data/MC of F per slice
  nucut_efficiency.root          the curves (data TH1D, MC MnvH1D with its bands)
                                 per slice and per nu_0, for overlays across runs
  nucut_efficiency_table.txt     per slice: F at each nu_0 and nu_0 +- 10 % (where
                                 the file has that threshold): data, MC, data/MC,
                                 the MC syst and the pull
                                 (F_data - F_MC) / sqrt(syst^2 + stat_data^2 + stat_MC^2)
The log prints the table at the four nu_0 and a count of |pull| > 1 and > 2.

Usage (GPVM, from plotting/, which holds myPlotStyle.h; py_classes/ comes from
setup_mat.sh):
    python3 nucut_efficiency.py <XSecInputs.root> [<output-dir>]
        [--sample inclusive|lowNu|highNu]   (default inclusive)
        [--emu-edges 1,2,3,...,20]          slice edges in GeV, each an E_mu bin
                                            edge of the file (default: every file
                                            bin in 1-20 GeV)
        [--label TEXT]                      run label on the plots (default: the
                                            tune / muon_shift / selection stamps)
        [--no-nucut]                        thresholds from the fine map only

Drafted 2026-11-20 (card A6).  Co-authored by Rob Fine and Claude (Anthropic).
"""

import argparse
from array import array
from collections import OrderedDict
import math
import os
import sys

import ROOT
from functions import *          # brings in PlotUtils (as recoil_emu_slices.py)
from errorMaps import *          # error_bands (OrderedDict)

ROOT.gROOT.SetBatch()
try:
    ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
    ROOT.myPlotStyle()
except Exception:
    pass
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

RESPONSE_BANDS = ["response_proton", "response_meson", "response_em", "response_other",
                  "response_low_neutron", "response_mid_neutron", "response_high_neutron"]
if not any(n in sum([list(v) for v in error_bands.values()], []) for n in RESPONSE_BANDS):
    error_bands["Response"] = RESPONSE_BANDS
BREAKDOWN_GROUPS = ["GENIE", "Model_Fits", "Response", "GEANT"]
BREAKDOWN_THRESHOLD = 0.001   # F is a ratio, so its uncertainties are small

plotter = PlotUtils.MnvPlotter()
plotter.error_summary_group_map.clear()
for group in error_bands:
    for error in error_bands[group]:
        plotter.error_summary_group_map[group].push_back(error)
plotter.SetLegendNColumns(2)

_keep = []   # Python-side references so ROOT drawables outlive SaveAs

# Staircase, mirroring LowNuBoundary.h: (nu_0, E_nu lo, E_nu hi)
STAIRCASE = [(0.3, 0.0, 3.0), (0.5, 3.0, 7.0), (1.0, 7.0, 12.0), (2.0, 12.0, 1.0e9)]
NU0 = [0.3, 0.5, 1.0, 2.0]
DELTA = 0.10           # the table also gives nu_0 * (1 -+ DELTA)
FINE, NUCUT = "emu_ehad_fine", "emu_ehad_nucut"
NUCUT_MAX = 2.4        # GeV: above it the nucut axis has one bin to 5 GeV
NU_MAX = 5.0
ZOOM_MAX = 1.2
EPS = 1e-6
EMU_PLOT_MAX = 20.0


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile")
    p.add_argument("outdir", nargs="?", default="./nucut_efficiency_plots")
    p.add_argument("--sample", default="inclusive", choices=["inclusive", "lowNu", "highNu"])
    p.add_argument("--emu-edges", default=None,
                   help="comma-separated E_mu slice edges in GeV (file bin edges)")
    p.add_argument("--label", default=None)
    p.add_argument("--no-nucut", action="store_true",
                   help="take every threshold from emu_ehad_fine")
    return p.parse_args()


def stamp(f, name):
    o = f.Get(name)
    return o.GetTitle() if o else None


def run_label(f):
    parts = []
    for key in ("tune", "muon_shift", "selection", "geant_hadron"):
        v = stamp(f, key)
        if v:
            parts.append("{0}: {1}".format(key, v.split(" (")[0]))
    return "   ".join(parts)


def get_hist(f, name, required=True):
    h = f.Get(name)
    if h:
        print("  using {0}".format(name))
        return h
    if required:
        print("  ERROR: '{0}' not found.".format(name))
        sys.exit(2)
    print("  NOTE: '{0}' not found.".format(name))
    return None


def slice_edges(xax, requested):
    """Slices as (lo, hi, bin_lo, bin_hi); every edge must be a file bin edge."""
    file_edges = [xax.GetBinLowEdge(b) for b in range(1, xax.GetNbins() + 2)]
    if requested is None:
        edges = [e for e in file_edges if 1.0 - EPS <= e <= 20.0 + EPS]
    else:
        edges = [float(x) for x in requested.split(",")]
        for e in edges:
            if min(abs(e - fe) for fe in file_edges) > EPS:
                print("  ERROR: E_mu edge {0:g} is not a bin edge of the file. File edges: {1}"
                      .format(e, ", ".join("{0:g}".format(x) for x in file_edges)))
                sys.exit(2)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append((lo, hi, xax.FindBin(lo + EPS), xax.FindBin(hi - EPS)))
    return out


def same_x_axis(a, b):
    xa, xb = a.GetXaxis(), b.GetXaxis()
    if xa.GetNbins() != xb.GetNbins():
        return False
    return all(abs(xa.GetBinLowEdge(i) - xb.GetBinLowEdge(i)) < EPS
               for i in range(1, xa.GetNbins() + 2))


def slice_tag(lo, hi):
    if float(lo).is_integer() and float(hi).is_integer():
        return "Emu_{0:02d}to{1:02d}".format(int(lo), int(hi))
    return "Emu_{0:g}to{1:g}".format(lo, hi).replace(".", "p")


def nu0_tag(v):
    return "nu0_{0:g}".format(v).replace(".", "p")


def nu0_lines(lo, hi):
    """Staircase nu_0 values that are the cut somewhere in E_mu in [lo, hi]."""
    return [nu0 for nu0, elo, ehi in STAIRCASE if lo + nu0 < ehi and hi + nu0 > elo]


def thresholds(fine_y, nucut_y):
    """[(nu_c, source, bin k)]: F(nu_c) sums bins 0..k of that source's nu axis."""
    out = []
    lo = 0.0
    if nucut_y is not None:
        for k in range(1, nucut_y.GetNbins() + 1):
            t = nucut_y.GetBinUpEdge(k)
            if t <= NUCUT_MAX + EPS:
                out.append((t, NUCUT, k))
        lo = NUCUT_MAX
    for k in range(1, fine_y.GetNbins() + 1):
        t = fine_y.GetBinUpEdge(k)
        if t > lo + EPS:
            out.append((t, FINE, k))
    return out


# --------------------------------------------------------------------------- #
# The cumulative fraction
# --------------------------------------------------------------------------- #
class Cumul(object):
    """Prefix sums of one nu projection, under- and overflow included."""

    def __init__(self, h, with_errors=True):
        n = h.GetNbinsX() + 2
        c = [h.GetBinContent(b) for b in range(n)]
        self.cum = []
        s = 0.0
        for v in c:
            s += v
            self.cum.append(s)
        self.W = s
        if with_errors:
            e2 = [h.GetBinError(b) ** 2 for b in range(n)]
            self.cum2 = []
            s = 0.0
            for v in e2:
                s += v
                self.cum2.append(s)
            self.S2 = s

    def frac(self, k):
        return self.cum[k] / self.W if self.W > 0 else 0.0

    def var(self, k):
        if self.W <= 0:
            return 0.0
        F = self.cum[k] / self.W
        return ((1.0 - F) ** 2 * self.cum2[k] + F ** 2 * (self.S2 - self.cum2[k])) / self.W ** 2


def bands_of(mnv):
    """OrderedDict (kind, name) -> ([universe TH1D], use_spread) of an MnvH1D."""
    out = OrderedDict()
    for kind in ("Vert", "Lat"):
        for bn in list(getattr(mnv, "Get%sErrorBandNames" % kind)()):
            band = getattr(mnv, "Get%sErrorBand" % kind)(str(bn))
            out[(kind, str(bn))] = ([band.GetHist(i) for i in range(band.GetNHists())],
                                    bool(band.GetUseSpreadError()))
    return out


def project(h2, name, blo, bhi):
    p = h2.ProjectionY(name, blo, bhi)
    p.SetDirectory(0)
    return p


def slice_curves(maps, blo, bhi, tag, ts):
    """F at every threshold for data, MC CV and every MC universe, in one slice.

    maps: {source: (data MnvH2D, mc MnvH2D)}.  Returns a dict, or None if empty.
    """
    cum = {}
    for src, (hd, hm) in maps.items():
        pd = project(hd, "pd_{0}_{1}".format(src, tag), blo, bhi)
        if hasattr(pd, "GetCVHistoWithStatError"):
            pd = pd.GetCVHistoWithStatError()
        pm = project(hm, "pm_{0}_{1}".format(src, tag), blo, bhi)
        cum[src] = {"data": Cumul(pd), "mc": Cumul(pm.GetCVHistoWithStatError()),
                    "univ": OrderedDict((key, ([Cumul(u, False) for u in us], spread))
                                        for key, (us, spread) in bands_of(pm).items())}
    ref = cum[FINE]
    if ref["data"].W <= 0 and ref["mc"].W <= 0:
        return None
    out = {"fd": [], "vd": [], "fm": [], "vm": [], "univ": OrderedDict(),
           "Nd": ref["data"].W, "Wm": ref["mc"].W}
    for key, (us, spread) in ref["univ"].items():
        out["univ"][key] = ([[] for _ in us], spread)
    for t, src, k in ts:
        c = cum[src]
        out["fd"].append(c["data"].frac(k))
        out["vd"].append(c["data"].var(k))
        out["fm"].append(c["mc"].frac(k))
        out["vm"].append(c["mc"].var(k))
        for key, (lists, _) in out["univ"].items():
            ucs = c["univ"][key][0]
            for iu, lst in enumerate(lists):
                lst.append(ucs[iu].frac(k))
    # closure: where both maps have a threshold, the same events give the same F
    out["closure"] = None
    if NUCUT in cum:
        fy = maps[FINE][0].GetYaxis()
        ny = maps[NUCUT][0].GetYaxis()
        worst = [0.0, 0.0, 0.0]
        nshared = 0
        for kf in range(1, fy.GetNbins() + 1):
            t = fy.GetBinUpEdge(kf)
            if t > NUCUT_MAX + EPS:
                continue
            kn = ny.FindBin(t - EPS)
            if abs(ny.GetBinUpEdge(kn) - t) > EPS:
                continue
            nshared += 1
            a, b = cum[FINE], cum[NUCUT]
            worst[0] = max(worst[0], abs(a["data"].frac(kf) - b["data"].frac(kn)))
            worst[1] = max(worst[1], abs(a["mc"].frac(kf) - b["mc"].frac(kn)))
            for key, (ucs, _) in a["univ"].items():
                for iu, uc in enumerate(ucs):
                    worst[2] = max(worst[2], abs(uc.frac(kf) - b["univ"][key][0][iu].frac(kn)))
        out["closure"] = (nshared, worst)
    return out


# --------------------------------------------------------------------------- #
# Histograms of F
# --------------------------------------------------------------------------- #
def centred_edges(ts):
    """Bin edges whose bin centres are the thresholds."""
    e = [ts[0] - 0.5 * (ts[1] - ts[0])]
    for a, b in zip(ts[:-1], ts[1:]):
        e.append(0.5 * (a + b))
    e.append(ts[-1] + 0.5 * (ts[-1] - ts[-2]))
    return e


def make_data(name, edges, vals, variances):
    h = ROOT.TH1D(name, name, len(edges) - 1, array('d', edges))
    h.SetDirectory(0)
    for i, (v, w) in enumerate(zip(vals, variances)):
        h.SetBinContent(i + 1, v)
        h.SetBinError(i + 1, math.sqrt(max(w, 0.0)))
    return h


def make_mc(name, edges, vals, variances, univ):
    """MnvH1D of F: CV with its binomial stat error, one band per input band."""
    h = make_data(name + "_cv", edges, vals, variances)
    m = PlotUtils.MnvH1D(h)
    m.SetName(name)
    m.SetDirectory(0)
    for (kind, bn), (lists, spread) in univ.items():
        vec = ROOT.std.vector("TH1D*")()
        hold = []                  # the band copies the universes when it is added
        for iu, lst in enumerate(lists):
            uh = make_data("{0}_{1}_u{2}".format(name, bn, iu), edges, lst, [0.0] * len(lst))
            hold.append(uh)
            vec.push_back(uh)
        getattr(m, "Add%sErrorBand" % kind)(bn, vec)
        try:
            getattr(m, "Get%sErrorBand" % kind)(bn).SetUseSpreadError(spread)
        except Exception:
            pass
    return m


def sys_frac(m):
    """Fractional systematic error (no stat) of an MnvH1D, 0 where the CV is 0."""
    err = m.GetTotalError(False, False, False)
    frac = err.Clone(m.GetName() + "_sysfrac")
    frac.SetDirectory(0)
    for b in range(0, frac.GetNbinsX() + 2):
        cv = m.GetBinContent(b)
        frac.SetBinContent(b, err.GetBinContent(b) / cv if cv > 0 else 0.0)
        frac.SetBinError(b, 0.0)
    return frac


def groups_present(m):
    present = set(m.GetVertErrorBandNames()) | set(m.GetLatErrorBandNames())
    return [g for g in BREAKDOWN_GROUPS if any(n in present for n in error_bands.get(g, []))]


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #
def two_pads(name):
    c = ROOT.TCanvas(name, "", 900, 900)
    top = ROOT.TPad(name + "_top", "", 0.0, 0.30, 1.0, 1.0)
    bot = ROOT.TPad(name + "_bot", "", 0.0, 0.0, 1.0, 0.30)
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.10)
    top.SetLeftMargin(0.16)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.32)
    bot.SetLeftMargin(0.16)
    top.Draw()
    bot.Draw()
    _keep.extend([c, top, bot])
    return c, top, bot


def vline(x, ymin, ymax, colour=ROOT.kRed + 2, style=2, width=2):
    ln = ROOT.TLine(x, ymin, x, ymax)
    ln.SetLineColor(colour)
    ln.SetLineStyle(style)
    ln.SetLineWidth(width)
    ln.Draw("SAME")
    _keep.append(ln)
    return ln


def ratio_range(r, xlo, xhi, minimum=0.8, maximum=1.2):
    vals = []
    for b in range(1, r.GetNbinsX() + 1):
        x = r.GetBinCenter(b)
        if xlo <= x <= xhi and r.GetBinContent(b) > 0 and r.GetBinError(b) < 0.1:
            vals.append(r.GetBinContent(b))
    lo = min([minimum] + [v - 0.03 for v in vals])
    hi = max([maximum] + [v + 0.03 for v in vals])
    return max(lo, 0.5), min(hi, 1.5)


def style_ratio_axes(h, ylo, yhi, xtitle):
    h.GetYaxis().SetRangeUser(ylo, yhi)
    h.GetYaxis().SetTitle("data / MC")
    h.GetYaxis().SetNdivisions(505)
    h.GetYaxis().SetTitleSize(0.11)
    h.GetYaxis().SetTitleOffset(0.55)
    h.GetYaxis().SetLabelSize(0.10)
    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleSize(0.12)
    h.GetXaxis().SetLabelSize(0.10)
    h.GetXaxis().SetTitleOffset(1.1)


def ratio_and_band(d, cv, frac, name):
    r = d.Clone(name)
    r.SetDirectory(0)
    r.Divide(cv)
    band = r.Clone(name + "_band")
    band.SetDirectory(0)
    for b in range(0, band.GetNbinsX() + 2):
        band.SetBinContent(b, 1.0)
        band.SetBinError(b, frac.GetBinContent(b))
    return r, band


def plot_curve(d, m, out_path, title, label, xlo, xhi, xtitle, ytitle, vlines=(),
               shade=None, ymax=None, smooth=True):
    """F (top) and data/MC (bottom).  d: TH1D, m: MnvH1D.  shade: (x0, x1) or None."""
    tag = os.path.splitext(os.path.basename(out_path))[0]
    errx = ROOT.gStyle.GetErrorX()
    ROOT.gStyle.SetErrorX(0.0 if smooth else errx)
    cv = m.GetCVHistoWithStatError()
    cv.SetDirectory(0)
    frac = sys_frac(m)
    band = cv.Clone(tag + "_bandabs")
    band.SetDirectory(0)
    for b in range(0, band.GetNbinsX() + 2):
        band.SetBinError(b, cv.GetBinContent(b) * frac.GetBinContent(b))
    c, top, bot = two_pads("c_" + tag)

    top.cd()
    band.SetTitle("")
    band.GetXaxis().SetRangeUser(xlo, xhi)
    band.GetXaxis().SetLabelSize(0)
    band.GetYaxis().SetTitle(ytitle)
    band.GetYaxis().SetTitleOffset(1.25)
    band.GetYaxis().SetTitleSize(0.045)
    if ymax is None:
        vis = [max(cv.GetBinContent(b), d.GetBinContent(b)) for b in range(1, cv.GetNbinsX() + 1)
               if xlo <= cv.GetBinCenter(b) <= xhi]
        ymax = (1.15 if smooth else 1.45) * max(vis) if vis and max(vis) > 0 else 1.0
    band.SetMinimum(0.0)
    band.SetMaximum(ymax)
    band.SetFillColorAlpha(ROOT.kRed - 9, 0.6)
    band.SetLineColor(ROOT.kRed + 1)
    band.SetMarkerSize(0)
    band.Draw("E3" if smooth else "E2")
    if shade is not None:
        box = ROOT.TBox(max(shade[0], xlo), 0.0, min(shade[1], xhi), ymax)
        box.SetFillColorAlpha(ROOT.kAzure - 9, 0.35)
        box.SetLineWidth(0)
        box.Draw("SAME")
        band.Draw("E3 SAME" if smooth else "E2 SAME")
        _keep.append(box)
    cv.SetLineColor(ROOT.kRed + 1)
    cv.SetLineWidth(2)
    cv.Draw("L HIST SAME" if smooth else "HIST ][ SAME")
    d.SetMarkerStyle(20)
    d.SetMarkerSize(0.6 if smooth else 0.9)
    d.SetMarkerColor(ROOT.kBlack)
    d.SetLineColor(ROOT.kBlack)
    d.Draw("PE SAME")
    for x in vlines:
        if xlo <= x <= xhi:
            vline(x, 0.0, ymax)
    leg = ROOT.TLegend(0.18, 0.66, 0.55, 0.86) if smooth else ROOT.TLegend(0.45, 0.70, 0.89, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(d, "Data (stat)", "PE")
    leg.AddEntry(band, "MC #pm syst (per-universe fractions)", "LF")
    if vlines:
        leg.AddEntry(_keep[-1] if any(xlo <= x <= xhi for x in vlines) else band,
                     "#nu_{0} in this slice", "L")
    if shade is not None:
        leg.AddEntry(box, "#nu_{0} is the staircase cut here", "F")
    leg.Draw()
    hdr = ROOT.TLatex()
    hdr.SetNDC()
    hdr.SetTextSize(0.028)
    hdr.SetTextAlign(31)
    hdr.DrawLatex(0.895, 0.955, title)
    if label:
        hdr.DrawLatex(0.895, 0.915, label)
    _keep.extend([cv, d, band, leg, hdr])

    bot.cd()
    r, rband = ratio_and_band(d, cv, frac, tag + "_ratio")
    ylo, yhi = ratio_range(r, xlo, xhi)
    style_ratio_axes(rband, ylo, yhi, xtitle)
    rband.GetXaxis().SetRangeUser(xlo, xhi)
    rband.SetTitle("")
    rband.SetFillColorAlpha(ROOT.kRed - 9, 0.6)
    rband.SetLineColor(ROOT.kRed - 9)
    rband.SetMarkerSize(0)
    rband.Draw("E3" if smooth else "E2")
    r.SetMarkerStyle(20)
    r.SetMarkerSize(0.5 if smooth else 0.8)
    r.SetLineColor(ROOT.kBlack)
    r.Draw("PE SAME")
    one = ROOT.TLine(xlo, 1.0, xhi, 1.0)
    one.SetLineStyle(2)
    one.Draw("SAME")
    for x in vlines:
        if xlo <= x <= xhi:
            vline(x, ylo, yhi)
    _keep.extend([r, rband, one])
    c.SaveAs(out_path)
    ROOT.gStyle.SetErrorX(errx)
    print("  wrote {0}".format(out_path))
    return r, rband


def plot_errsum(m, out_path, title, xtitle, xlo, xhi, group=""):
    c = ROOT.TCanvas("c_" + os.path.basename(out_path), "", 1000, 750)
    m.GetXaxis().SetTitle(xtitle)
    m.GetXaxis().SetRangeUser(xlo, xhi)
    try:
        if group:
            plotter.DrawErrorSummary(m, "TR", False, True, BREAKDOWN_THRESHOLD, False,
                                     group, True, "", False)
        else:
            plotter.DrawErrorSummary(m, "TR", True, True, 0.00001, False, "", True, "", True)
    except Exception as e:
        print("  WARN DrawErrorSummary failed ({0})".format(e))
        return
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.03)
    tex.DrawLatex(0.14, 0.92, title + ("   " + group + " knobs" if group else ""))
    _keep.append(tex)
    c.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_ratio_grid(ratios, out_path, xlo, xhi):
    n = len(ratios)
    if n == 0:
        return
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    errx = ROOT.gStyle.GetErrorX()
    ROOT.gStyle.SetErrorX(0.0)
    c = ROOT.TCanvas("c_grid_" + os.path.basename(out_path), "", 1600, 400 * nrow)
    c.Divide(ncol, nrow, 0.002, 0.002)
    for i, (lo, hi, r, rb) in enumerate(ratios):
        c.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetBottomMargin(0.14)
        rb2 = rb.Clone(rb.GetName() + "_g")
        r2 = r.Clone(r.GetName() + "_g")
        rb2.SetTitle("")
        rb2.GetXaxis().SetRangeUser(xlo, xhi)
        rb2.GetYaxis().SetRangeUser(0.7, 1.3)
        rb2.GetXaxis().SetTitle("#nu_{c}  [GeV]")
        rb2.GetYaxis().SetTitle("data/MC of F(#nu < #nu_{c})")
        for ax in (rb2.GetXaxis(), rb2.GetYaxis()):
            ax.SetLabelSize(0.05)
            ax.SetTitleSize(0.055)
        rb2.GetXaxis().SetTitleOffset(1.0)
        rb2.GetYaxis().SetTitleOffset(1.1)
        rb2.Draw("E3")
        r2.SetMarkerSize(0.4)
        r2.Draw("PE SAME")
        one = ROOT.TLine(xlo, 1.0, xhi, 1.0)
        one.SetLineStyle(2)
        one.Draw("SAME")
        for x in nu0_lines(lo, hi):
            if xlo <= x <= xhi:
                vline(x, 0.7, 1.3)
        lab = ROOT.TLatex()
        lab.SetNDC()
        lab.SetTextSize(0.06)
        lab.DrawLatex(0.17, 0.92, "{0:g} < E_{{#mu}} < {1:g} GeV".format(lo, hi))
        _keep.extend([rb2, r2, one, lab])
    c.SaveAs(out_path)
    ROOT.gStyle.SetErrorX(errx)
    print("  wrote {0}".format(out_path))


# --------------------------------------------------------------------------- #
def main():
    a = parse_args()
    for d in (a.outdir, os.path.join(a.outdir, "breakdowns")):
        if not os.path.isdir(d):
            os.makedirs(d)
    print("nucut_efficiency.py v1.0: {0}".format(a.infile))
    f = ROOT.TFile.Open(a.infile)
    if not f or f.IsZombie():
        print("  ERROR: cannot open {0}".format(a.infile))
        sys.exit(2)
    print("  git_describe: {0}".format(stamp(f, "git_describe")))
    label = a.label if a.label is not None else run_label(f)
    if label:
        print("  label: {0}".format(label))

    maps = OrderedDict()
    maps[FINE] = (get_hist(f, "selection_data_{0}_{1}".format(FINE, a.sample)),
                  get_hist(f, "selection_mc_{0}_{1}".format(FINE, a.sample)))
    if not a.no_nucut:
        hd = get_hist(f, "selection_data_{0}_{1}".format(NUCUT, a.sample), required=False)
        hm = get_hist(f, "selection_mc_{0}_{1}".format(NUCUT, a.sample), required=False)
        if hd and hm:
            if same_x_axis(hm, maps[FINE][1]):
                maps[NUCUT] = (hd, hm)
            else:
                print("  WARN: emu_ehad_nucut has a different E_mu axis; using the fine map only")
        else:
            print("  NOTE: no emu_ehad_nucut in this file: thresholds every 0.1 GeV only")
    nv = len(maps[FINE][1].GetVertErrorBandNames())
    nl = len(maps[FINE][1].GetLatErrorBandNames())
    print("  MC bands: {0} vertical, {1} lateral".format(nv, nl))
    if nv + nl == 0:
        print("  WARN: the MC map has no error bands; F has no systematic band.")

    ts = thresholds(maps[FINE][1].GetYaxis(),
                    maps[NUCUT][1].GetYaxis() if NUCUT in maps else None)
    tvals = [t for t, _, _ in ts]
    print("  {0} thresholds, {1:g} to {2:g} GeV ({3})".format(
        len(ts), tvals[0], tvals[-1],
        "0.01 GeV steps to 0.6, 0.05 to 2.4, 0.1 to 5" if NUCUT in maps else "0.1 GeV steps"))
    edges = centred_edges(tvals)

    def index_of(v):
        for i, t in enumerate(tvals):
            if abs(t - v) < EPS:
                return i
        return None

    slices = slice_edges(maps[FINE][1].GetXaxis(), a.emu_edges)
    fout = ROOT.TFile(os.path.join(a.outdir, "nucut_efficiency.root"), "RECREATE")
    per_slice = []
    ratios_full, ratios_zoom = [], []
    worst_closure = None
    xt = "#nu_{c}  [GeV]"
    yt = "F(#nu < #nu_{c})  (total incl. #nu > 5 GeV)"
    for lo, hi, blo, bhi in slices:
        tag = slice_tag(lo, hi)
        cur = slice_curves(maps, blo, bhi, tag, ts)
        if cur is None:
            print("  {0}: empty, skipped".format(tag))
            continue
        if cur["closure"] is not None:
            n, w = cur["closure"]
            if worst_closure is None:
                worst_closure = [n, list(w)]
            else:
                worst_closure[1] = [max(x, y) for x, y in zip(worst_closure[1], w)]
        d = make_data("data_" + tag, edges, cur["fd"], cur["vd"])
        m = make_mc("mc_" + tag, edges, cur["fm"], cur["vm"], cur["univ"])
        fout.cd()
        d.Write()
        m.Write("mc_" + tag)
        per_slice.append((lo, hi, tag, cur, d, m))
        title = "{0:g} < E_{{#mu}} < {1:g} GeV:  fraction of events with #nu < #nu_{{c}}".format(lo, hi)
        lines = nu0_lines(lo, hi)
        r, rb = plot_curve(d, m, os.path.join(a.outdir, "nucut_eff_{0}.png".format(tag)),
                           title, label, 0.0, NU_MAX, xt, yt, lines, ymax=1.05)
        ratios_full.append((lo, hi, r, rb))
        r, rb = plot_curve(d, m, os.path.join(a.outdir, "nucut_eff_{0}_zoom.png".format(tag)),
                           title, label, 0.0, ZOOM_MAX, xt, yt, lines)
        ratios_zoom.append((lo, hi, r, rb))
        if nv + nl > 0:
            stitle = "{0:g} < E_{{#mu}} < {1:g} GeV, F(#nu < #nu_{{c}})".format(lo, hi)
            plot_errsum(m, os.path.join(a.outdir, "nucut_eff_{0}_errsum.png".format(tag)),
                        stitle, xt, 0.0, NU_MAX)
            for g in groups_present(m):
                plot_errsum(m, os.path.join(a.outdir, "breakdowns",
                                            "nucut_eff_{0}_{1}.png".format(tag, g)),
                            stitle, xt, 0.0, NU_MAX, group=g)
    if not per_slice:
        print("  ERROR: every slice is empty")
        sys.exit(2)
    if worst_closure is not None:
        n, w = worst_closure
        flag = "OK" if max(w) < 1e-9 else "WARN"
        print("  closure {0}: fine vs nucut map at {1} shared thresholds, max |dF| "
              "data {2:.2g}, MC CV {3:.2g}, universes {4:.2g}".format(flag, n, w[0], w[1], w[2]))

    plot_ratio_grid(ratios_full, os.path.join(a.outdir, "nucut_eff_ratio_grid.png"), 0.0, NU_MAX)
    plot_ratio_grid(ratios_zoom, os.path.join(a.outdir, "nucut_eff_ratio_grid_zoom.png"),
                    0.0, ZOOM_MAX)

    # F at each nu_0 against E_mu
    eedges = [per_slice[0][0]] + [s[1] for s in per_slice]
    contiguous = all(abs(per_slice[i][1] - per_slice[i + 1][0]) < EPS
                     for i in range(len(per_slice) - 1))
    if not contiguous:
        print("  NOTE: the non-empty slices are not contiguous; the vs-E_mu plots join them")
    for nu0, elo, ehi in STAIRCASE:
        i = index_of(nu0)
        if i is None:
            continue
        ntag = nu0_tag(nu0)
        fd = [s[3]["fd"][i] for s in per_slice]
        vd = [s[3]["vd"][i] for s in per_slice]
        fm = [s[3]["fm"][i] for s in per_slice]
        vm = [s[3]["vm"][i] for s in per_slice]
        univ = OrderedDict()
        for key, (lists, spread) in per_slice[0][3]["univ"].items():
            univ[key] = ([[s[3]["univ"][key][0][iu][i] for s in per_slice]
                          for iu in range(len(lists))], spread)
        d = make_data("data_" + ntag, eedges, fd, vd)
        m = make_mc("mc_" + ntag, eedges, fm, vm, univ)
        fout.cd()
        d.Write()
        m.Write("mc_" + ntag)
        shade = (max(elo - nu0, eedges[0]), min(ehi - nu0, eedges[-1]))
        shade = shade if shade[1] > shade[0] else None
        xe = "E_{#mu}  [GeV]"
        title = "F(#nu < {0:g} GeV) against E_{{#mu}}".format(nu0)
        plot_curve(d, m, os.path.join(a.outdir, "nucut_eff_vs_Emu_{0}.png".format(ntag)),
                   title, label, eedges[0], min(eedges[-1], EMU_PLOT_MAX), xe,
                   "F(#nu < {0:g} GeV)".format(nu0), shade=shade, smooth=False)
        if nv + nl > 0:
            plot_errsum(m, os.path.join(a.outdir, "nucut_eff_vs_Emu_{0}_errsum.png".format(ntag)),
                        title, xe, eedges[0], eedges[-1])
            for g in groups_present(m):
                plot_errsum(m, os.path.join(a.outdir, "breakdowns",
                                            "nucut_eff_vs_Emu_{0}_{1}.png".format(ntag, g)),
                            title, xe, eedges[0], eedges[-1], group=g)
    fout.Close()

    # the table: F at each nu_0 and nu_0 * (1 -+ DELTA)
    targets = []
    for nu0 in NU0:
        for v in (nu0 * (1 - DELTA), nu0, nu0 * (1 + DELTA)):
            if index_of(v) is not None:
                targets.append((nu0, v))
    rows = []
    for lo, hi, tag, cur, d, m in per_slice:
        frac = sys_frac(m)
        for nu0, v in targets:
            i = index_of(v)
            fdv, fmv = cur["fd"][i], cur["fm"][i]
            sd, sm = math.sqrt(cur["vd"][i]), math.sqrt(cur["vm"][i])
            ss = fmv * frac.GetBinContent(i + 1)
            tot = math.sqrt(ss ** 2 + sd ** 2 + sm ** 2)
            pull = (fdv - fmv) / tot if tot > 0 else 0.0
            rows.append((lo, hi, nu0, v, fdv, sd, fmv, sm, ss, fdv / fmv if fmv > 0 else 0.0, pull))
    tpath = os.path.join(a.outdir, "nucut_efficiency_table.txt")
    with open(tpath, "w") as t:
        t.write("# nucut_efficiency.py v1.0: {0}\n# {1}\n".format(a.infile, label))
        t.write("# F = N(nu < nu_c) / N(all) per E_mu slice; MC = all selected MC; "
                "pull = (F_data - F_MC) / sqrt(syst^2 + stat_data^2 + stat_MC^2)\n")
        t.write("{0:>6s} {1:>6s} {2:>5s} {3:>6s} {4:>9s} {5:>9s} {6:>9s} {7:>9s} {8:>9s} "
                "{9:>8s} {10:>7s}\n".format("Emu_lo", "Emu_hi", "nu0", "nu_c", "F_data", "stat_d",
                                            "F_MC", "stat_MC", "syst_MC", "data/MC", "pull"))
        for row in rows:
            t.write("{0:6g} {1:6g} {2:5g} {3:6.3f} {4:9.5f} {5:9.5f} {6:9.5f} {7:9.5f} {8:9.5f} "
                    "{9:8.4f} {10:7.2f}\n".format(*row))
    print("  wrote {0}".format(tpath))

    print("")
    print("  F(nu < nu_0) per slice (MC = all selected MC; syst and stat relative to F):")
    for nu0 in NU0:
        sel = [r for r in rows if r[2] == nu0 and abs(r[3] - nu0) < EPS]
        if not sel:
            continue
        print("  nu_0 = {0:g} GeV".format(nu0))
        print("    {0:>12s} {1:>8s} {2:>8s} {3:>8s} {4:>7s} {5:>7s} {6:>6s}".format(
            "E_mu [GeV]", "F_data", "F_MC", "data/MC", "syst", "stat_d", "pull"))
        for lo, hi, _, _, fdv, sd, fmv, sm, ss, rat, pull in sel:
            print("    {0:>5g} - {1:<5g} {2:8.4f} {3:8.4f} {4:8.3f} {5:6.1%} {6:6.1%} {7:6.1f}".format(
                lo, hi, fdv, fmv, rat, ss / fmv if fmv > 0 else 0.0, sd / fdv if fdv > 0 else 0.0,
                pull))
    n1 = sum(1 for r in rows if abs(r[10]) > 1)
    n2 = sum(1 for r in rows if abs(r[10]) > 2)
    print("  |pull| > 1 at {0} and > 2 at {1} of {2} (slice, threshold) points in the table"
          .format(n1, n2, len(rows)))
    print("  done: {0}".format(a.outdir))


if __name__ == "__main__":
    main()
