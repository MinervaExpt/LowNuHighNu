"""
recoil_emu_slices.py   (v1.0 -- 2026-11-13)

Recoil (nu) distributions in slices of muon energy, from the raw event-loop
output (XSecInputs_*.root).  Kevin's first study (2026-11-06): recoil
projections in bins of E_mu (not E_nu, which mixes muon and recoil effects),
with their systematics.

Input: <prefix>_emu_ehad_fine_<sample>, the (E_mu, nu) map that every event
loop since 2026-10 writes: x = E_mu on the analysis E_nu binning, y = nu on
the 0.1 GeV ehad_fine binning (0-5 GeV); MC with its synced vertical bands,
data stat-only by design.

Per E_mu slice it writes:
  recoil_<slice>.png               data vs MC (POT-normalized), selected
                                   background shown separately, MC systematic
                                   band; ratio panel data/MC with the band
                                   around 1.  Log y.
  recoil_<slice>_areanorm.png      the same, data and MC each normalized to unit
                                   area over 0-5 GeV; the band is shape-only
                                   (the covariance of area-normalized universes).
                                   Linear y.
  recoil_<slice>_errsum.png        grouped fractional systematic breakdown
  recoil_<slice>_errsum_shape.png  the same, shape-only (area-normalized)
and over all slices:
  recoil_ratio_grid.png            data/MC per slice with the systematic band
  recoil_ratio_grid_areanorm.png   the same, area-normalized, shape-only band
  recoil_Emu_slices.root           the projections (data TH1D, MC and
                                   background MnvH1D with their bands, POT-scaled)
                                   for overlays across runs later
  a printed table per slice: events, data/MC, the fraction above 5 GeV.

Dashed lines mark the staircase nu_0 values that are a boundary somewhere in
the slice: an event is low-nu when nu < nu_0(E_nu) with E_nu = E_mu + nu, so in
an E_mu slice the cut can sit at more than one nu_0.

Usage (from plotting/, which holds myPlotStyle.h):
    python3 recoil_emu_slices.py <XSecInputs.root> [<output-dir>]
        [--sample inclusive|lowNu|highNu]   (default inclusive)
        [--emu-edges 1,2,3,...,20]          slice edges in GeV; each must be an
                                            E_mu bin edge of the file
                                            (default: every file bin in 1-20 GeV)
        [--rebin N]                         merge N nu bins (default 1, 0.1 GeV)
        [--label TEXT]                      run label on the plots
                                            (default: the file's tune /
                                            muon_shift / selection stamps)
"""

import argparse
import os
import sys

import ROOT
from functions import *          # brings in PlotUtils (as boundary_diagnostic.py)
from errorMaps import *          # error_bands (OrderedDict)

ROOT.gROOT.SetBatch()
try:
    ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
    ROOT.myPlotStyle()
except Exception:
    pass
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

plotter = PlotUtils.MnvPlotter()
plotter.error_summary_group_map.clear()
for group in error_bands:
    for error in error_bands[group]:
        plotter.error_summary_group_map[group].push_back(error)
plotter.SetLegendNColumns(2)

_keep = []   # Python-side references so ROOT drawables outlive SaveAs

# Staircase, mirroring LowNuBoundary.h: (nu_0, E_nu lo, E_nu hi)
STAIRCASE = [(0.3, 0.0, 3.0), (0.5, 3.0, 7.0), (1.0, 7.0, 12.0), (2.0, 12.0, 1.0e9)]
VAR = "emu_ehad_fine"
NU_MAX = 5.0    # top of the ehad_fine axis (GeV)


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile")
    p.add_argument("outdir", nargs="?", default="./recoil_emu_slices_plots")
    p.add_argument("--sample", default="inclusive", choices=["inclusive", "lowNu", "highNu"])
    p.add_argument("--emu-edges", default=None,
                   help="comma-separated E_mu slice edges in GeV (file bin edges)")
    p.add_argument("--rebin", type=int, default=1)
    p.add_argument("--label", default=None)
    return p.parse_args()


def read_pot(f):
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


def stamp(f, name):
    o = f.Get(name)
    return o.GetTitle() if o else None


def run_label(f):
    parts = []
    for key in ("tune", "muon_shift", "selection"):
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
        print("  ERROR: '{0}' not found. Keys containing '{1}':".format(name, VAR))
        for k in sorted(k.GetName() for k in f.GetListOfKeys() if VAR in k.GetName()):
            print("          {0}".format(k))
        sys.exit(2)
    print("  NOTE: '{0}' not found; the background is not drawn.".format(name))
    return None


def slice_edges(xax, requested):
    """Slices as (lo, hi, bin_lo, bin_hi); every edge must be a file bin edge."""
    file_edges = [xax.GetBinLowEdge(b) for b in range(1, xax.GetNbins() + 2)]
    if requested is None:
        edges = [e for e in file_edges if 1.0 - 1e-6 <= e <= 20.0 + 1e-6]
    else:
        edges = [float(x) for x in requested.split(",")]
        for e in edges:
            if min(abs(e - fe) for fe in file_edges) > 1e-6:
                print("  ERROR: E_mu edge {0:g} is not a bin edge of the file. File edges: {1}"
                      .format(e, ", ".join("{0:g}".format(x) for x in file_edges)))
                sys.exit(2)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append((lo, hi, xax.FindBin(lo + 1e-6), xax.FindBin(hi - 1e-6)))
    return out


def slice_tag(lo, hi):
    if float(lo).is_integer() and float(hi).is_integer():
        return "Emu_{0:02d}to{1:02d}".format(int(lo), int(hi))
    return "Emu_{0:g}to{1:g}".format(lo, hi).replace(".", "p")


def nu0_lines(lo, hi):
    """Staircase nu_0 values that are the cut somewhere in E_mu in [lo, hi]."""
    out = []
    for nu0, elo, ehi in STAIRCASE:
        if lo + nu0 < ehi and hi + nu0 >= elo:
            out.append(nu0)
    return out


# --------------------------------------------------------------------------- #
# Drawing
# --------------------------------------------------------------------------- #
def frac_band(cv, frac):
    """Clone of cv whose errors are cv * frac (frac: fractional TH1)."""
    band = cv.Clone(cv.GetName() + "_band")
    band.SetDirectory(0)
    for b in range(0, band.GetNbinsX() + 2):
        band.SetBinError(b, cv.GetBinContent(b) * frac.GetBinContent(b))
    return band


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


def draw_nu0(lines, ymin, ymax):
    for nu0 in lines:
        ln = ROOT.TLine(nu0, ymin, nu0, ymax)
        ln.SetLineColor(ROOT.kRed + 2)
        ln.SetLineStyle(2)
        ln.SetLineWidth(2)
        ln.Draw("SAME")
        _keep.append(ln)


def ratio_with_band(data, mc_cv, frac, name):
    """(ratio TH1, band TH1 around 1) ; ratio carries data and MC stat errors."""
    r = data.Clone(name)
    r.SetDirectory(0)
    r.Divide(mc_cv)
    band = r.Clone(name + "_band")
    band.SetDirectory(0)
    for b in range(0, band.GetNbinsX() + 2):
        band.SetBinContent(b, 1.0)
        band.SetBinError(b, frac.GetBinContent(b) if frac is not None else 0.0)
    return r, band


def style_ratio_axes(h, ylo, yhi):
    h.GetYaxis().SetRangeUser(ylo, yhi)
    h.GetYaxis().SetTitle("data / MC")
    h.GetYaxis().SetNdivisions(505)
    h.GetYaxis().SetTitleSize(0.11)
    h.GetYaxis().SetTitleOffset(0.55)
    h.GetYaxis().SetLabelSize(0.10)
    h.GetXaxis().SetTitle("#nu = E_{had} (reco)  [GeV]")
    h.GetXaxis().SetTitleSize(0.12)
    h.GetXaxis().SetLabelSize(0.10)
    h.GetXaxis().SetTitleOffset(1.1)


def plot_slice(data, mc_mnv, bg_mnv, lo, hi, out_path, label, area_norm):
    """One slice: overlay (top) and data/MC (bottom).

    data: TH1D (stat); mc_mnv: POT-scaled MnvH1D with bands; bg_mnv or None.
    """
    tag = os.path.splitext(os.path.basename(out_path))[0]
    mc = mc_mnv.GetCVHistoWithStatError()
    mc.SetDirectory(0)
    d = data.Clone(tag + "_data")
    d.SetDirectory(0)
    bg = bg_mnv.GetCVHistoWithStatError() if bg_mnv is not None else None
    nb = mc.GetNbinsX()
    has_bands = len(mc_mnv.GetVertErrorBandNames()) + len(mc_mnv.GetLatErrorBandNames()) > 0
    # sys-only fractional; shape-only when area-normalized
    frac = mc_mnv.GetTotalError(False, True, area_norm) if has_bands else None

    if area_norm:
        mi, di = mc.Integral(1, nb), d.Integral(1, nb)
        if mi > 0:
            if bg is not None:
                bg.Scale(1.0 / mi)
            mc.Scale(1.0 / mi)
        if di > 0:
            d.Scale(1.0 / di)
        ytitle = "fraction of events / {0:g} GeV".format(mc.GetXaxis().GetBinWidth(1))
    else:
        ytitle = "events / {0:g} GeV".format(mc.GetXaxis().GetBinWidth(1))

    band = frac_band(mc, frac) if frac is not None else None
    c, top, bot = two_pads("c_" + tag)

    top.cd()
    if not area_norm and mc.GetMaximum() > 0:
        top.SetLogy()
    frame = band if band is not None else mc
    title = "{0:g} < E_{{#mu}} < {1:g} GeV{2}".format(
        lo, hi, ",  area-normalized (shape-only band)" if area_norm else "")
    frame.SetTitle(title)
    frame.GetYaxis().SetTitle(ytitle)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0)
    hmax = max(mc.GetMaximum(), d.GetMaximum())
    if area_norm or hmax <= 0:
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.3 * hmax if hmax > 0 else 1.0)
    else:
        pos = [mc.GetBinContent(b) for b in range(1, nb + 1) if mc.GetBinContent(b) > 0]
        frame.SetMinimum(0.3 * min(pos))
        frame.SetMaximum(5.0 * hmax)
    if band is not None:
        band.SetFillColorAlpha(ROOT.kRed - 9, 0.6)
        band.SetLineColor(ROOT.kRed + 1)
        band.SetMarkerSize(0)
        band.Draw("E2")
        mc.Draw("HIST SAME")
    else:
        mc.Draw("HIST")
    mc.SetLineColor(ROOT.kRed + 1)
    mc.SetLineWidth(2)
    if bg is not None:
        bg.SetFillColorAlpha(ROOT.kGray + 1, 0.7)
        bg.SetLineColor(ROOT.kGray + 2)
        bg.Draw("HIST SAME")
    d.SetMarkerStyle(20)
    d.SetMarkerSize(0.9)
    d.SetMarkerColor(ROOT.kBlack)
    d.SetLineColor(ROOT.kBlack)
    d.Draw("PE SAME")
    top.Update()
    lines = nu0_lines(lo, hi)
    ylo, yhi = frame.GetMinimum(), frame.GetMaximum()
    draw_nu0(lines, ylo, yhi)

    leg = ROOT.TLegend(0.55, 0.66, 0.89, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(d, "Data (stat)", "PE")
    leg.AddEntry(band if band is not None else mc,
                 "MC #pm syst" if band is not None else "MC (CV)", "LF" if band is not None else "L")
    if bg is not None:
        leg.AddEntry(bg, "MC: selected background", "F")
    if lines:
        leg.AddEntry(_keep[-1], "#nu_{0} in this slice", "L")
    leg.Draw()
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.028)
    hdr = ROOT.TLatex()
    hdr.SetNDC()
    hdr.SetTextSize(0.028)
    hdr.SetTextAlign(31)       # right-aligned: the axis exponent sits top-left
    hdr.DrawLatex(0.895, 0.955, title)
    if label:
        hdr.DrawLatex(0.895, 0.915, label)
    _keep.append(hdr)
    if not area_norm:
        mi, di = mc.Integral(1, nb), d.Integral(1, nb)
        tex.DrawLatex(0.16, 0.07,
                      "data/MC (0-{0:g} GeV) = {1:.3f}".format(NU_MAX, di / mi if mi > 0 else 0.0))
    _keep.extend([mc, d, bg, band, leg, tex])

    bot.cd()
    r, rband = ratio_with_band(d, mc, frac, tag + "_ratio")
    style_ratio_axes(rband, 0.5, 1.5)
    rband.SetTitle("")
    rband.SetFillColorAlpha(ROOT.kRed - 9, 0.6)
    rband.SetLineColor(ROOT.kRed - 9)
    rband.SetMarkerSize(0)
    rband.Draw("E2")
    r.SetMarkerStyle(20)
    r.SetMarkerSize(0.8)
    r.SetLineColor(ROOT.kBlack)
    r.Draw("PE SAME")
    one = ROOT.TLine(r.GetXaxis().GetXmin(), 1.0, r.GetXaxis().GetXmax(), 1.0)
    one.SetLineStyle(2)
    one.Draw("SAME")
    draw_nu0(lines, 0.5, 1.5)
    _keep.extend([r, rband, one])
    c.SaveAs(out_path)
    print("  wrote {0}".format(out_path))
    return r, rband


def plot_errsum(mc_mnv, lo, hi, out_path, area_norm):
    c = ROOT.TCanvas("c_" + os.path.basename(out_path), "", 1000, 750)
    try:
        mc_mnv.GetXaxis().SetTitle("#nu = E_{had} (reco)  [GeV]")
        plotter.DrawErrorSummary(mc_mnv, "TR", True, True, 0.00001, area_norm,
                                 "", True, "", True)
    except Exception as e:
        print("  WARN DrawErrorSummary failed ({0})".format(e))
        return
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.03)
    tex.DrawLatex(0.14, 0.92, "{0:g} < E_{{#mu}} < {1:g} GeV{2}".format(
        lo, hi, ", shape-only" if area_norm else ""))
    _keep.append(tex)
    c.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_ratio_grid(ratios, out_path, title):
    """All slices' data/MC panels on one page (the band: sys, or shape-only)."""
    n = len(ratios)
    if n == 0:
        return
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    c = ROOT.TCanvas("c_grid_" + os.path.basename(out_path), "", 1600, 400 * nrow)
    c.Divide(ncol, nrow, 0.002, 0.002)
    for i, (lo, hi, r, rb) in enumerate(ratios):
        c.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetBottomMargin(0.14)
        rb2 = rb.Clone(rb.GetName() + "_g")
        r2 = r.Clone(r.GetName() + "_g")
        rb2.SetTitle("{0:g} < E_{{#mu}} < {1:g} GeV".format(lo, hi))
        rb2.GetYaxis().SetRangeUser(0.5, 1.5)
        for ax in (rb2.GetXaxis(), rb2.GetYaxis()):
            ax.SetLabelSize(0.05)
            ax.SetTitleSize(0.055)
        rb2.GetXaxis().SetTitleOffset(1.0)
        rb2.GetYaxis().SetTitleOffset(1.1)
        rb2.Draw("E2")
        r2.Draw("PE SAME")
        one = ROOT.TLine(r2.GetXaxis().GetXmin(), 1.0, r2.GetXaxis().GetXmax(), 1.0)
        one.SetLineStyle(2)
        one.Draw("SAME")
        draw_nu0(nu0_lines(lo, hi), 0.5, 1.5)
        lab = ROOT.TLatex()
        lab.SetNDC()
        lab.SetTextSize(0.06)
        lab.DrawLatex(0.17, 0.92, "{0:g} < E_{{#mu}} < {1:g} GeV".format(lo, hi))
        _keep.extend([rb2, r2, one, lab])
    c.SaveAs(out_path)
    print("  wrote {0}  ({1})".format(out_path, title))


# --------------------------------------------------------------------------- #
def main():
    a = parse_args()
    if not os.path.isdir(a.outdir):
        os.makedirs(a.outdir)
    print("recoil_emu_slices.py v1.0: {0}".format(a.infile))
    f = ROOT.TFile.Open(a.infile)
    if not f or f.IsZombie():
        print("  ERROR: cannot open {0}".format(a.infile))
        sys.exit(2)
    print("  git_describe: {0}".format(stamp(f, "git_describe")))
    label = a.label if a.label is not None else run_label(f)
    if label:
        print("  label: {0}".format(label))
    mc_pot, data_pot = read_pot(f)
    mc_scale = data_pot / mc_pot

    h_data = get_hist(f, "selection_data_{0}_{1}".format(VAR, a.sample))
    h_mc = get_hist(f, "selection_mc_{0}_{1}".format(VAR, a.sample))
    h_bg = get_hist(f, "bg_{0}_{1}".format(VAR, a.sample), required=False)
    try:
        nv = len(h_mc.GetVertErrorBandNames())
        nl = len(h_mc.GetLatErrorBandNames())
    except Exception:
        nv, nl = 0, 0
    print("  MC bands: {0} vertical, {1} lateral".format(nv, nl))
    if nv + nl == 0:
        print("  WARN: the MC map has no error bands; plots are CV-only.")

    slices = slice_edges(h_mc.GetXaxis(), a.emu_edges)
    fout = ROOT.TFile(os.path.join(a.outdir, "recoil_Emu_slices.root"), "RECREATE")
    ratios_abs, ratios_shape = [], []
    rows = []
    for lo, hi, blo, bhi in slices:
        tag = slice_tag(lo, hi)
        d = h_data.ProjectionY("data_" + tag, blo, bhi)
        if hasattr(d, "GetCVHistoWithStatError"):
            d = d.GetCVHistoWithStatError()
        d.SetName("data_" + tag)
        d.SetDirectory(0)
        m = h_mc.ProjectionY("mc_" + tag, blo, bhi)      # MnvH1D, bands kept
        m.Scale(mc_scale)
        b = None
        if h_bg is not None:
            b = h_bg.ProjectionY("bg_" + tag, blo, bhi)
            b.Scale(mc_scale)
        if a.rebin > 1:
            d.Rebin(a.rebin)
            m.Rebin(a.rebin)          # in place; MnvH1D::Rebin carries the bands
            if b is not None:
                b.Rebin(a.rebin)
        nb = m.GetNbinsX()
        dn, mn = d.Integral(1, nb), m.Integral(1, nb)
        d_of = d.GetBinContent(nb + 1) / (dn + d.GetBinContent(nb + 1)) if dn > 0 else 0.0
        m_of = m.GetBinContent(nb + 1) / (mn + m.GetBinContent(nb + 1)) if mn > 0 else 0.0
        bgf = b.Integral(1, nb) / mn if (b is not None and mn > 0) else float("nan")
        rows.append((lo, hi, dn, mn, dn / mn if mn > 0 else 0.0, bgf, d_of, m_of))
        if dn <= 0 and mn <= 0:
            print("  {0}: empty, skipped".format(tag))
            continue
        fout.cd()
        d.Write()
        m.Write("mc_" + tag)
        if b is not None:
            b.Write("bg_" + tag)

        r, rb = plot_slice(d, m, b, lo, hi, os.path.join(a.outdir, "recoil_" + tag + ".png"),
                           label, area_norm=False)
        ratios_abs.append((lo, hi, r, rb))
        r, rb = plot_slice(d, m, b, lo, hi,
                           os.path.join(a.outdir, "recoil_" + tag + "_areanorm.png"),
                           label, area_norm=True)
        ratios_shape.append((lo, hi, r, rb))
        if nv + nl > 0:
            plot_errsum(m, lo, hi, os.path.join(a.outdir, "recoil_" + tag + "_errsum.png"), False)
            plot_errsum(m, lo, hi, os.path.join(a.outdir, "recoil_" + tag + "_errsum_shape.png"),
                        True)
        _keep.extend([d, m, b])
    fout.Close()
    plot_ratio_grid(ratios_abs, os.path.join(a.outdir, "recoil_ratio_grid.png"),
                    "data/MC, POT-normalized")
    plot_ratio_grid(ratios_shape, os.path.join(a.outdir, "recoil_ratio_grid_areanorm.png"),
                    "data/MC, area-normalized")

    print("")
    print("  Per slice (nu in 0-{0:g} GeV; MC POT-scaled):".format(NU_MAX))
    print("    {0:>12s} {1:>11s} {2:>11s} {3:>8s} {4:>8s} {5:>12s} {6:>12s}".format(
        "E_mu [GeV]", "data", "MC", "data/MC", "bg/MC", "data >5 GeV", "MC >5 GeV"))
    for lo, hi, dn, mn, r, bgf, dof, mof in rows:
        print("    {0:>5g} - {1:<5g} {2:>11.1f} {3:>11.1f} {4:>8.3f} {5:>8.3f} {6:>11.2%} {7:>11.2%}"
              .format(lo, hi, dn, mn, r, bgf, dof, mof))
    print("  done: {0}".format(a.outdir))


if __name__ == "__main__":
    main()
