"""
cv_diagnostics_plots.py   (v1.0 -- 2026-11-20)

Plots of the central-value diagnostics of patch 0035: what the selected MC is
made of.  Two parts, from the raw event-loop output (XSecInputs_*.root), merged:

(1) What fills the lowest recoil (Kevin, 2026-11-13: "are we doing the right
    thing at the very low nu?").  The selected MC by GENIE channel on the
    (E_mu, nu) map, diag_chan_{qe,2p2h,res,dis,coh,other,bg}_emu_ehad_fine
    (signal by interaction type; all background in one), with the data map
    selection_data_emu_ehad_fine_inclusive.  Per E_mu slice (Emu_<lo>to<hi>):
      chan_<slice>.png          the recoil distribution stacked by channel,
                                POT-normalized, data on top; data/MC below.
                                0-5 GeV and the overflow (nu > 5) as a last
                                bin; log y
      chan_<slice>_zoom.png     the same for nu < 2 GeV, linear y
      chan_<slice>_frac.png     each channel's share of the selected MC
                                against nu (stacked to 1)
    and over all slices:
      chan_stack_grid_zoom.png  the zoomed stacks, one pad per slice
      chan_frac_grid_zoom.png   the shares for nu < 2 GeV, one pad per slice
      chan_table.txt            per slice and nu range: data, MC, data/MC and
                                each channel's share
    The log prints, per slice, data/MC and the shares for nu < 0.3 GeV and
    0.3-1 GeV.

(2) The selected background by true category, for the background-subtraction
    page: diag_bg_{nc,numubar,nue,ccnumu,other}_enu_ehad_<sample> (neutral
    current; CC anti-nu_mu, the wrong sign; CC nu_e or anti-nu_e; CC nu_mu
    outside the signal volume; anything else), against reco E_nu, as a share
    of all selected MC (selection_mc_enu_ehad_<sample>), for the low-nu and
    inclusive samples:
      bg_frac_vs_Enu_<sample>.png   each category's share and the total
      bg_table.txt                  per E_nu region (0-3, 3-7, 7-12, 12-20,
                                    above 20 GeV): data, selected MC, the
                                    background share by category
    printed in the log as well.

Both parts are central values only: these histograms carry no error bands.
Every count includes the under- and overflow of the axis it is summed over.
The MC is scaled to the data POT.

Usage (GPVM, from plotting/, which holds myPlotStyle.h; py_classes/ comes from
setup_mat.sh):
    python3 cv_diagnostics_plots.py <XSecInputs.root> [<output-dir>]
        [--emu-edges 1,2,3,...,20]   slice edges in GeV, each an E_mu bin edge
                                     (default: every file bin in 1-20 GeV)
        [--label TEXT]               run label on the plots (default: the
                                     file's stamps)

Drafted 2026-11-20 (card A15's channel breakdown; card A11's input).
Co-authored by Rob Fine and Claude (Anthropic).
"""

import argparse
from array import array
import os
import sys

import ROOT
from functions import *          # brings in PlotUtils (as recoil_emu_slices.py)

ROOT.gROOT.SetBatch()
try:
    ROOT.gROOT.ProcessLine(".L myPlotStyle.h")
    ROOT.myPlotStyle()
except Exception:
    pass
ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

_keep = []

CHANNELS = [("qe", "QE", ROOT.kAzure + 1), ("2p2h", "2p2h", ROOT.kSpring - 5),
            ("res", "Resonant", ROOT.kOrange - 3), ("dis", "DIS", ROOT.kViolet - 4),
            ("coh", "Coherent", ROOT.kPink + 1), ("other", "Other signal", ROOT.kGray),
            ("bg", "Background", ROOT.kGray + 2)]
BG_CATS = [("nc", "NC", ROOT.kAzure + 1), ("numubar", "CC #bar{#nu}_{#mu} (wrong sign)", ROOT.kOrange - 3),
           ("nue", "CC #nu_{e} / #bar{#nu}_{e}", ROOT.kSpring - 5),
           ("ccnumu", "CC #nu_{#mu} outside the volume", ROOT.kViolet - 4),
           ("other", "other", ROOT.kGray + 1)]
STAIRCASE = [(0.3, 0.0, 3.0), (0.5, 3.0, 7.0), (1.0, 7.0, 12.0), (2.0, 12.0, 1.0e9)]
NU_MAX = 5.0
OVERFLOW_WIDTH = 0.5
ZOOM_MAX = 2.0
NU_RANGES = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5), (0.5, 1.0), (1.0, 2.0),
             (2.0, 5.0), (5.0, None)]
ENU_REGIONS = [(0.0, 3.0), (3.0, 7.0), (7.0, 12.0), (12.0, 20.0), (20.0, None)]
ENU_PLOT_MAX = 20.0
EPS = 1e-6


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile")
    p.add_argument("outdir", nargs="?", default="./cv_diagnostics_plots")
    p.add_argument("--emu-edges", default=None)
    p.add_argument("--label", default=None)
    return p.parse_args()


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


def get(f, name):
    h = f.Get(name)
    if not h:
        print("  ERROR: '{0}' not found (patch 0035 writes it; older files have none)".format(name))
        sys.exit(2)
    return h


def cv(h):
    c = h.GetCVHistoWithStatError() if hasattr(h, "GetCVHistoWithStatError") else h.Clone()
    c.SetDirectory(0)
    return c


def slice_edges(xax, requested):
    file_edges = [xax.GetBinLowEdge(b) for b in range(1, xax.GetNbins() + 2)]
    if requested is None:
        edges = [e for e in file_edges if 1.0 - EPS <= e <= 20.0 + EPS]
    else:
        edges = [float(x) for x in requested.split(",")]
        for e in edges:
            if min(abs(e - fe) for fe in file_edges) > EPS:
                print("  ERROR: E_mu edge {0:g} is not a bin edge of the file".format(e))
                sys.exit(2)
    return [(lo, hi, xax.FindBin(lo + EPS), xax.FindBin(hi - EPS))
            for lo, hi in zip(edges[:-1], edges[1:])]


def slice_tag(lo, hi):
    if float(lo).is_integer() and float(hi).is_integer():
        return "Emu_{0:02d}to{1:02d}".format(int(lo), int(hi))
    return "Emu_{0:g}to{1:g}".format(lo, hi).replace(".", "p")


def nu0_lines(lo, hi):
    return [nu0 for nu0, elo, ehi in STAIRCASE if lo + nu0 < ehi and hi + nu0 > elo]


def with_overflow_bin(h, name):
    """TH1D copy with one more visible bin holding the overflow (nu > NU_MAX)."""
    ax = h.GetXaxis()
    nb = h.GetNbinsX()
    edges = [ax.GetBinLowEdge(b) for b in range(1, nb + 2)] + [ax.GetBinUpEdge(nb) + OVERFLOW_WIDTH]
    out = ROOT.TH1D(name, "", len(edges) - 1, array('d', edges))
    out.SetDirectory(0)
    for b in range(0, nb + 2):
        out.SetBinContent(b, h.GetBinContent(b))
        out.SetBinError(b, h.GetBinError(b))
    out.SetBinContent(nb + 2, 0.0)
    return out


def range_sum(h, lo, hi):
    """Sum of h over [lo, hi) (hi None: to the overflow), under/overflow included."""
    ax = h.GetXaxis()
    b0 = 0 if lo <= ax.GetXmin() + EPS else ax.FindBin(lo + EPS)
    b1 = h.GetNbinsX() + 1 if hi is None else ax.FindBin(hi - EPS)
    return sum(h.GetBinContent(b) for b in range(b0, b1 + 1))


def vline(x, ymin, ymax):
    ln = ROOT.TLine(x, ymin, x, ymax)
    ln.SetLineColor(ROOT.kRed + 2)
    ln.SetLineStyle(2)
    ln.SetLineWidth(2)
    ln.Draw("SAME")
    _keep.append(ln)
    return ln


def overflow_marker(ymin, ymax, size=0.03):
    ln = ROOT.TLine(NU_MAX, ymin, NU_MAX, ymax)
    ln.SetLineColor(ROOT.kGray + 2)
    ln.SetLineStyle(3)
    ln.Draw("SAME")
    t = ROOT.TLatex()
    t.SetTextSize(size)
    t.SetTextColor(ROOT.kGray + 2)
    t.SetTextAngle(90)
    t.DrawLatex(NU_MAX + 0.35 * OVERFLOW_WIDTH, ymin + 0.05 * (ymax - ymin), "#nu > 5")
    _keep.extend([ln, t])


def header(title, label, size=0.028):
    hdr = ROOT.TLatex()
    hdr.SetNDC()
    hdr.SetTextSize(size)
    hdr.SetTextAlign(31)
    hdr.DrawLatex(0.895, 0.955, title)
    if label:
        hdr.DrawLatex(0.895, 0.915, label)
    _keep.append(hdr)


# --------------------------------------------------------------------------- #
# (1) channels
# --------------------------------------------------------------------------- #
def stack_of(chans, name):
    """Cumulative histograms: element i is the sum of channels 0..i.

    Drawn from the last to the first, they make the stack (no THStack: simpler
    to keep alive from Python, and it draws the same everywhere)."""
    cums = []
    run = None
    for (key, title, colour), h in zip(CHANNELS, chans):
        c = h.Clone("{0}_{1}".format(name, key))
        c.SetDirectory(0)
        if run is not None:
            c.Add(run)
        c.SetFillColor(colour)
        c.SetFillStyle(1001)
        c.SetLineColor(colour)
        c.SetLineWidth(1)
        cums.append(c)
        run = c
    _keep.extend(cums)
    return cums


def draw_stack(cums):
    for c in reversed(cums):
        c.Draw("HIST SAME")


def plot_stack(data, chans, total, lo, hi, out_path, label, zoom):
    tag = os.path.splitext(os.path.basename(out_path))[0]
    c = ROOT.TCanvas("c_" + tag, "", 900, 900)
    top = ROOT.TPad(tag + "_t", "", 0.0, 0.30, 1.0, 1.0)
    bot = ROOT.TPad(tag + "_b", "", 0.0, 0.0, 1.0, 0.30)
    top.SetBottomMargin(0.02)
    top.SetTopMargin(0.10)
    top.SetLeftMargin(0.16)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.32)
    bot.SetLeftMargin(0.16)
    top.Draw()
    bot.Draw()
    _keep.extend([c, top, bot])
    xmax = ZOOM_MAX if zoom else NU_MAX + OVERFLOW_WIDTH
    top.cd()
    st = stack_of(chans, "st_" + tag)
    frame = total.Clone("frame_" + tag)
    frame.Reset()
    frame.SetDirectory(0)
    frame.GetXaxis().SetRangeUser(0.0, xmax - EPS)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetYaxis().SetTitle("events / {0:g} GeV".format(total.GetXaxis().GetBinWidth(1)))
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.GetYaxis().SetTitleSize(0.045)
    vis = [max(total.GetBinContent(b), data.GetBinContent(b)) for b in range(1, total.GetNbinsX() + 1)
           if total.GetBinCenter(b) < xmax]
    ymax = max(vis) if vis else 1.0
    if zoom:
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.35 * ymax)
    else:
        top.SetLogy()
        pos = [total.GetBinContent(b) for b in range(1, total.GetNbinsX() + 1) if total.GetBinContent(b) > 0]
        frame.SetMinimum(0.3 * min(pos) if pos else 0.1)
        frame.SetMaximum(8.0 * ymax)
    frame.Draw("AXIS")
    draw_stack(st)
    d = data.Clone("d_" + tag)
    d.SetMarkerStyle(20)
    d.SetMarkerSize(0.8)
    d.SetLineColor(ROOT.kBlack)
    d.Draw("PE SAME")
    frame.Draw("AXIS SAME")
    ylo, yhi = frame.GetMinimum(), frame.GetMaximum()
    if not zoom:
        overflow_marker(ylo, yhi)
    lines = nu0_lines(lo, hi)
    for x in lines:
        if x < xmax:
            vline(x, ylo, yhi)
    leg = ROOT.TLegend(0.55, 0.50, 0.89, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(d, "Data (stat)", "PE")
    for (key, title, colour), h in reversed(list(zip(CHANNELS, st))):
        leg.AddEntry(h, title, "F")
    if lines:
        leg.AddEntry(_keep[-1], "#nu_{0} in this slice", "L")
    leg.Draw()
    header("{0:g} < E_{{#mu}} < {1:g} GeV:  selected MC by channel (CV)".format(lo, hi), label)
    _keep.extend([frame, d, leg])

    bot.cd()
    r = data.Clone("r_" + tag)
    r.SetDirectory(0)
    r.Divide(total)
    r.SetTitle("")
    r.GetXaxis().SetRangeUser(0.0, xmax - EPS)
    r.GetYaxis().SetRangeUser(0.6, 1.4)
    r.GetYaxis().SetTitle("data / MC")
    r.GetYaxis().SetNdivisions(505)
    r.GetYaxis().SetTitleSize(0.11)
    r.GetYaxis().SetTitleOffset(0.55)
    r.GetYaxis().SetLabelSize(0.10)
    r.GetXaxis().SetTitle("#nu = E_{had} (reco)  [GeV]" + ("" if zoom else "  (last bin: #nu > 5)"))
    r.GetXaxis().SetTitleSize(0.12)
    r.GetXaxis().SetLabelSize(0.10)
    r.GetXaxis().SetTitleOffset(1.1)
    r.SetMarkerStyle(20)
    r.SetMarkerSize(0.8)
    r.SetLineColor(ROOT.kBlack)
    r.Draw("PE")
    one = ROOT.TLine(0.0, 1.0, xmax, 1.0)
    one.SetLineStyle(2)
    one.Draw("SAME")
    for x in lines:
        if x < xmax:
            vline(x, 0.6, 1.4)
    _keep.extend([r, one])
    c.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def fractions(chans, total, name):
    out = []
    for (key, title, colour), h in zip(CHANNELS, chans):
        fr = h.Clone("{0}_{1}".format(name, key))
        fr.SetDirectory(0)
        for b in range(0, fr.GetNbinsX() + 2):
            t = total.GetBinContent(b)
            fr.SetBinContent(b, h.GetBinContent(b) / t if t > 0 else 0.0)
            fr.SetBinError(b, 0.0)
        out.append(fr)
    return out


def draw_frac(fr, lo, hi, xmax, small=False):
    st = stack_of(fr, "stf_{0:g}_{1:g}_{2}".format(lo, hi, len(_keep)))
    frame = fr[0].Clone("ff_{0}".format(len(_keep)))
    frame.Reset()
    frame.SetDirectory(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.0)
    frame.GetXaxis().SetRangeUser(0.0, xmax - EPS)
    frame.GetXaxis().SetTitle("#nu = E_{had} (reco)  [GeV]")
    frame.GetYaxis().SetTitle("share of selected MC")
    if small:
        for ax in (frame.GetXaxis(), frame.GetYaxis()):
            ax.SetLabelSize(0.05)
            ax.SetTitleSize(0.055)
    frame.Draw("AXIS")
    draw_stack(st)
    frame.Draw("AXIS SAME")
    for x in nu0_lines(lo, hi):
        if x < xmax:
            vline(x, 0.0, 1.0)
    _keep.append(frame)
    return st


def plot_frac(fr, lo, hi, out_path, label):
    c = ROOT.TCanvas("c_" + os.path.basename(out_path), "", 1000, 750)
    ROOT.gPad.SetTopMargin(0.12)
    ROOT.gPad.SetRightMargin(0.27)
    st = draw_frac(fr, lo, hi, NU_MAX + OVERFLOW_WIDTH)
    overflow_marker(0.0, 1.0)
    leg = ROOT.TLegend(0.74, 0.35, 0.99, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for (key, title, colour), h in reversed(list(zip(CHANNELS, st))):
        leg.AddEntry(h, title, "F")
    leg.Draw()
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.03)
    tex.DrawLatex(0.12, 0.95, "{0:g} < E_{{#mu}} < {1:g} GeV:  share of the selected MC by channel".format(lo, hi))
    if label:
        tex.DrawLatex(0.12, 0.91, label)
    _keep.extend([leg, tex])
    c.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def plot_grid(items, out_path, kind):
    n = len(items)
    if n == 0:
        return
    ncol = 4
    nrow = (n + ncol - 1) // ncol
    c = ROOT.TCanvas("c_" + os.path.basename(out_path), "", 1600, 400 * nrow)
    c.Divide(ncol, nrow, 0.002, 0.002)
    for i, it in enumerate(items):
        c.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetBottomMargin(0.14)
        lo, hi = it[0], it[1]
        if kind == "frac":
            draw_frac(it[2], lo, hi, ZOOM_MAX, small=True)
        else:
            data, chans, total = it[2], it[3], it[4]
            st = stack_of(chans, "stg_{0}".format(i))
            frame = total.Clone("fg_{0}".format(i))
            frame.Reset()
            frame.SetDirectory(0)
            frame.GetXaxis().SetRangeUser(0.0, ZOOM_MAX - EPS)
            vis = [max(total.GetBinContent(b), data.GetBinContent(b))
                   for b in range(1, total.GetNbinsX() + 1) if total.GetBinCenter(b) < ZOOM_MAX]
            frame.SetMinimum(0.0)
            frame.SetMaximum(1.3 * max(vis) if vis else 1.0)
            frame.GetXaxis().SetTitle("#nu  [GeV]")
            for ax in (frame.GetXaxis(), frame.GetYaxis()):
                ax.SetLabelSize(0.05)
                ax.SetTitleSize(0.055)
            frame.Draw("AXIS")
            draw_stack(st)
            d = data.Clone("dg_{0}".format(i))
            d.SetMarkerStyle(20)
            d.SetMarkerSize(0.5)
            d.SetMarkerColor(ROOT.kBlack)
            d.SetLineColor(ROOT.kBlack)
            d.Draw("PE SAME")
            frame.Draw("AXIS SAME")
            for x in nu0_lines(lo, hi):
                if x < ZOOM_MAX:
                    vline(x, 0.0, frame.GetMaximum())
            _keep.extend([frame, d])
        lab = ROOT.TLatex()
        lab.SetNDC()
        lab.SetTextSize(0.06)
        lab.DrawLatex(0.17, 0.92, "{0:g} < E_{{#mu}} < {1:g} GeV".format(lo, hi))
        _keep.append(lab)
    c.SaveAs(out_path)
    print("  wrote {0}".format(out_path))


def channels_part(f, a, label, mc_scale):
    print("  (1) channels on (E_mu, nu)")
    h_data = get(f, "selection_data_emu_ehad_fine_inclusive")
    h_chan = [get(f, "diag_chan_{0}_emu_ehad_fine".format(k)) for k, _, _ in CHANNELS]
    rows = []
    stacks, fracs = [], []
    for lo, hi, blo, bhi in slice_edges(h_chan[0].GetXaxis(), a.emu_edges):
        tag = slice_tag(lo, hi)
        d = cv(h_data.ProjectionY("d_" + tag, blo, bhi))
        chans = []
        for (k, _, _), h in zip(CHANNELS, h_chan):
            p = cv(h.ProjectionY("c_{0}_{1}".format(k, tag), blo, bhi))
            p.Scale(mc_scale)
            chans.append(p)
        total = chans[0].Clone("tot_" + tag)
        total.SetDirectory(0)
        for p in chans[1:]:
            total.Add(p)
        if total.Integral(0, total.GetNbinsX() + 1) <= 0 and d.Integral(0, d.GetNbinsX() + 1) <= 0:
            print("  {0}: empty, skipped".format(tag))
            continue
        # rows for the table, from the raw projections (under/overflow included)
        for nlo, nhi in NU_RANGES:
            dv, tv = range_sum(d, nlo, nhi), range_sum(total, nlo, nhi)
            shares = [range_sum(p, nlo, nhi) / tv if tv > 0 else 0.0 for p in chans]
            rows.append((lo, hi, nlo, nhi, dv, tv, dv / tv if tv > 0 else 0.0, shares))
        d = with_overflow_bin(d, "dov_" + tag)
        chans = [with_overflow_bin(p, p.GetName() + "_ov") for p in chans]
        total = with_overflow_bin(total, "totov_" + tag)
        fr = fractions(chans, total, "fr_" + tag)
        plot_stack(d, chans, total, lo, hi, os.path.join(a.outdir, "chan_{0}.png".format(tag)),
                   label, zoom=False)
        plot_stack(d, chans, total, lo, hi, os.path.join(a.outdir, "chan_{0}_zoom.png".format(tag)),
                   label, zoom=True)
        plot_frac(fr, lo, hi, os.path.join(a.outdir, "chan_{0}_frac.png".format(tag)), label)
        stacks.append((lo, hi, d, chans, total))
        fracs.append((lo, hi, fr))
    plot_grid(stacks, os.path.join(a.outdir, "chan_stack_grid_zoom.png"), "stack")
    plot_grid(fracs, os.path.join(a.outdir, "chan_frac_grid_zoom.png"), "frac")

    tpath = os.path.join(a.outdir, "chan_table.txt")
    with open(tpath, "w") as t:
        t.write("# cv_diagnostics_plots.py v1.0: {0}\n# {1}\n".format(a.infile, label))
        t.write("# selected MC by channel, MC scaled to the data POT; shares of the selected MC\n")
        t.write("{0:>6s} {1:>6s} {2:>5s} {3:>5s} {4:>11s} {5:>11s} {6:>9s}".format(
            "Emu_lo", "Emu_hi", "nu_lo", "nu_hi", "data", "MC", "data/MC"))
        for k, _, _ in CHANNELS:
            t.write(" {0:>9s}".format(k))
        t.write("\n")
        for lo, hi, nlo, nhi, dv, tv, r, shares in rows:
            t.write("{0:6g} {1:6g} {2:5g} {3:>5s} {4:11.1f} {5:11.2f} {6:9.6f}".format(
                lo, hi, nlo, "inf" if nhi is None else "{0:g}".format(nhi), dv, tv, r))
            for s in shares:
                t.write(" {0:9.6f}".format(s))
            t.write("\n")
    print("  wrote {0}".format(tpath))

    print("")
    print("  Shares of the selected MC by channel, and data/MC (CV, POT-normalized):")
    hdr = "    {0:>12s} {1:>9s} {2:>8s}".format("E_mu [GeV]", "nu [GeV]", "data/MC")
    for k, _, _ in CHANNELS:
        hdr += " {0:>6s}".format(k)
    print(hdr)
    for lo, hi in sorted(set((r[0], r[1]) for r in rows)):
        for nlo, nhi in ((0.0, 0.3), (0.3, 1.0)):
            sel = [r for r in rows if r[0] == lo and r[1] == hi and nlo - EPS <= r[2] and
                   r[3] is not None and r[3] <= nhi + EPS]
            dv = sum(r[4] for r in sel)
            tv = sum(r[5] for r in sel)
            shares = [sum(r[5] * r[7][i] for r in sel) / tv if tv > 0 else 0.0
                      for i in range(len(CHANNELS))]
            line = "    {0:>5g} - {1:<5g} {2:>4g}-{3:<4g} {4:8.3f}".format(lo, hi, nlo, nhi,
                                                                     dv / tv if tv > 0 else 0.0)
            for s in shares:
                line += " {0:6.1%}".format(s)
            print(line)


# --------------------------------------------------------------------------- #
# (2) background categories
# --------------------------------------------------------------------------- #
def background_part(f, a, label, mc_scale):
    print("")
    print("  (2) background by category against reco E_nu")
    rows = []
    for sample in ("lowNu", "inclusive"):
        h_mc = get(f, "selection_mc_enu_ehad_{0}".format(sample))
        h_data = get(f, "selection_data_enu_ehad_{0}".format(sample))
        ny = h_mc.GetNbinsY() + 1
        mc = cv(h_mc.ProjectionX("bmc_" + sample, 0, ny))
        mc.Scale(mc_scale)
        data = cv(h_data.ProjectionX("bdata_" + sample, 0, ny))
        cats = []
        for k, _, _ in BG_CATS:
            p = cv(get(f, "diag_bg_{0}_enu_ehad_{1}".format(k, sample)).ProjectionX(
                "bg_{0}_{1}".format(k, sample), 0, ny))
            p.Scale(mc_scale)
            cats.append(p)
        tot = cats[0].Clone("bgtot_" + sample)
        tot.SetDirectory(0)
        for p in cats[1:]:
            tot.Add(p)
        for elo, ehi in ENU_REGIONS:
            mv, dv = range_sum(mc, elo, ehi), range_sum(data, elo, ehi)
            bv = range_sum(tot, elo, ehi)
            rows.append((sample, elo, ehi, dv, mv, bv / mv if mv > 0 else 0.0,
                         [range_sum(p, elo, ehi) / mv if mv > 0 else 0.0 for p in cats]))
        # the plot: shares against E_nu
        c = ROOT.TCanvas("c_bg_" + sample, "", 1000, 750)
        ROOT.gPad.SetTopMargin(0.12)
        ROOT.gPad.SetLeftMargin(0.13)
        frs = []
        for (k, title, colour), p in zip(BG_CATS + [("total", "all background", ROOT.kBlack)], cats + [tot]):
            fr = p.Clone("bgfr_{0}_{1}".format(k, sample))
            fr.SetDirectory(0)
            for b in range(0, fr.GetNbinsX() + 2):
                m = mc.GetBinContent(b)
                fr.SetBinContent(b, 100.0 * p.GetBinContent(b) / m if m > 0 else 0.0)
                fr.SetBinError(b, 0.0)
            fr.SetLineColor(colour)
            fr.SetLineWidth(3 if k == "total" else 2)
            frs.append((fr, title))
        ymax = max(fr.GetBinContent(b) for fr, _ in frs for b in range(1, fr.GetNbinsX() + 1)
                   if fr.GetBinCenter(b) < ENU_PLOT_MAX)
        frame = frs[-1][0]
        frame.SetTitle("")
        frame.GetXaxis().SetRangeUser(0.0, ENU_PLOT_MAX - EPS)
        frame.GetXaxis().SetTitle("E_{#nu} (reco)  [GeV]")
        frame.GetYaxis().SetTitle("share of the selected MC  [%]")
        frame.SetMinimum(0.0)
        frame.SetMaximum(1.3 * ymax if ymax > 0 else 1.0)
        frame.Draw("HIST ][")
        for fr, _ in frs[:-1]:
            fr.Draw("HIST ][ SAME")
        frame.Draw("HIST ][ SAME")
        leg = ROOT.TLegend(0.45, 0.62, 0.89, 0.87)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for fr, title in reversed(frs):
            leg.AddEntry(fr, title, "L")
        leg.Draw()
        tex = ROOT.TLatex()
        tex.SetNDC()
        tex.SetTextSize(0.03)
        tex.DrawLatex(0.13, 0.95, "selected background by true category, {0} sample (CV)".format(sample))
        if label:
            tex.DrawLatex(0.13, 0.91, label)
        _keep.extend([c, leg, tex] + [fr for fr, _ in frs])
        out_path = os.path.join(a.outdir, "bg_frac_vs_Enu_{0}.png".format(sample))
        c.SaveAs(out_path)
        print("  wrote {0}".format(out_path))

    tpath = os.path.join(a.outdir, "bg_table.txt")
    with open(tpath, "w") as t:
        t.write("# cv_diagnostics_plots.py v1.0: {0}\n# {1}\n".format(a.infile, label))
        t.write("# reco E_nu regions; MC scaled to the data POT; shares of all selected MC\n")
        t.write("{0:>9s} {1:>5s} {2:>5s} {3:>11s} {4:>11s} {5:>9s}".format(
            "sample", "E_lo", "E_hi", "data", "MC", "bg"))
        for k, _, _ in BG_CATS:
            t.write(" {0:>9s}".format(k))
        t.write("\n")
        for sample, elo, ehi, dv, mv, bf, shares in rows:
            t.write("{0:>9s} {1:5g} {2:>5s} {3:11.1f} {4:11.2f} {5:9.6f}".format(
                sample, elo, "inf" if ehi is None else "{0:g}".format(ehi), dv, mv, bf))
            for s in shares:
                t.write(" {0:9.6f}".format(s))
            t.write("\n")
    print("  wrote {0}".format(tpath))
    print("")
    print("  Background as a share of the selected MC, by reco E_nu region:")
    hdr = "    {0:>9s} {1:>9s} {2:>7s}".format("sample", "E_nu", "all bg")
    for k, _, _ in BG_CATS:
        hdr += " {0:>8s}".format(k)
    print(hdr)
    for sample, elo, ehi, dv, mv, bf, shares in rows:
        line = "    {0:>9s} {1:>4g}-{2:<4s} {3:7.1%}".format(
            sample, elo, "inf" if ehi is None else "{0:g}".format(ehi), bf)
        for s in shares:
            line += " {0:8.2%}".format(s)
        print(line)


def main():
    a = parse_args()
    if not os.path.isdir(a.outdir):
        os.makedirs(a.outdir)
    print("cv_diagnostics_plots.py v1.0: {0}".format(a.infile))
    f = ROOT.TFile.Open(a.infile)
    if not f or f.IsZombie():
        print("  ERROR: cannot open {0}".format(a.infile))
        sys.exit(2)
    print("  git_describe: {0}".format(stamp(f, "git_describe")))
    label = a.label if a.label is not None else run_label(f)
    mc_pot = f.Get("mc_pot").GetBinContent(1) if f.Get("mc_pot") else 1.0
    data_pot = f.Get("data_pot").GetBinContent(1) if f.Get("data_pot") else 1.0
    mc_scale = data_pot / mc_pot if mc_pot > 0 else 1.0
    print("  POT: mc={0:.4g}  data={1:.4g}  ->  mcScale = {2:.4g}".format(mc_pot, data_pot, mc_scale))
    channels_part(f, a, label, mc_scale)
    background_part(f, a, label, mc_scale)
    print("  done: {0}".format(a.outdir))


if __name__ == "__main__":
    main()
