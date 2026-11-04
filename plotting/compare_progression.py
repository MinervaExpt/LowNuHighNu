#!/usr/bin/env python
"""compare_progression.py -- N-way comparison of processed low-nu results (the ME1D progression).

Generalizes study_thesisComparison.py (thesis vs one processedHists file) to N
labelled inputs, plotted together with ratios to a reference (the first input,
or --ref).

Usage (GPVM, after `source setup_mat.sh`; run from plotting/ or anywhere):
  python compare_progression.py --list <file.root>          # list a file's objects (do this first on the thesis file)
  python compare_progression.py -o <plotdir> LABEL=FILE [LABEL=FILE ...]
        [--playlist ME1D] [--ref LABEL] [--xrange 2 20] [--ratio-range LO HI] [--ratio-errors]
        [--per-pot] [--pot LABEL=POT ...] [--thesis-kincuts MC_MERGED.root DATA_MERGED.root]

Example (the Friday progression; quote labels that contain spaces):
  python compare_progression.py -o /exp/minerva/data/users/finer/MATAna/runs/progression_ME1D \\
    "thesis=/exp/minerva/data/users/finer/highNu/analysisHists/thesisProcessing_noFluxConstraint_yesMuonCVShift_2020-06-13.root" \\
    "Jan 2026=/exp/minerva/data/users/finer/MATAna/2026-02_development/processedHists_2026-02-19.root" \\
    "Oct code=<processedHists of the October-code run>" "v4.4.1=<...wk45...>" "v4.4.1 P9=<...wk45-p9...>"

Inputs.  Each FILE is either a processHists.py output ("current" format,
playlist directories named ME1D) or the thesis-era processed file ("thesis"
format, playlist directories named minervame1D); the format is detected from
the directory names.  For each input the script reads, per sample (inclusive,
lowNu, highNu): the data rate, the efficiency numerator and denominator and the
efficiency, from the --playlist directory; and, from the combined directories:
the low-nu flux, sigma*_low-nu, and the CC-inclusive and high-nu cross sections.
The combined quantities cover every playlist the file was made from: the
legend says "[N playlists]" when that is more than one (the thesis file is
likely all of ME).  Where to find each quantity is the NAMES table below.

The thesis names were confirmed on 2026-11-04 (--list on the thesis file, and a
per-POT check against the thesis mergedHists files).  In the thesis file each
playlist's reco histogram h_ENu_VS_nu_<sample>_reco exists in two cycles under
the same name: ;1 is the data rate and ;2 the MC efficiency numerator (a plain
Get returns ;2), so NAMES gives the cycle explicitly.  A quantity that is
missing from an input is skipped for that input, with a note.

Outputs in <plotdir>:
  <quantity>_<sample>.png   overlay (top) and ratio to the reference (bottom)
  progression.root          every curve and ratio drawn (TH1D), for later reading
  summary.md                per quantity: mean / min / max of the per-bin ratio to
                            the reference over --xrange, and for rates the ratio of
                            the sums.  This is the statistics-free preview of a change.

Notes.
  - Ratios are taken only in bins whose edges match the reference's (the thesis
    E_nu binning has a 22 GeV edge; 0-20 GeV is common to both).
  - Overlays show the CV with its total (stat + syst) error.  Ratio panels show no
    error bars by default: points of a progression share data and MC, so their
    errors are highly correlated, and uncorrelated error bars would overstate the
    scatter.  --ratio-errors draws them anyway (uncorrelated, total errors).
  - Display units: rates per GeV in 10^4 events (as in study_thesisComparison.py),
    or per 10^20 POT with --per-pot; flux per GeV; cross sections in 10^-38 cm^2
    per nucleon; the efficiency as is.  Ratios do not depend on these choices.
  - --per-pot divides the rates by each input's data POT: from <PL>/data_pot in
    processedHists files made by processHists.py from 2026-11-03 on, or --pot LABEL=VALUE.
  - --thesis-kincuts adds the thesis efficiency denominator with the kinematic
    cuts (h_ENu_VS_nu_<sample>_truth_addKinematicCuts in the thesis mergedHists MC
    file), scaled by POT_data / POT_MC from the two files' Meta trees, as
    study_thesisComparison.py does.
"""
import argparse
import math
import os
import re
import sys

SAMPLES = ['inclusive', 'lowNu', 'highNu']

# Where each quantity lives.  {pl} = playlist directory of the format (ME1D or
# minervame1D); {sd} = sample.  The first name found is used.
#   per_sample: one plot per sample, read from the playlist directory
#   kind:       'rate' (event counts; per GeV, POT-normalizable), 'eff', 'flux', 'xsec'
NAMES = [
    # key              title                         per_sample kind    current                                 thesis
    ('dataRate',       'Data event rate',            True,  'rate', ['{pl}/selection_data_enu_ehad_{sd}'],  ['{pl}/h_ENu_VS_nu_{sd}_reco;1']),   # confirmed: cycle 1 = data
    ('effNumerator',   'Efficiency numerator',       True,  'rate', ['{pl}/effnum_enu_ehad_{sd}'],          ['{pl}/h_ENu_VS_nu_{sd}_reco;2']),   # confirmed: cycle 2 = MC
    ('effDenominator', 'Efficiency denominator',     True,  'rate', ['{pl}/effdenom_enu_ehad_{sd}'],        ['{pl}/h_ENu_VS_nu_{sd}_truth']),    # confirmed
    ('eff',            'Efficiency',                 True,  'eff',  ['{pl}/eff_{sd}_{pl}'],                 ['{pl}/eff_{sd}_{pl}']),             # confirmed
    ('fluxLowNu',      'Low-#nu flux',               False, 'flux', ['flux/flux_lowNu_ME'],                 ['flux/flux_lowNu_ME']),             # confirmed
    ('xSectionLowNu',  '#sigma*_{low-#nu}',          False, 'xsec', ['xSections/xSection_lowNu_ME'],        ['xSections/xSection_lowNu_ME']),    # confirmed
    ('xSectionIncl',   'CC-inclusive #sigma',        False, 'xsec', ['xSections/xSection_inclusive_ME'],    ['xSections/xSection_inclusive_ME']),  # confirmed
    ('xSectionHighNu', '#sigma_{high-#nu}',          False, 'xsec', ['xSections/xSection_highNu_ME'],       ['xSections/xSection_highNu_ME']),   # confirmed
]

Y_TITLES = {
    'rate': 'Events / GeV (10^{4})',
    'rate_pot': 'Events / GeV / 10^{20} POT',
    'eff': 'Efficiency',
    'flux': '#Phi / GeV',
    'xsec': '#sigma (10^{-38} cm^{2}/nucleon)',
}
COLORS = [1, 600 + 1, 632 + 1, 416 + 2, 880 + 1, 800 + 7, 432 + 2, 616 + 2]   # kBlack, kBlue+1, kRed+1, kGreen+2, kViolet+1, kOrange+7, kCyan+2, kMagenta+2
MARKERS = [20, 21, 22, 23, 33, 34, 29, 47]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('inputs', nargs='*', metavar='LABEL=FILE')
    ap.add_argument('-o', '--outdir', help='plot directory (created if needed)')
    ap.add_argument('--list', metavar='FILE', help='list every object in FILE (recursively) and exit')
    ap.add_argument('--playlist', default='ME1D', help='playlist for the per-playlist quantities (default ME1D)')
    ap.add_argument('--ref', help='label of the reference input (default: the first)')
    ap.add_argument('--xrange', nargs=2, type=float, default=[2.0, 20.0], metavar=('LO', 'HI'),
                    help='E_nu range drawn and summarized, GeV (default 2 20)')
    ap.add_argument('--ratio-range', nargs=2, type=float, metavar=('LO', 'HI'), help='y range of the ratio panels')
    ap.add_argument('--ratio-errors', action='store_true', help='draw uncorrelated error bars on the ratios')
    ap.add_argument('--per-pot', action='store_true', help='normalize the rates to each input\'s data POT')
    ap.add_argument('--pot', action='append', default=[], metavar='LABEL=POT', help='data POT of an input (with --per-pot)')
    ap.add_argument('--thesis-kincuts', nargs=2, metavar=('MC_MERGED', 'DATA_MERGED'),
                    help='thesis mergedHists MC and data files: adds the kinematic-cut denominator')
    args = ap.parse_args()
    if args.list:
        return args
    if not args.outdir or len(args.inputs) < 2:
        ap.error('need -o PLOTDIR and at least two LABEL=FILE inputs (or --list FILE)')
    specs = []
    for spec in args.inputs:
        label, sep, path = spec.partition('=')
        if not sep or not label or not path:
            ap.error('inputs are LABEL=FILE (got %r)' % spec)
        if label in [l for l, _ in specs]:
            ap.error('duplicate label %r' % label)
        if not os.path.isfile(path):
            ap.error('no such file: %s' % path)
        specs.append((label, path))
    args.specs = specs
    labels = [l for l, _ in specs]
    if args.ref is None:
        args.ref = labels[0]
    elif args.ref not in labels:
        ap.error('--ref %r is not one of the labels %s' % (args.ref, labels))
    pots = {}
    for spec in args.pot:
        label, sep, value = spec.partition('=')
        if not sep or label not in labels:
            ap.error('--pot takes LABEL=VALUE with a known label (got %r)' % spec)
        pots[label] = float(value)
    args.pots = pots
    if not re.match(r'^ME1[A-Z]$', args.playlist):
        ap.error('--playlist must look like ME1D')
    return args


ARGS = parse_args()
sys.argv = sys.argv[:1]   # this script's options are not for ROOT

import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
try:
    from ROOT import PlotUtils  # noqa: F401,E402  (MnvH1D / MnvH2D)
except Exception as exc:  # pragma: no cover
    print('WARNING: could not import PlotUtils (%s)' % exc)


def walk(directory, prefix=''):
    """[(path, class name)] for every object under directory, recursively."""
    out = []
    for k in directory.GetListOfKeys():
        name = prefix + k.GetName()
        cls = ROOT.TClass.GetClass(k.GetClassName())
        if cls and cls.InheritsFrom('TDirectory'):
            out.append((name + '/', k.GetClassName()))
            out.extend(walk(k.ReadObj(), name + '/'))
        else:
            out.append(('%s;%d' % (name, k.GetCycle()), k.GetClassName()))
    return out


def list_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        sys.exit('ERROR: cannot open %s' % path)
    for name, cls in walk(f):
        print('%-12s %s' % (cls, name))
    return 0


def detect_format(f, playlist):
    names = [k.GetName() for k in f.GetListOfKeys()]
    thesis_pl = 'minervame1' + playlist[-1]
    if thesis_pl in names:
        return 'thesis', thesis_pl, sum(1 for n in names if re.match(r'^minervame1[A-Z]$', n))
    if playlist in names:
        return 'current', playlist, sum(1 for n in names if re.match(r'^ME1[A-Z]$', n))
    return None, None, 0


def to_1d(obj, name):
    """CV histogram (TH1D, total error) of a 1D or 2D object; 2D is projected onto x (E_nu) over all y."""
    if obj.InheritsFrom('TH2'):
        obj = obj.ProjectionX(name + '_px')
    if hasattr(obj, 'GetCVHistoWithError'):
        h = ROOT.TH1D(obj.GetCVHistoWithError())
    else:
        h = ROOT.TH1D(obj)
    h.SetName(name)
    h.SetTitle('')
    return h


def edges(h):
    ax = h.GetXaxis()
    return [ax.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)]


def ratio_hist(num, den, name, with_errors):
    """num / den in the bins of num whose edges match a bin of den; other bins left empty (content 0, error 0)."""
    r = num.Clone(name)
    r.Reset()
    den_bins = {}
    de = edges(den)
    for j in range(1, den.GetNbinsX() + 1):
        den_bins[(round(de[j - 1], 9), round(de[j], 9))] = j
    ne = edges(num)
    matched = []
    for i in range(1, num.GetNbinsX() + 1):
        j = den_bins.get((round(ne[i - 1], 9), round(ne[i], 9)))
        if j is None:
            continue
        a, b = num.GetBinContent(i), den.GetBinContent(j)
        if b == 0:
            continue
        rv = a / b
        r.SetBinContent(i, rv)
        if with_errors and a != 0:
            ea, eb = num.GetBinError(i), den.GetBinError(j)
            r.SetBinError(i, abs(rv) * math.sqrt((ea / a) ** 2 + (eb / b) ** 2))
        else:
            r.SetBinError(i, 0.0)
        matched.append(i)
    return r, matched


def meta_pot(path):
    """Sum of POT_Used over the Meta tree of a thesis-era mergedHists file."""
    from array import array
    chain = ROOT.TChain('Meta')
    chain.Add(path)
    used = array('d', [0.0])
    chain.SetBranchAddress('POT_Used', used)
    total = 0.0
    for i in range(chain.GetEntries()):
        chain.GetEntry(i)
        total += used[0]
    return total


def style(h, k):
    h.SetLineColor(COLORS[k % len(COLORS)])
    h.SetMarkerColor(COLORS[k % len(COLORS)])
    h.SetMarkerStyle(MARKERS[k % len(MARKERS)])
    h.SetMarkerSize(1.1)
    h.SetLineWidth(2)


def in_range(h, i, lo, hi):
    return h.GetXaxis().GetBinLowEdge(i) >= lo - 1e-9 and h.GetXaxis().GetBinUpEdge(i) <= hi + 1e-9


def main():
    if ARGS.list:
        return list_file(ARGS.list)

    style_header = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'myPlotStyle.h')
    if os.path.isfile(style_header):
        ROOT.gROOT.ProcessLine('.L %s' % style_header)
        ROOT.myPlotStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    inputs = []   # dicts: label, path, file, fmt, pl, nplaylists, pot
    for label, path in ARGS.specs:
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            sys.exit('ERROR: cannot open %s' % path)
        fmt, pl, npl = detect_format(f, ARGS.playlist)
        if fmt is None:
            sys.exit('ERROR: %s has neither a %s nor a minervame1%s directory: not a processed file for this playlist'
                     % (path, ARGS.playlist, ARGS.playlist[-1]))
        pot = ARGS.pots.get(label)
        if pot is None:
            h = f.Get('%s/data_pot' % pl)
            if h:
                pot = h.GetBinContent(1)
        stamps = []
        for name in ('processHists_git_describe', '%s/git_describe' % pl, '%s/tune' % pl):
            o = f.Get(name)
            if o:
                stamps.append('%s=%s' % (name.split('/')[-1], o.GetTitle()))
        print('%-14s %-7s %d playlist(s)  POT %s  %s\n    %s' % (label, fmt, npl, '%.4g' % pot if pot else 'unknown',
                                                               ' '.join(stamps), path))
        inputs.append(dict(label=label, path=path, file=f, fmt=fmt, pl=pl, npl=npl, pot=pot))

    if ARGS.per_pot:
        missing = [d['label'] for d in inputs if not d['pot']]
        if missing:
            sys.exit('ERROR: --per-pot needs the data POT of %s: give --pot LABEL=VALUE (processedHists made '
                     'before 2026-11-03 carry no POT)' % ', '.join(missing))
    elif all(d['pot'] for d in inputs) and len(set(round(d['pot'] / inputs[0]['pot'], 4) for d in inputs)) > 1:
        print('NOTE: the inputs have different data POT; rates are not comparable without --per-pot')

    kincuts = None
    if ARGS.thesis_kincuts:
        mc_path, data_path = ARGS.thesis_kincuts
        mc_pot, data_pot = meta_pot(mc_path), meta_pot(data_path)
        print('thesis kinematic-cut denominator: POT data %.4g / MC %.4g = %.4g' % (data_pot, mc_pot, data_pot / mc_pot))
        kincuts = dict(file=ROOT.TFile.Open(mc_path), scale=data_pot / mc_pot, pot=data_pot)

    os.makedirs(ARGS.outdir, exist_ok=True)
    out_root = ROOT.TFile(os.path.join(ARGS.outdir, 'progression.root'), 'recreate')
    lo, hi = ARGS.xrange
    summary = ['# Progression summary (%s, %g-%g GeV)' % (ARGS.playlist, lo, hi), '',
               'Reference: **%s**. Per-bin ratios input / reference over bins inside the range whose edges '
               'match the reference; "sum" = ratio of the summed contents (rates only).' % ARGS.ref, '']

    for key, title, per_sample, kind, cur_names, thesis_names in NAMES:
        for sd in (SAMPLES if per_sample else [None]):
            tag = key + ('_' + sd if sd else '')
            curves = []   # (label, TH1D)
            for d in inputs:
                cands = cur_names if d['fmt'] == 'current' else thesis_names
                obj, used = None, None
                for c in cands:
                    n = c.format(pl=d['pl'], sd=sd)
                    obj = d['file'].Get(n)
                    if obj:
                        used = n
                        break
                if not obj:
                    print('  %-26s %-14s not found (tried %s)' % (tag, d['label'],
                          ', '.join(c.format(pl=d['pl'], sd=sd) for c in cands)))
                    continue
                h = to_1d(obj, '%s__%s' % (tag, re.sub(r'\W', '_', d['label'])))
                if kind == 'rate':
                    if ARGS.per_pot:
                        h.Scale(1e20 / d['pot'], 'width')
                    else:
                        h.Scale(1e-4, 'width')
                elif kind == 'flux':
                    h.Scale(1.0, 'width')
                elif kind == 'xsec':
                    h.Scale(1e38)
                label = d['label']
                if not per_sample and d['npl'] > 1:
                    label += ' [%d playlists]' % d['npl']
                curves.append((d['label'], label, h))
            if kincuts is not None and key == 'effDenominator':
                o = kincuts['file'].Get('h_ENu_VS_nu_%s_truth_addKinematicCuts' % sd)
                if o:
                    h = to_1d(o, '%s__thesis_kincuts' % tag)
                    h.Scale(kincuts['scale'])
                    if ARGS.per_pot:
                        h.Scale(1e20 / kincuts['pot'], 'width')
                    else:
                        h.Scale(1e-4, 'width')
                    curves.append(('thesis + kin. cuts', 'thesis + kin. cuts', h))
                else:
                    print('  %-26s thesis + kin. cuts not found' % tag)
            ref = [c for c in curves if c[0] == ARGS.ref]
            if len(curves) < 2 or not ref:
                print('  %-26s skipped (%s)' % (tag, 'no reference' if not ref else 'fewer than two inputs'))
                continue
            ref_h = ref[0][2]
            draw(tag, title, kind, sd, curves, ref_h, out_root, summary)

    out_root.Close()
    with open(os.path.join(ARGS.outdir, 'summary.md'), 'w') as fh:
        fh.write('\n'.join(summary) + '\n')
    print('\n'.join(summary))
    print('Wrote %s/{*.png, progression.root, summary.md}' % ARGS.outdir)
    return 0


def draw(tag, title, kind, sd, curves, ref_h, out_root, summary):
    lo, hi = ARGS.xrange
    canvas = ROOT.TCanvas('c_' + tag, '', 900, 900)
    top = ROOT.TPad('top_' + tag, '', 0.0, 0.32, 1.0, 1.0)
    bot = ROOT.TPad('bot_' + tag, '', 0.0, 0.0, 1.0, 0.32)
    top.SetBottomMargin(0.02)
    top.SetLeftMargin(0.14)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.32)
    bot.SetLeftMargin(0.14)
    top.Draw()
    bot.Draw()

    top.cd()
    ymax = 0.0
    for _, _, h in curves:
        for i in range(1, h.GetNbinsX() + 1):
            if in_range(h, i, lo, hi):
                ymax = max(ymax, h.GetBinContent(i) + h.GetBinError(i))
    legend = ROOT.TLegend(0.50, 0.66, 0.93, 0.89)   # headroom: the y axis goes to 1.6x the highest point
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    legend.SetNColumns(2 if len(curves) > 4 else 1)
    ytitle = Y_TITLES['rate_pot' if (kind == 'rate' and ARGS.per_pot) else kind]
    for k, (_, label, h) in enumerate(curves):
        style(h, k)
        h.GetXaxis().SetRangeUser(lo, hi)
        h.GetXaxis().SetLabelSize(0)
        h.GetYaxis().SetTitle(ytitle)
        h.GetYaxis().SetTitleSize(0.05)
        h.GetYaxis().SetTitleOffset(1.3)
        h.SetMinimum(0.0)
        h.SetMaximum(1.6 * ymax if ymax > 0 else 1.0)
        h.Draw('PE' if k == 0 else 'PE SAME')
        legend.AddEntry(h, label, 'pe')
        out_root.cd()
        h.Write()
    legend.Draw()
    caption = ROOT.TLatex()
    caption.SetNDC()
    caption.SetTextSize(0.045)
    caption.DrawLatex(0.16, 0.93, '%s%s, %s' % (title, ' (%s)' % sd if sd else '', ARGS.playlist))

    bot.cd()
    ratios = []
    rmin, rmax = 1e9, -1e9
    ref_label = ARGS.ref
    summary.append('**%s**' % tag)
    summary.append('')
    summary.append('| input | mean ratio | min | max | sum ratio |')
    summary.append('|---|---|---|---|---|')
    for k, (label, legend_label, h) in enumerate(curves):
        if h is ref_h:
            continue   # identically 1: the dashed line
        r, matched = ratio_hist(h, ref_h, 'ratio_' + h.GetName(), ARGS.ratio_errors)
        style(r, k)
        vals = [r.GetBinContent(i) for i in matched if in_range(r, i, lo, hi)]
        if vals:
            rmin, rmax = min(rmin, min(vals)), max(rmax, max(vals))
            sums = ''
            if kind == 'rate':
                num = sum(h.GetBinContent(i) * h.GetBinWidth(i) for i in matched if in_range(h, i, lo, hi))
                ref_bins = [i for i in range(1, ref_h.GetNbinsX() + 1) if in_range(ref_h, i, lo, hi)]
                den = sum(ref_h.GetBinContent(i) * ref_h.GetBinWidth(i) for i in ref_bins)
                sums = '%.4f' % (num / den) if den else ''
            summary.append('| %s | %.4f | %.4f | %.4f | %s |' % (legend_label, sum(vals) / len(vals), min(vals), max(vals), sums))
        else:
            summary.append('| %s | (no matching bins) | | | |' % legend_label)
        ratios.append(r)
    summary.append('')
    if ARGS.ratio_range:
        rlo, rhi = ARGS.ratio_range
    elif rmax > rmin:
        pad = 0.1 * (rmax - rmin) + 0.02
        rlo, rhi = max(0.0, rmin - pad), rmax + pad
    else:
        rlo, rhi = 0.5, 1.5
    if not ARGS.ratio_range:
        rlo, rhi = min(rlo, 0.95), max(rhi, 1.05)   # keep the line at 1 in view
    for k, r in enumerate(ratios):
        r.GetXaxis().SetRangeUser(lo, hi)
        r.SetMinimum(rlo)
        r.SetMaximum(rhi)
        r.GetXaxis().SetTitle('E_{#nu} (GeV)')
        r.GetXaxis().SetTitleSize(0.11)
        r.GetXaxis().SetLabelSize(0.10)
        r.GetYaxis().SetTitle('ratio to %s' % ref_label)
        r.GetYaxis().SetTitleSize(0.085)
        r.GetYaxis().SetTitleOffset(0.75)
        r.GetYaxis().SetLabelSize(0.09)
        r.GetYaxis().SetNdivisions(505)
        r.GetYaxis().SetNoExponent(True)
        r.Draw(('PE' if ARGS.ratio_errors else 'P') + ('' if k == 0 else ' SAME'))
        out_root.cd()
        r.Write()
    line = ROOT.TLine(lo, 1.0, hi, 1.0)
    line.SetLineStyle(2)
    line.Draw()

    canvas.Print(os.path.join(ARGS.outdir, tag + '.png'))
    canvas.Close()


if __name__ == '__main__':
    sys.exit(main())
