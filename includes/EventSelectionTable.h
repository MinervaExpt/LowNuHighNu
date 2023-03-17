#ifndef EventSelectionPlots_h
#define EventSelectionPlots_h

#include <iostream>
#include <stdexcept>

#include "Cuts.h"
#include "SignalDefinition.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"

//======================================================================
// Declare Functions
//======================================================================
void PlotEPTGraph(std::map<ECuts, double> EffNum, double EffDen,
                  std::map<ECuts, double> PurDen);

void PrintEffPurTable_Cut(const std::string name, const double n_sig,
                          const double n_bg, const double prev_n_sig,
                          const double n_all_sig, const double n_all_bg,
                          const double n_data = -1.,
                          const double mc_scale = -1.);

void PrintEffPurTable(const EventCount signal, const EventCount background,
                      const EventCount data, const double data_pot,
                      const double mc_pot);
//======================================================================

void PlotEPTGraph(std::map<ECuts, double> EffNum, double EffDen,
                  std::map<ECuts, double> PurDen) {
  gStyle->SetOptTitle(1);

  // Prepare eff and pur vectors
  std::vector<double> vCuts, vEff, vPur;
  double i = 0.5;
  for (auto c : kCutsVector) {
    vCuts.push_back(i);
    i += 1;

    double eff = EffNum[c] / EffDen;
    vEff.push_back(eff);

    double pur = EffNum[c] / PurDen[c];
    vPur.push_back(pur);
  }
  std::vector<double> vCutserr(kCutsVector.size(), 0.0);
  std::vector<double> vEfferr(kCutsVector.size(), 0.0);
  std::vector<double> vPurerr(kCutsVector.size(), 0.0);

  // make the tgraphs
  TGraphErrors *gEff = new TGraphErrors(vCuts.size(), &(vCuts[0]), &(vEff[0]),
                                        &(vCutserr[0]), &(vEfferr[0]));
  TGraphErrors *gPur = new TGraphErrors(vCuts.size(), &(vCuts[0]), &(vPur[0]),
                                        &(vCutserr[0]), &(vPurerr[0]));

  TCanvas *c = new TCanvas("c", "comp", 200, 10, 700, 500);
  TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
  p1->SetGrid();

  // label the axes
  TAxis *ax = gEff->GetHistogram()->GetXaxis();
  Double_t x1 = ax->GetBinLowEdge(1);
  Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
  gEff->GetHistogram()->GetXaxis()->Set(kCutsVector.size(), x1, x2);
  int j = 0;
  for (auto c : kCutsVector) {
    gEff->GetHistogram()->GetXaxis()->SetBinLabel(j + 1, GetCutName(c).c_str());
    j += 1;
  }

  gEff->SetTitle("Efficiency/Purity vs Cut");
  gEff->GetYaxis()->SetTitle("\%");
  gEff->GetYaxis()->SetRangeUser(0, 1);

  p1->Draw();
  p1->cd();
  gEff->SetMarkerStyle(21);
  gEff->SetMarkerSize(1);
  gEff->SetMarkerColor(kRed);
  gEff->SetLineColor(kRed);
  gEff->Draw("ALP");
  gPad->Update();

  gPur->SetMarkerStyle(21);
  gPur->SetMarkerSize(1);
  gPur->Draw("LP");
  gPad->Update();

  // Legend
  TLegend *leg = new TLegend(0.50, 0.7, 0.7, 0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(0.036);
  leg->AddEntry(gEff, " Efficiency", "LEP");
  leg->AddEntry(gPur, " Purity", "LEP");
  leg->Draw();
  gPad->Update();

  c->cd();
  c->Print("EffPur_graph.png");
  // c->Print("EffPur_graph.eps");
  // c->Print("EffPur_graph.pdf");
}

// only needs mc_scale, not individual pots
// %6.2f
// 6 == Pad with spaces for minimum length of 6
// .2 == round to 2 decimals places
// f == double (i and d are both int)
void PrintEffPurTable_Cut(const std::string name, const double n_sig,
                          const double n_bg, const double prev_n_sig,
                          const double n_all_sig, const double n_all_bg,
                          const double n_data, const double mc_scale) {
  const double n_all = n_sig + n_bg;
  const double sig_eff = 100. * n_sig / n_all_sig;
  const double prev_sig_eff = 100. * prev_n_sig / n_all_sig;
  const double cut_sig_eff = 100. * sig_eff / prev_sig_eff;
  const double pur = 100. * n_sig / n_all;
  const double bg_eff = 100. * n_bg / n_all_bg;
  const double all_eff = 100. * n_all / (n_all_sig + n_all_bg);
  // double cut_eff = 100.*(1-(prev_n_sig - n_sig)/n_all_sig); // old definition
  if (n_data < 0. || mc_scale < 0.) {
    printf(
        " %s   & %6.2f     & "
        "%6.2f\\% & %6.2f\\% & %6.2f\\% & %6.2f & %6.2f\\% & "
        "%6.2f & %6.2f\\% & %s & %s & %s \\\\ \\hline \n",
        name.c_str(), n_sig, sig_eff, cut_sig_eff, pur, n_bg, bg_eff, n_all,
        all_eff, "NA", "NA", "NA");
  } else {
    const double n_all_pot_norm = n_all * mc_scale;
    const double data_mc_ratio = n_data / n_all_pot_norm;
    printf(
        " %s   & %6.2f     & "
        "%6.2f\\% & %6.2f\\% & %6.2f\\% & %6.2f & %6.2f\\% & "
        "%6.2f     & %6.2f\\% & %6.2f     & %6.2f & %6.2f \\\\ \\hline \n",
        name.c_str(), n_sig, sig_eff, cut_sig_eff, pur, n_bg, bg_eff, n_all,
        all_eff, n_all_pot_norm, n_data, data_mc_ratio);
  }
}

void PrintEffPurTable(const EventCount signal, const EventCount background,
                      const EventCount data, const double data_pot,
                      const double mc_pot) {
  const double mc_scale = data_pot / mc_pot;

  std::cout << "\\input{preamble}" << std::endl;
  std::cout << "\\begin{document}" << std::endl;
  std::cout << "\\begin{landscape}" << std::endl;
  std::cout << "\\begin{sidewaystable}[h]" << std::endl;
  std::cout << "\\footnotesize" << std::endl;
  // std::cout << "\\caption{\\today. DataPOT: " << data_pot << ". MCPOT: " <<
  // mc_pot << ".}" << std::endl;
  printf("\\caption{\\today. DataPOT: %0.2f. MCPOT: %0.2f.}", data_pot, mc_pot);
  std::cout << "\\begin{tabular}{|*{12}{l|}}" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout
      << " & \\multicolumn{4}{c|}{Signal} & \\multicolumn{2}{c|}{Background} & "
         "\\multicolumn{2}{c|}{Total} & \\multicolumn{3}{c|}{Data} \\\\ "
      << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "& N     & Eff     & Cut Eff & Pur    & N         & Eff     & N "
               "        & Eff     & N MC (scale) & N Data    & Data/MC \\\\";
  std::cout << "\\hline" << std::endl;

  double n_all_sig = -99999.;
  double n_all_bg = -99999.;
  try {
    n_all_sig = signal.at(kNoCuts);
    n_all_bg  = background.at(kNoCuts);
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << " " << GetCutName(kNoCuts) << "\n";
    std::cerr << "This usually means that no events passed this cut.\n";
  }
  double prev_n_sig = n_all_sig;
  for (auto i_cut : kCutsVector) {
    if (IsPrecut(i_cut)) {
      try {
        PrintEffPurTable_Cut(GetCutName(i_cut), signal.at(i_cut),
                             background.at(i_cut), prev_n_sig, n_all_sig,
                             n_all_bg);
      } catch (const std::out_of_range &oor) {
        std::cerr << "Out of Range error: " << oor.what() << "  " << i_cut
                  << "\n";
        std::cerr << "This usually means that no events passed this cut.\n";
      }
    } else {
      try {
        PrintEffPurTable_Cut(GetCutName(i_cut), signal.at(i_cut),
                             background.at(i_cut), prev_n_sig, n_all_sig,
                             n_all_bg, data.at(i_cut), mc_scale);
      } catch (const std::out_of_range &oor) {
        std::cerr << "Out of Range error: " << oor.what() << " "
                  << GetCutName(i_cut) << "\n";
        std::cerr << "This usually means that no events passed this cut.\n";
      }
    }
    prev_n_sig = signal.at(i_cut);
  }
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{sidewaystable}" << std::endl;
  std::cout << "\\end{landscape}" << std::endl;
  std::cout << "\\end{document}" << std::endl;
}

#endif  // EventSelectionPlots_h
