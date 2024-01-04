#ifndef MYPLOTSTYLE_H
#define MYPLOTSTYLE_H

#include "TStyle.h"
#include "TColor.h"
#include "lab2rgb.h"

#include "TList.h"
#include "TPad.h"
#include "TH1.h"
#include "TCanvas.h"

//======================================================================
void reMax(TPad* pad, double min=0.1, double mult=1.4)
{
  TH1 *firstHist  = 0;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  Double_t runningMax=-9e99;//Hparam.ymax;
  while (( obj=next() )) {
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      if (!firstHist) firstHist=(TH1*)obj;
      TH1* curHist=(TH1*)obj;
      const double thisMax=curHist->GetBinContent(curHist->GetMaximumBin());
      if (thisMax > runningMax) {
        runningMax=thisMax;
      }
    }
  }
  firstHist->GetYaxis()->SetRangeUser(min, mult*runningMax);
  pad->Update();
}

//======================================================================
void setLabPalette(double lStart=0.1, double lEnd=0.9)
{
  // A colour palette based on the La*b* colour space, and stolen from http://davidad.net/colorviz/ 
  const int nCols=20;
  gStyle->SetNumberContours(999);
  Double_t stops[nCols];
  Double_t reds[nCols];
  Double_t greens[nCols];
  Double_t blues[nCols];

  for(int i=0; i<nCols; ++i){
    float frac=float(i)/nCols;
    stops[i]=frac;

    const double l = lStart+frac*(lEnd-lStart);
    double r=l*(1-0.1)+0.1;
    const double angle = frac*TAU + TAU/2;
    bool wasClipped=true;

    unsigned char red, green, blue;
    while(wasClipped && r>0){
      const double a = sin(angle)*r;
      const double b = cos(angle)*r;
      wasClipped=lab2rgb(&red, &green, &blue,
                         l, a, b, 
                         true);
      r-=0.01;
    }
    reds[i]=red/255.;
    greens[i]=green/255.;
    blues[i]=blue/255.;
  }

  TColor::CreateGradientColorTable(nCols, stops, reds, greens, blues, 999);
}

//======================================================================
void setRedPalette()
{
  // White -> red
  const int NRGBs = 9, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
  Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.99, 0.98, 0.94, 0.80, 0.65, 0.40 };
  Double_t green[NRGBs] = { 0.96, 0.88, 0.73, 0.57, 0.42, 0.23, 0.09, 0.06, 0.00 };
  Double_t blue[NRGBs]  = { 0.94, 0.82, 0.63, 0.45, 0.29, 0.17, 0.11, 0.08, 0.05 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

//======================================================================
void setCorrelationPalette()
{
  // A colour palette that goes blue->white->red, useful for
  // correlation matrices
  const int NRGBs = 3, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.50, 1.00};
  Double_t red[NRGBs]   = { 0.00, 1.00, 1.00};
  Double_t green[NRGBs] = { 0.00, 1.00, 0.00};
  Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}

//==========
void setRainbowToWhitePalette()
{
  // Matt Strait's colour palette that fades out to white at zero
  const int NRGBs = 7, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.05, 0.23, 0.45, 0.60, 0.85, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.00, 1.00, 1.00, 0.33 };
  Double_t green[NRGBs] = { 1.00, 1.00, 0.30, 0.40, 1.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 1.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}


//==========
void setDeepSeaPalette()
{
  const int NRGBs = 9, NCont = 999;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.125, 0.25, 0.375, 0.50, 0.625, 0.75, 0.875, 1.00 };
  Double_t red[NRGBs]   = {  0./255.,  9./255., 13./255., 17./255., 24./255.,  32./255.,  27./255.,  25./255.,  29./255.};
  Double_t green[NRGBs] = {  0./255.,  0./255.,  0./255.,  2./255., 37./255.,  74./255., 113./255., 160./255., 221./255.};
  Double_t blue[NRGBs]  = { 28./255., 42./255., 59./255., 78./255., 98./255., 129./255., 154./255., 184./255., 221./255.};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
}


//==========
void setCoarsePalette()
{
  const int NRGBs = 9, NCont = 9;
  gStyle->SetNumberContours(NCont);
  Double_t stops[NRGBs] = { 0.00, 0.125, 0.25, 0.375, 0.50, 0.625, 0.75, 0.875, 1.00 };
  Double_t red[NRGBs]   = {  27./255.,  217./255., 117./255., 231./255., 102./255.,  230./255.,  166./255.,  102./255., 255./255.};
  Double_t green[NRGBs] = {  158./255.,  95./255.,  112./255.,  41./255., 166./255.,  171./255., 118./255., 102./255., 255./255.};
  Double_t blue[NRGBs]  = { 119./255., 2./255., 179./255., 138./255., 30./255., 2./255., 29./255., 102./255., 255./255.};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
}

//==========
void setCustomRochesterPalette()
{
  const int NRGBs = 20, NCont = 20;
  gStyle->SetNumberContours(NCont);
  // Keep first 12 bins as bottom color, then start gradient
  Double_t stops[NRGBs] = {0.000, 0.053, 0.105, 0.158, 0.211, 0.263, 0.316, 0.368, 0.421, 0.474, 0.526, 0.579, 0.632, 0.684, 0.737, 0.789, 0.842, 0.895, 0.947, 1.000};
  Double_t red[NRGBs]   = {  0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 22./255., 44./255., 68./255., 93./255., 124./255., 160./255., 192./255., 237./255.};
  Double_t green[NRGBs] = {  0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 16./255., 41./255., 67./255., 93./255., 125./255., 162./255., 194./255., 241./255.};
  Double_t blue[NRGBs]  = { 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 100./255., 99./255., 99./255., 93./255.,  68./255.,  44./255.,  26./255.,  74./255.};
  // Keep first 8 bins as bottom color, then start gradient
  //Double_t stops[NRGBs] = {0.000, 0.053, 0.105, 0.158, 0.211, 0.263, 0.316, 0.368, 0.421, 0.474, 0.526, 0.579, 0.632, 0.684, 0.737, 0.789, 0.842, 0.895, 0.947, 1.000};
  //Double_t red[NRGBs]   = {0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255.,0./255., 22./255., 43./255., 65./255., 86./255., 108./255., 129./255., 151./255., 172./255., 194./255., 215./255., 237./255.};
  //Double_t green[NRGBs] = {0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255., 0./255.,0./255., 22./255., 44./255., 66./255., 88./255., 110./255., 131./255., 153./255., 175./255., 197./255., 219./255., 241./255.};
  //Double_t blue[NRGBs]  = {97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255., 97./255.,97./255., 91./255., 84./255., 78./255., 71./255., 65./255., 58./255., 52./255., 45./255., 39./255., 32./255., 26./255.};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  
}

//======================================================================
void myPlotStyle()
{
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.15);

  // use large fonts

  gStyle->SetLabelFont(42);
  gStyle->SetTitleFont(42);

  gStyle->SetTextSize(0.08);

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleX(0.25);
  gStyle->SetTitleFontSize(0.08);
  //gStyle->SetTitleOffset(1.0, "Y");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.1, "Z");
  //gStyle->SetLabelOffset(0.02, "X");
  gStyle->SetLabelOffset(0.01, "X");
  gStyle->SetLabelOffset(0.02, "Y");
  gStyle->SetTitleAlign(23);

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTickLength(0.01, "Y");
  gStyle->SetTickLength(0.02, "X");

  gStyle->SetNdivisions(505, "XYZ");
  gStyle->SetStripDecimals(false);

  TGaxis::SetMaxDigits(2);
}

#endif
