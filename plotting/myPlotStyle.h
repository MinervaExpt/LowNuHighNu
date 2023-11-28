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

//==========
void TColor::SetPalette(Int_t ncolors, Int_t *colors, Float_t alpha)
{
   Int_t i;

   static Int_t paletteType = 0;

   Int_t palette[50] = {19,18,17,16,15,14,13,12,11,20,
                        21,22,23,24,25,26,27,28,29,30, 8,
                        31,32,33,34,35,36,37,38,39,40, 9,
                        41,42,43,44,45,47,48,49,46,50, 2,
                         7, 6, 5, 4, 3, 2,1};

   // set default palette (pad type)
   if (ncolors <= 0) {
      ncolors = 50;
      fgPalette.Set(ncolors);
      for (i=0;i<ncolors;i++) fgPalette.fArray[i] = palette[i];
      paletteType = 1;
      return;
   }

   // set Rainbow Color map. Kept for backward compatibility.
   if (ncolors == 1 && colors == 0) {
      ncolors = 50;
      fgPalette.Set(ncolors);
      for (i=0;i<ncolors-1;i++) fgPalette.fArray[i] = 51+i;
      fgPalette.fArray[ncolors-1] = kRed; // the last color of this palette is red
      paletteType = 2;
      return;
   }

   // High quality palettes (255 levels)
   if (colors == 0 && ncolors>50) {

      if (!fgPalettesList.fN) fgPalettesList.Set(62);        // Right now 62 high quality palettes
      Int_t Idx = (Int_t)fgPalettesList.fArray[ncolors-51];  // High quality palettes indices start at 51

      // This high quality palette has already been created. Reuse it.
      if (Idx > 0) {
         Double_t alphas = 10*(fgPalettesList.fArray[ncolors-51]-Idx);
         Bool_t same_alpha = TMath::Abs(alpha-alphas) < 0.0001;
         if (paletteType == ncolors && same_alpha) return; // The current palette is already this one.
         fgPalette.Set(255); // High quality palettes have 255 entries
         for (i=0;i<255;i++) fgPalette.fArray[i] = Idx+i;
         paletteType = ncolors;

         // restore the palette transparency if needed
          if (alphas>0 && !same_alpha) {
             TColor *ca;
             for (i=0;i<255;i++) {
                ca = gROOT->GetColor(Idx+i);
                ca->SetAlpha(alpha);
             }
             fgPalettesList.fArray[paletteType-51] = (Double_t)Idx+alpha/10.;
          }
         return;
      }

      TColor::InitializeColors();
      Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
      Double_t custom_stops[9] = { 0.0000, 0.6500, 0.7000, 0.7500, 0.8000, 0.8500, 0.9000, 0.9500, 1.0000};

      switch (ncolors) {
      // Deep Sea
      case 51:
         {
            Double_t red[9]   = {  0./255.,  9./255., 13./255., 17./255., 24./255.,  32./255.,  27./255.,  25./255.,  29./255.};
            Double_t green[9] = {  0./255.,  0./255.,  0./255.,  2./255., 37./255.,  74./255., 113./255., 160./255., 221./255.};
            Double_t blue[9]  = { 28./255., 42./255., 59./255., 78./255., 98./255., 129./255., 154./255., 184./255., 221./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Grey Scale
      case 52:
         {
            Double_t red[9]   = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Double_t green[9] = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Double_t blue[9]  = { 0./255., 32./255., 64./255., 96./255., 128./255., 160./255., 192./255., 224./255., 255./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Dark Body Radiator
      case 53:
         {
            Double_t red[9]   = { 0./255., 45./255., 99./255., 156./255., 212./255., 230./255., 237./255., 234./255., 242./255.};
            Double_t green[9] = { 0./255.,  0./255.,  0./255.,  45./255., 101./255., 168./255., 238./255., 238./255., 243./255.};
            Double_t blue[9]  = { 0./255.,  1./255.,  1./255.,   3./255.,   9./255.,   8./255.,  11./255.,  95./255., 230./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Two-color hue (dark blue through neutral gray to bright yellow)
      case 54:
         {
            Double_t red[9]   = {  0./255.,  22./255., 44./255., 68./255., 93./255., 124./255., 160./255., 192./255., 237./255.};
            Double_t green[9] = {  0./255.,  16./255., 41./255., 67./255., 93./255., 125./255., 162./255., 194./255., 241./255.};
            Double_t blue[9]  = { 97./255., 100./255., 99./255., 99./255., 93./255.,  68./255.,  44./255.,  26./255.,  74./255.};
            Idx = TColor::CreateGradientColorTable(9, custom_stops, red, green, blue, 255, alpha);
         }
         break;

      // Rain Bow
      case 55:
         {
            Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
            Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
            Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Inverted Dark Body Radiator
      case 56:
         {
            Double_t red[9]   = { 242./255., 234./255., 237./255., 230./255., 212./255., 156./255., 99./255., 45./255., 0./255.};
            Double_t green[9] = { 243./255., 238./255., 238./255., 168./255., 101./255.,  45./255.,  0./255.,  0./255., 0./255.};
            Double_t blue[9]  = { 230./255.,  95./255.,  11./255.,   8./255.,   9./255.,   3./255.,  1./255.,  1./255., 0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Bird
      case 57:
         {
            Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
            Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
            Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cubehelix
      case 58:
         {
            Double_t red[9]   = { 0.0000, 0.0956, 0.0098, 0.2124, 0.6905, 0.9242, 0.7914, 0.7596, 1.0000};
            Double_t green[9] = { 0.0000, 0.1147, 0.3616, 0.5041, 0.4577, 0.4691, 0.6905, 0.9237, 1.0000};
            Double_t blue[9]  = { 0.0000, 0.2669, 0.3121, 0.1318, 0.2236, 0.6741, 0.9882, 0.9593, 1.0000};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Green Red Violet
      case 59:
         {
            Double_t red[9]   = {13./255., 23./255., 25./255., 63./255., 76./255., 104./255., 137./255., 161./255., 206./255.};
            Double_t green[9] = {95./255., 67./255., 37./255., 21./255.,  0./255.,  12./255.,  35./255.,  52./255.,  79./255.};
            Double_t blue[9]  = { 4./255.,  3./255.,  2./255.,  6./255., 11./255.,  22./255.,  49./255.,  98./255., 208./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Blue Red Yellow
      case 60:
         {
            Double_t red[9]   = {0./255.,  61./255.,  89./255., 122./255., 143./255., 160./255., 185./255., 204./255., 231./255.};
            Double_t green[9] = {0./255.,   0./255.,   0./255.,   0./255.,  14./255.,  37./255.,  72./255., 132./255., 235./255.};
            Double_t blue[9]  = {0./255., 140./255., 224./255., 144./255.,   4./255.,   5./255.,   6./255.,   9./255.,  13./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Ocean
      case 61:
         {
            Double_t red[9]   = { 14./255.,  7./255.,  2./255.,  0./255.,  5./255.,  11./255.,  55./255., 131./255., 229./255.};
            Double_t green[9] = {105./255., 56./255., 26./255.,  1./255., 42./255.,  74./255., 131./255., 171./255., 229./255.};
            Double_t blue[9]  = {  2./255., 21./255., 35./255., 60./255., 92./255., 113./255., 160./255., 185./255., 229./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Color Printable On Grey
      case 62:
         {
            Double_t red[9]   = { 0./255.,   0./255.,   0./255.,  70./255., 148./255., 231./255., 235./255., 237./255., 244./255.};
            Double_t green[9] = { 0./255.,   0./255.,   0./255.,   0./255.,   0./255.,  69./255.,  67./255., 216./255., 244./255.};
            Double_t blue[9]  = { 0./255., 102./255., 228./255., 231./255., 177./255., 124./255., 137./255.,  20./255., 244./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Alpine
      case 63:
         {
            Double_t red[9]   = { 50./255., 56./255., 63./255., 68./255.,  93./255., 121./255., 165./255., 192./255., 241./255.};
            Double_t green[9] = { 66./255., 81./255., 91./255., 96./255., 111./255., 128./255., 155./255., 189./255., 241./255.};
            Double_t blue[9]  = { 97./255., 91./255., 75./255., 65./255.,  77./255., 103./255., 143./255., 167./255., 217./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Aquamarine
      case 64:
         {
            Double_t red[9]   = { 145./255., 166./255., 167./255., 156./255., 131./255., 114./255., 101./255., 112./255., 132./255.};
            Double_t green[9] = { 158./255., 178./255., 179./255., 181./255., 163./255., 154./255., 144./255., 152./255., 159./255.};
            Double_t blue[9]  = { 190./255., 199./255., 201./255., 192./255., 176./255., 169./255., 160./255., 166./255., 190./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Army
      case 65:
         {
            Double_t red[9]   = { 93./255.,   91./255.,  99./255., 108./255., 130./255., 125./255., 132./255., 155./255., 174./255.};
            Double_t green[9] = { 126./255., 124./255., 128./255., 129./255., 131./255., 121./255., 119./255., 153./255., 173./255.};
            Double_t blue[9]  = { 103./255.,  94./255.,  87./255.,  85./255.,  80./255.,  85./255., 107./255., 120./255., 146./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Atlantic
      case 66:
         {
            Double_t red[9]   = { 24./255., 40./255., 69./255.,  90./255., 104./255., 114./255., 120./255., 132./255., 103./255.};
            Double_t green[9] = { 29./255., 52./255., 94./255., 127./255., 150./255., 162./255., 159./255., 151./255., 101./255.};
            Double_t blue[9]  = { 29./255., 52./255., 96./255., 132./255., 162./255., 181./255., 184./255., 186./255., 131./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Aurora
      case 67:
         {
            Double_t red[9]   = { 46./255., 38./255., 61./255., 92./255., 113./255., 121./255., 132./255., 150./255., 191./255.};
            Double_t green[9] = { 46./255., 36./255., 40./255., 69./255., 110./255., 135./255., 131./255.,  92./255.,  34./255.};
            Double_t blue[9]  = { 46./255., 80./255., 74./255., 70./255.,  81./255., 105./255., 165./255., 211./255., 225./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Avocado
      case 68:
         {
            Double_t red[9]   = { 0./255.,  4./255., 12./255.,  30./255.,  52./255., 101./255., 142./255., 190./255., 237./255.};
            Double_t green[9] = { 0./255., 40./255., 86./255., 121./255., 140./255., 172./255., 187./255., 213./255., 240./255.};
            Double_t blue[9]  = { 0./255.,  9./255., 14./255.,  18./255.,  21./255.,  23./255.,  27./255.,  35./255., 101./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Beach
      case 69:
         {
            Double_t red[9]   = { 198./255., 206./255., 206./255., 211./255., 198./255., 181./255., 161./255., 171./255., 244./255.};
            Double_t green[9] = { 103./255., 133./255., 150./255., 172./255., 178./255., 174./255., 163./255., 175./255., 244./255.};
            Double_t blue[9]  = {  49./255.,  54./255.,  55./255.,  66./255.,  91./255., 130./255., 184./255., 224./255., 244./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Black Body
      case 70:
         {
            Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
            Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
            Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Blue Green Yellow
      case 71:
         {
            Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
            Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
            Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Brown Cyan
      case 72:
         {
            Double_t red[9]   = { 68./255., 116./255., 165./255., 182./255., 189./255., 180./255., 145./255., 111./255.,  71./255.};
            Double_t green[9] = { 37./255.,  82./255., 135./255., 178./255., 204./255., 225./255., 221./255., 202./255., 147./255.};
            Double_t blue[9]  = { 16./255.,  55./255., 105./255., 147./255., 196./255., 226./255., 232./255., 224./255., 178./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // CMYK
      case 73:
         {
            Double_t red[9]   = {  61./255.,  99./255., 136./255., 181./255., 213./255., 225./255., 198./255., 136./255., 24./255.};
            Double_t green[9] = { 149./255., 140./255.,  96./255.,  83./255., 132./255., 178./255., 190./255., 135./255., 22./255.};
            Double_t blue[9]  = { 214./255., 203./255., 168./255., 135./255., 110./255., 100./255., 111./255., 113./255., 22./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Candy
      case 74:
         {
            Double_t red[9]   = { 76./255., 120./255., 156./255., 183./255., 197./255., 180./255., 162./255., 154./255., 140./255.};
            Double_t green[9] = { 34./255.,  35./255.,  42./255.,  69./255., 102./255., 137./255., 164./255., 188./255., 197./255.};
            Double_t blue[9]  = { 64./255.,  69./255.,  78./255., 105./255., 142./255., 177./255., 205./255., 217./255., 198./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cherry
      case 75:
         {
            Double_t red[9]   = { 37./255., 102./255., 157./255., 188./255., 196./255., 214./255., 223./255., 235./255., 251./255.};
            Double_t green[9] = { 37./255.,  29./255.,  25./255.,  37./255.,  67./255.,  91./255., 132./255., 185./255., 251./255.};
            Double_t blue[9]  = { 37./255.,  32./255.,  33./255.,  45./255.,  66./255.,  98./255., 137./255., 187./255., 251./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Coffee
      case 76:
         {
            Double_t red[9]   = { 79./255., 100./255., 119./255., 137./255., 153./255., 172./255., 192./255., 205./255., 250./255.};
            Double_t green[9] = { 63./255.,  79./255.,  93./255., 103./255., 115./255., 135./255., 167./255., 196./255., 250./255.};
            Double_t blue[9]  = { 51./255.,  59./255.,  66./255.,  61./255.,  62./255.,  70./255., 110./255., 160./255., 250./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Dark Rain Bow
      case 77:
         {
            Double_t red[9]   = {  43./255.,  44./255., 50./255.,  66./255., 125./255., 172./255., 178./255., 155./255., 157./255.};
            Double_t green[9] = {  63./255.,  63./255., 85./255., 101./255., 138./255., 163./255., 122./255.,  51./255.,  39./255.};
            Double_t blue[9]  = { 121./255., 101./255., 58./255.,  44./255.,  47./255.,  55./255.,  57./255.,  44./255.,  43./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Dark Terrain
      case 78:
         {
            Double_t red[9]   = {  0./255., 41./255., 62./255., 79./255., 90./255., 87./255., 99./255., 140./255., 228./255.};
            Double_t green[9] = {  0./255., 57./255., 81./255., 93./255., 85./255., 70./255., 71./255., 125./255., 228./255.};
            Double_t blue[9]  = { 95./255., 91./255., 91./255., 82./255., 60./255., 43./255., 44./255., 112./255., 228./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Fall
      case 79:
         {
            Double_t red[9]   = { 49./255., 59./255., 72./255., 88./255., 114./255., 141./255., 176./255., 205./255., 222./255.};
            Double_t green[9] = { 78./255., 72./255., 66./255., 57./255.,  59./255.,  75./255., 106./255., 142./255., 173./255.};
            Double_t blue[9]  = { 78./255., 55./255., 46./255., 40./255.,  39./255.,  39./255.,  40./255.,  41./255.,  47./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Fruit Punch
      case 80:
         {
            Double_t red[9]   = { 243./255., 222./255., 201./255., 185./255., 165./255., 158./255., 166./255., 187./255., 219./255.};
            Double_t green[9] = {  94./255., 108./255., 132./255., 135./255., 125./255.,  96./255.,  68./255.,  51./255.,  61./255.};
            Double_t blue[9]  = {   7./255.,  9./255.,   12./255.,  19./255.,  45./255.,  89./255., 118./255., 146./255., 118./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Fuchsia
      case 81:
         {
            Double_t red[9]   = { 19./255., 44./255., 74./255., 105./255., 137./255., 166./255., 194./255., 206./255., 220./255.};
            Double_t green[9] = { 19./255., 28./255., 40./255.,  55./255.,  82./255., 110./255., 159./255., 181./255., 220./255.};
            Double_t blue[9]  = { 19./255., 42./255., 68./255.,  96./255., 129./255., 157./255., 188./255., 203./255., 220./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Grey Yellow
      case 82:
         {
            Double_t red[9]   = { 33./255., 44./255., 70./255.,  99./255., 140./255., 165./255., 199./255., 211./255., 216./255.};
            Double_t green[9] = { 38./255., 50./255., 76./255., 105./255., 140./255., 165./255., 191./255., 189./255., 167./255.};
            Double_t blue[9]  = { 55./255., 67./255., 97./255., 124./255., 140./255., 166./255., 163./255., 129./255.,  52./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Green Brown Terrain
      case 83:
         {
            Double_t red[9]   = { 0./255., 33./255., 73./255., 124./255., 136./255., 152./255., 159./255., 171./255., 223./255.};
            Double_t green[9] = { 0./255., 43./255., 92./255., 124./255., 134./255., 126./255., 121./255., 144./255., 223./255.};
            Double_t blue[9]  = { 0./255., 43./255., 68./255.,  76./255.,  73./255.,  64./255.,  72./255., 114./255., 223./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Green Pink
      case 84:
         {
            Double_t red[9]   = {  5./255.,  18./255.,  45./255., 124./255., 193./255., 223./255., 205./255., 128./255., 49./255.};
            Double_t green[9] = { 48./255., 134./255., 207./255., 230./255., 193./255., 113./255.,  28./255.,   0./255.,  7./255.};
            Double_t blue[9]  = {  6./255.,  15./255.,  41./255., 121./255., 193./255., 226./255., 208./255., 130./255., 49./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Island
      case 85:
         {
            Double_t red[9]   = { 180./255., 106./255., 104./255., 135./255., 164./255., 188./255., 189./255., 165./255., 144./255.};
            Double_t green[9] = {  72./255., 126./255., 154./255., 184./255., 198./255., 207./255., 205./255., 190./255., 179./255.};
            Double_t blue[9]  = {  41./255., 120./255., 158./255., 188./255., 194./255., 181./255., 145./255., 100./255.,  62./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Lake
      case 86:
         {
            Double_t red[9]   = {  57./255.,  72./255.,  94./255., 117./255., 136./255., 154./255., 174./255., 192./255., 215./255.};
            Double_t green[9] = {   0./255.,  33./255.,  68./255., 109./255., 140./255., 171./255., 192./255., 196./255., 209./255.};
            Double_t blue[9]  = { 116./255., 137./255., 173./255., 201./255., 200./255., 201./255., 203./255., 190./255., 187./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Light Temperature
      case 87:
         {
            Double_t red[9]   = {  31./255.,  71./255., 123./255., 160./255., 210./255., 222./255., 214./255., 199./255., 183./255.};
            Double_t green[9] = {  40./255., 117./255., 171./255., 211./255., 231./255., 220./255., 190./255., 132./255.,  65./255.};
            Double_t blue[9]  = { 234./255., 214./255., 228./255., 222./255., 210./255., 160./255., 105./255.,  60./255.,  34./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Light Terrain
      case 88:
         {
            Double_t red[9]   = { 123./255., 108./255., 109./255., 126./255., 154./255., 172./255., 188./255., 196./255., 218./255.};
            Double_t green[9] = { 184./255., 138./255., 130./255., 133./255., 154./255., 175./255., 188./255., 196./255., 218./255.};
            Double_t blue[9]  = { 208./255., 130./255., 109./255.,  99./255., 110./255., 122./255., 150./255., 171./255., 218./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Mint
      case 89:
         {
            Double_t red[9]   = { 105./255., 106./255., 122./255., 143./255., 159./255., 172./255., 176./255., 181./255., 207./255.};
            Double_t green[9] = { 252./255., 197./255., 194./255., 187./255., 174./255., 162./255., 153./255., 136./255., 125./255.};
            Double_t blue[9]  = { 146./255., 133./255., 144./255., 155./255., 163./255., 167./255., 166./255., 162./255., 174./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Neon
      case 90:
         {
            Double_t red[9]   = { 171./255., 141./255., 145./255., 152./255., 154./255., 159./255., 163./255., 158./255., 177./255.};
            Double_t green[9] = { 236./255., 143./255., 100./255.,  63./255.,  53./255.,  55./255.,  44./255.,  31./255.,   6./255.};
            Double_t blue[9]  = {  59./255.,  48./255.,  46./255.,  44./255.,  42./255.,  54./255.,  82./255., 112./255., 179./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Pastel
      case 91:
         {
            Double_t red[9]   = { 180./255., 190./255., 209./255., 223./255., 204./255., 228./255., 205./255., 152./255.,  91./255.};
            Double_t green[9] = {  93./255., 125./255., 147./255., 172./255., 181./255., 224./255., 233./255., 198./255., 158./255.};
            Double_t blue[9]  = { 236./255., 218./255., 160./255., 133./255., 114./255., 132./255., 162./255., 220./255., 218./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Pearl
      case 92:
         {
            Double_t red[9]   = { 225./255., 183./255., 162./255., 135./255., 115./255., 111./255., 119./255., 145./255., 211./255.};
            Double_t green[9] = { 205./255., 177./255., 166./255., 135./255., 124./255., 117./255., 117./255., 132./255., 172./255.};
            Double_t blue[9]  = { 186./255., 165./255., 155./255., 135./255., 126./255., 130./255., 150./255., 178./255., 226./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Pigeon
      case 93:
         {
            Double_t red[9]   = { 39./255., 43./255., 59./255., 63./255., 80./255., 116./255., 153./255., 177./255., 223./255.};
            Double_t green[9] = { 39./255., 43./255., 59./255., 74./255., 91./255., 114./255., 139./255., 165./255., 223./255.};
            Double_t blue[9]  = { 39./255., 50./255., 59./255., 70./255., 85./255., 115./255., 151./255., 176./255., 223./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Plum
      case 94:
         {
            Double_t red[9]   = { 0./255., 38./255., 60./255., 76./255., 84./255., 89./255., 101./255., 128./255., 204./255.};
            Double_t green[9] = { 0./255., 10./255., 15./255., 23./255., 35./255., 57./255.,  83./255., 123./255., 199./255.};
            Double_t blue[9]  = { 0./255., 11./255., 22./255., 40./255., 63./255., 86./255.,  97./255.,  94./255.,  85./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Red Blue
      case 95:
         {
            Double_t red[9]   = { 94./255., 112./255., 141./255., 165./255., 167./255., 140./255.,  91./255.,  49./255.,  27./255.};
            Double_t green[9] = { 27./255.,  46./255.,  88./255., 135./255., 166./255., 161./255., 135./255.,  97./255.,  58./255.};
            Double_t blue[9]  = { 42./255.,  52./255.,  81./255., 106./255., 139./255., 158./255., 155./255., 137./255., 116./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Rose
      case 96:
         {
            Double_t red[9]   = { 30./255., 49./255., 79./255., 117./255., 135./255., 151./255., 146./255., 138./255., 147./255.};
            Double_t green[9] = { 63./255., 60./255., 72./255.,  90./255.,  94./255.,  94./255.,  68./255.,  46./255.,  16./255.};
            Double_t blue[9]  = { 18./255., 28./255., 41./255.,  56./255.,  62./255.,  63./255.,  50./255.,  36./255.,  21./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Rust
      case 97:
         {
            Double_t red[9]   = {  0./255., 30./255., 63./255., 101./255., 143./255., 152./255., 169./255., 187./255., 230./255.};
            Double_t green[9] = {  0./255., 14./255., 28./255.,  42./255.,  58./255.,  61./255.,  67./255.,  74./255.,  91./255.};
            Double_t blue[9]  = { 39./255., 26./255., 21./255.,  18./255.,  15./255.,  14./255.,  14./255.,  13./255.,  13./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Sandy Terrain
      case 98:
         {
            Double_t red[9]   = { 149./255., 140./255., 164./255., 179./255., 182./255., 181./255., 131./255., 87./255., 61./255.};
            Double_t green[9] = {  62./255.,  70./255., 107./255., 136./255., 144./255., 138./255., 117./255., 87./255., 74./255.};
            Double_t blue[9]  = {  40./255.,  38./255.,  45./255.,  49./255.,  49./255.,  49./255.,  38./255., 32./255., 34./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Sienna
      case 99:
         {
            Double_t red[9]   = { 99./255., 112./255., 148./255., 165./255., 179./255., 182./255., 183./255., 183./255., 208./255.};
            Double_t green[9] = { 39./255.,  40./255.,  57./255.,  79./255., 104./255., 127./255., 148./255., 161./255., 198./255.};
            Double_t blue[9]  = { 15./255.,  16./255.,  18./255.,  33./255.,  51./255.,  79./255., 103./255., 129./255., 177./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Solar
      case 100:
         {
            Double_t red[9]   = { 99./255., 116./255., 154./255., 174./255., 200./255., 196./255., 201./255., 201./255., 230./255.};
            Double_t green[9] = {  0./255.,   0./255.,   8./255.,  32./255.,  58./255.,  83./255., 119./255., 136./255., 173./255.};
            Double_t blue[9]  = {  5./255.,   6./255.,   7./255.,   9./255.,   9./255.,  14./255.,  17./255.,  19./255.,  24./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // South West
      case 101:
         {
            Double_t red[9]   = { 82./255., 106./255., 126./255., 141./255., 155./255., 163./255., 142./255., 107./255.,  66./255.};
            Double_t green[9] = { 62./255.,  44./255.,  69./255., 107./255., 135./255., 152./255., 149./255., 132./255., 119./255.};
            Double_t blue[9]  = { 39./255.,  25./255.,  31./255.,  60./255.,  73./255.,  68./255.,  49./255.,  72./255., 188./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Starry Night
      case 102:
         {
            Double_t red[9]   = { 18./255., 29./255., 44./255.,  72./255., 116./255., 158./255., 184./255., 208./255., 221./255.};
            Double_t green[9] = { 27./255., 46./255., 71./255., 105./255., 146./255., 177./255., 189./255., 190./255., 183./255.};
            Double_t blue[9]  = { 39./255., 55./255., 80./255., 108./255., 130./255., 133./255., 124./255., 100./255.,  76./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Sunset
      case 103:
         {
            Double_t red[9]   = { 0./255., 48./255., 119./255., 173./255., 212./255., 224./255., 228./255., 228./255., 245./255.};
            Double_t green[9] = { 0./255., 13./255.,  30./255.,  47./255.,  79./255., 127./255., 167./255., 205./255., 245./255.};
            Double_t blue[9]  = { 0./255., 68./255.,  75./255.,  43./255.,  16./255.,  22./255.,  55./255., 128./255., 245./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Temperature Map
      case 104:
         {
            Double_t red[9]   = {  34./255.,  70./255., 129./255., 187./255., 225./255., 226./255., 216./255., 193./255., 179./255.};
            Double_t green[9] = {  48./255.,  91./255., 147./255., 194./255., 226./255., 229./255., 196./255., 110./255.,  12./255.};
            Double_t blue[9]  = { 234./255., 212./255., 216./255., 224./255., 206./255., 110./255.,  53./255.,  40./255.,  29./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Thermometer
      case 105:
         {
            Double_t red[9]   = {  30./255.,  55./255., 103./255., 147./255., 174./255., 203./255., 188./255., 151./255., 105./255.};
            Double_t green[9] = {   0./255.,  65./255., 138./255., 182./255., 187./255., 175./255., 121./255.,  53./255.,   9./255.};
            Double_t blue[9]  = { 191./255., 202./255., 212./255., 208./255., 171./255., 140./255.,  97./255.,  57./255.,  30./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Valentine
      case 106:
         {
            Double_t red[9]   = { 112./255., 97./255., 113./255., 125./255., 138./255., 159./255., 178./255., 188./255., 225./255.};
            Double_t green[9] = {  16./255., 17./255.,  24./255.,  37./255.,  56./255.,  81./255., 110./255., 136./255., 189./255.};
            Double_t blue[9]  = {  38./255., 35./255.,  46./255.,  59./255.,  78./255., 103./255., 130./255., 152./255., 201./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Visible Spectrum
      case 107:
         {
            Double_t red[9]   = { 18./255.,  72./255.,   5./255.,  23./255.,  29./255., 201./255., 200./255., 98./255., 29./255.};
            Double_t green[9] = {  0./255.,   0./255.,  43./255., 167./255., 211./255., 117./255.,   0./255.,  0./255.,  0./255.};
            Double_t blue[9]  = { 51./255., 203./255., 177./255.,  26./255.,  10./255.,   9./255.,   8./255.,  3./255.,  0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Water Melon
      case 108:
         {
            Double_t red[9]   = { 19./255., 42./255., 64./255.,  88./255., 118./255., 147./255., 175./255., 187./255., 205./255.};
            Double_t green[9] = { 19./255., 55./255., 89./255., 125./255., 154./255., 169./255., 161./255., 129./255.,  70./255.};
            Double_t blue[9]  = { 19./255., 32./255., 47./255.,  70./255., 100./255., 128./255., 145./255., 130./255.,  75./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Cool
      case 109:
         {
            Double_t red[9]   = {  33./255.,  31./255.,  42./255.,  68./255.,  86./255., 111./255., 141./255., 172./255., 227./255.};
            Double_t green[9] = { 255./255., 175./255., 145./255., 106./255.,  88./255.,  55./255.,  15./255.,   0./255.,   0./255.};
            Double_t blue[9]  = { 255./255., 205./255., 202./255., 203./255., 208./255., 205./255., 203./255., 206./255., 231./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Copper
      case 110:
         {
            Double_t red[9]   = { 0./255., 25./255., 50./255., 79./255., 110./255., 145./255., 181./255., 201./255., 254./255.};
            Double_t green[9] = { 0./255., 16./255., 30./255., 46./255.,  63./255.,  82./255., 101./255., 124./255., 179./255.};
            Double_t blue[9]  = { 0./255., 12./255., 21./255., 29./255.,  39./255.,  49./255.,  61./255.,  74./255., 103./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Gist Earth
      case 111:
         {
            Double_t red[9]   = { 0./255., 13./255.,  30./255.,  44./255.,  72./255., 120./255., 156./255., 200./255., 247./255.};
            Double_t green[9] = { 0./255., 36./255.,  84./255., 117./255., 141./255., 153./255., 151./255., 158./255., 247./255.};
            Double_t blue[9]  = { 0./255., 94./255., 100./255.,  82./255.,  56./255.,  66./255.,  76./255., 131./255., 247./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      // Viridis
      case 112:
         {
            Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
            Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
            Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};
            Idx = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
         }
         break;

      default:
         ::Error("SetPalette", "Unknown palette number %d", ncolors);
         return;
      }
      paletteType = ncolors;
      if (Idx>0) fgPalettesList.fArray[paletteType-51] = (Double_t)Idx;
      else       fgPalettesList.fArray[paletteType-51] = 0.;
      if (alpha > 0.) fgPalettesList.fArray[paletteType-51] += alpha/10.;
      return;
   }

   // set user defined palette
   if (colors)  {
      fgPalette.Set(ncolors);
      for (i=0;i<ncolors;i++) fgPalette.fArray[i] = colors[i];
   } else {
      fgPalette.Set(TMath::Min(50,ncolors));
      for (i=0;i<TMath::Min(50,ncolors);i++) fgPalette.fArray[i] = palette[i];
   }
   paletteType = 3;
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

  //gStyle->SetPalette(54); //kBlueYellow
  gStyle->SetPalette(57); //kBird

  TGaxis::SetMaxDigits(2);
}

#endif
