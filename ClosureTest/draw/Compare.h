#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TLine.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>

using namespace std;

///////////////////////////

//Disegnare e salvare un istogramma
void DrawAndSave(const char *output_dir, const char *nameFile ,TH1F* h1, double xmin, double xmax, const char *XTitle, const char *YTitle, TPaveText *text){
  
  char * nameSaved ;
  nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, output_dir);
  //unisco szStringa3 con szStinga2
  strcat(nameSaved, nameFile);
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1    -> SetTitle(" ");
  //  h1    -> SetLineColor(kColor);
  h1    -> SetStats(0);
  h1    -> SetXTitle(XTitle);
  h1    -> SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.4);
  h1    -> GetXaxis()->SetRangeUser(xmin, xmax);
  h1    -> Draw();
  text->Draw();
  canvas ->SaveAs(nameSaved);
  canvas ->Destructor();
}

//////////////////////////
//Disegnare e salvare il confronto tra due istogrammi
void DrawAndSave(const char *output_dir, const char *nameFile ,TH1D* h1, TH1D* h2, int xmin, int xmax, const char *XTitle, const char *YTitle, TLegend* leg){
   
  char * nameSaved ;
  nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, output_dir);
  //unisco szStringa3 con szStinga2
  strcat(nameSaved, nameFile);
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->SetStats(0);
  h1->SetTitle(" ");
  h1->SetXTitle(XTitle);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1->SetLineColor(kBlue);
  h1->SetMarkerColor(kBlue);
  h1->SetMarkerStyle(2);
  h2->SetLineColor(kRed);
  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(3);
  h1->Draw();
  h2->Draw("same");
  leg->Draw();
  canvas ->SaveAs(nameSaved);
  canvas ->Destructor();
}
///////////////////////////////////////////////
//Disegnare e salvare il confronto tra due istogrammi
void DrawAndSave(const char *output_dir, const char *nameFile ,TH1D* h1, TH1D* h2, int xmin, int xmax, const char *XTitle, const char *YTitle, int h1Color, int h2Color){
   
  char * nameSaved ;
  nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, output_dir);
  //unisco szStringa3 con szStinga2
  strcat(nameSaved, nameFile);
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->SetStats(0);
  h1->SetTitle(" ");
  h1->SetXTitle(XTitle);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1->SetLineColor(h1Color);
  h1->SetMarkerColor(h1Color);
  h1->SetMarkerStyle(2);
  h2->SetLineColor(h2Color);
  h2->SetMarkerColor(h2Color);
  h2->SetMarkerStyle(3);
  h1->Draw();
  h2->Draw("same");
  //leg->Draw();
  canvas ->SaveAs(nameSaved);
  canvas ->Destructor();
  }
///////////////////////////////////
//Disegnare e salvare il confronto tra 4 istogrammi
void DrawAndSave(const char *output_dir, const char *nameFile ,TH1D* h1, TH1D* h2, TH1D *h3, TH1D *h4, int xmin, int xmax, const char *XTitle, const char *YTitle, TLegend* leg){
  char * nameSaved ;
   nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  strcpy(nameSaved, output_dir);
  strcat(nameSaved, nameFile);
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->SetStats(0);
  h1->SetTitle(" ");
  h1->SetXTitle(XTitle);
  h1->SetYTitle(YTitle);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  h3->SetLineColor(kGreen+1);
  h4->SetLineColor(kViolet);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  leg->Draw();
  canvas ->SaveAs(nameSaved);
  canvas ->Destructor();
}

///////////////////////////////////
//Disegnare e salvare il confronto tra due istogrammi con il rapporto sotto
//void DrawRatioAndSave(const char *output_dir, const char *nameFile,  TH1D* h1, TH1D* h2, int xmin, int xmax, const char *XTitle, const char *YTitle, TLegend *leg,  bool ){
void DrawRatioAndSave(const char *output_dir, const char *nameFile,  TH1D* h1, TH1D* h2, int xmin, int xmax, const char *XTitle, const char *YTitle,  bool ){
    char * nameSaved ;
    nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
    //copio stringa 1 in stringa nuova  
    strcpy(nameSaved, output_dir);
    //unisco szStringa3 con szStinga2
    strcat(nameSaved, nameFile);
    
    TCanvas *c = new TCanvas("c","",800,800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    h1->SetStats(0);          // No statistics on upper plot
    h1 -> SetTitle(" ");
    h1 -> SetYTitle(YTitle);
    h1 -> GetYaxis()->SetTitleOffset(1.55);
    h1 -> SetLineColor(kBlue);
    h1 -> SetMinimum(-0.001);
    double ymax;
    if (h1->GetMaximum() > h2->GetMaximum() ) ymax = h1->GetMaximum() ; 
    if (h2->GetMaximum() > h1->GetMaximum() ) ymax = h2->GetMaximum() ; 
    h1 -> SetMaximum(ymax+200);
    h1 -> GetXaxis()->SetRangeUser(xmin,xmax);
    // h1 -> SetMinimum(-0.001);
    h2 ->SetLineColor(kRed);
    h1 ->Draw();               // Draw h1
    h2 ->Draw("same");         // Draw h2 on top of h1

    // leg -> Draw();   
 // lower plot will be in pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    int npoints = 15;
    float x [npoints];
    float ex [npoints];
    float y [npoints];
    float yu95[npoints];
    float yd95[npoints];
    for(Int_t ii = 0; ii < npoints; ii++){
      x[ii] = ii * 1000 - 1010;
      ex[ii] = 0;
      y[ii] = 1;  
      yd95[ii] = 0.05;
      yu95[ii] = 0.05;
    }
    TGraph* Line = new TGraph (npoints, x, y);
    Line->SetLineColor(kRed);
    TGraphAsymmErrors *Band95 = new TGraphAsymmErrors(npoints, x, y, ex, ex, yd95, yu95);
    Band95->GetXaxis()->SetLimits(xmin, xmax);
    Band95->GetYaxis()->SetRangeUser(0.8, 1.2);
    Band95->SetTitle("");
    Band95->GetXaxis()->SetTitle(XTitle);
    Band95->GetXaxis()->SetTitleSize(20);
    Band95->GetXaxis()->SetTitleFont(43);
    Band95->GetXaxis()->SetTitleOffset(4.);
    Band95->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Band95->GetXaxis()->SetLabelSize(15);
    Band95->GetYaxis()->SetTitle("Ratio");
    Band95->GetYaxis()->SetTitleOffset(1.0);
    Band95->GetYaxis()->SetNdivisions(505);
    Band95->GetYaxis()->SetTitleSize(20);
    Band95->GetYaxis()->SetTitleFont(43);
    Band95->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Band95->GetYaxis()->SetLabelSize(15);
    Band95->SetFillColor(8);
    Band95->SetLineColor(10);
    if( true) Band95->Draw("a3");
    Line -> Draw("l");
    // Define the ratio plot
    TH1D *h3 = (TH1D*)h2->Clone("h3");
    h3->Divide(h1);
    h3->SetLineColor(kBlack);
    h3->SetMarkerStyle(21);
    h3->Draw("same ep");       // Draw the ratio plot
    c -> SaveAs(nameSaved);
    c->Destructor();
  }

///////////////////////////////////////////

void DrawPullAndSave(const char *output_dir, const char *nameFile, TH1F* h1, TH1F *h2, const char *XTitle, const char *YTitle){

  bool verbose = false;  

  char * nameSaved ;
  nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, output_dir);
  //unisco szStringa3 con szStinga2
  strcat(nameSaved, nameFile);
  
  //calcolo il pull
  int size_h1 = h1 -> GetSize();
  int size_h2 = h2 -> GetSize();
  
  if(size_h1 == size_h2){
    TH1D *H_Pull = new TH1D("H_Pull"," ", 1000, -5, 5);
    TH2D *H_Mass_Pull = new TH2D("H_Mass_Pull"," ",1000,0,10000,100, -5, 5);
    
    for(int nbin = 0; nbin<size_h1; nbin++){
      double binCenter          = h1->GetXaxis()->GetBinCenter(nbin);
      double binContent_h1  = h1 -> GetBinContent(nbin);
      double binError_h1       = h1 -> GetBinError(nbin);
      double binContent_h2  = h2 -> GetBinContent(nbin);
      double binError_h2       = h2 -> GetBinError(nbin);
      if(binError_h1 ==0 && binError_h2 == 0) continue;
      double error_diff = sqrt( binError_h1 * binError_h1 + binError_h2* binError_h2 );
      double pull = (binContent_h2 - binContent_h1 ) / error_diff ;

      if(verbose){
	cout<<"nbin "<<nbin<<endl;
	cout<<"binCenter " <<binCenter<<endl;
	cout<<"binContent_h1 " <<binContent_h1<<endl;
	cout<<"binError_h1 " <<binError_h1<<endl;
	cout<<"binContent_h2 " <<binContent_h2<<endl;
	cout<<"binError_h2 " <<binError_h2<<endl;
	cout<<"error_diff " <<error_diff<<endl;
	cout<<"pull " <<pull<<endl;
      }

      H_Pull -> Fill(pull);
      H_Mass_Pull -> Fill(binCenter, pull);
    }
    
      //      TCanvas *canvas1 = new TCanvas("canvas1","",800,800);
      //H_Pull ->SetMarkerSize(1);
      //H_Pull ->SetMarkerStyle(21);
      //H_Pull->Draw("ep");
      //canvas -> SaveAs("H_Pull.png");
           
    TCanvas *Canv = new TCanvas("Canv","",800,800);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    h1->SetStats(0);          // No statistics on upper plot
    h1 -> SetTitle(" ");
    h1 -> SetYTitle(YTitle);
    h1 -> GetYaxis()->SetTitleOffset(1.55);
    h1 -> SetLineColor(kBlue);
    h1 -> SetMinimum(-0.01);
    double ymax;
    if (h1->GetMaximum() > h2->GetMaximum() ) ymax = h1->GetMaximum() ; 
    if (h2->GetMaximum() > h1->GetMaximum() ) ymax = h2->GetMaximum() ; 
    h1 -> SetMaximum(ymax+200);
    //    h1 -> GetXaxis()->SetRangeUser(xmin,xmax);
    h2 ->SetLineColor(kRed);
    h1 ->Draw();             
    h2 ->Draw("same");   
    // lower plot will be in pad
    Canv->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    H_Mass_Pull->SetStats(0);      // No statistics on lower plot
    H_Mass_Pull->SetTitle(""); // Remove the ratio title  
    H_Mass_Pull->SetLineColor(kBlack);
    //    H_Mass_Pull ->GetXaxis() -> SetRangeUser(xmin,xmax);
    H_Mass_Pull ->GetYaxis() -> SetRangeUser(-4,4);
    H_Mass_Pull ->SetMarkerSize(1);
    H_Mass_Pull ->SetMarkerStyle(21);
    // X axis ratio plot settings
    H_Mass_Pull -> SetXTitle(XTitle);
    H_Mass_Pull->GetXaxis()->SetTitleSize(20);
    H_Mass_Pull->GetXaxis()->SetTitleFont(43);
    H_Mass_Pull->GetXaxis()->SetTitleOffset(4.);
    H_Mass_Pull->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    H_Mass_Pull->GetXaxis()->SetLabelSize(15);
    
    // Y axis ratio plot settings
    H_Mass_Pull->GetYaxis()->SetTitle("Pull");
    H_Mass_Pull->GetYaxis()->SetNdivisions(512);
    H_Mass_Pull->GetYaxis()->SetTitleSize(20);
    H_Mass_Pull->GetYaxis()->SetTitleFont(43);
    H_Mass_Pull->GetYaxis()->SetTitleOffset(1.55);
    H_Mass_Pull->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    H_Mass_Pull->GetYaxis()->SetLabelSize(15);
    H_Mass_Pull->Draw("p");
    TF1 *f1 = new TF1("f1","0",0,10000);
    f1->Draw("same");
    Canv -> SaveAs(nameSaved);
    delete H_Pull;
    delete H_Mass_Pull;
    delete Canv;
  }
}

///////////////////////////////////////////////////////////

 void Normalizer(TH1D* h1){
   //  char HistoName[200] = GetTitle(h1);
  int factor=1;
  h1 -> Scale(factor/h1->Integral());
  cout<<" Integral normalized to: "<<h1->Integral() << endl;

}

////////////////////////////////////////////////////////////
















