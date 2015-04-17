#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>

using namespace std;

///////////////////////////

void DrawAndSave(const char *namefile ,TH1D* h1, int xmin, int xmax){
   
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetLineColor(kBlue);
  //  h1->SetXTitle("R");
  h1->Draw();
  
  canvas ->SaveAs(namefile);
  canvas ->Destructor();
}

//////////////////////////
void DrawAndSave(const char *namefile ,TH1D* h1, TH1D* h2, int xmin, int xmax){
   
  //  cout<<"working"<<endl;
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->Draw();
  h2->Draw("same");
  
  canvas ->SaveAs(namefile);
  canvas ->Destructor();
}

///////////////////////////////////

void DrawRatioAndSave(const char *namefile, const char *XTitle, TH1D* h1, TH1D* h2, int xmin, int xmax){

 TCanvas *c = new TCanvas("c","",800,800);
    // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
 
  h1->SetStats(0);          // No statistics on upper plot
  //  h1 -> SetYTitle("Normalized to unit");
  h1 -> SetYTitle("Events");
  h1 -> SetTitle(" ");
  h1 -> GetYaxis()->SetTitleOffset(1.55);
  h1 -> SetLineColor(kBlue);
  h1 -> SetMinimum(-0.01);
  //h1 -> SetMaximum(0.12);
  h1 -> GetXaxis()->SetRangeUser(xmin,xmax);
  h2 ->SetLineColor(kRed);
  h1 ->Draw();               // Draw h1
  h2 ->Draw("same");         // Draw h2 on top of h1

  //  TLegend *leg1 = new TLegend( 0.9, 0.9, 0.6, 0.6);
  //leg1 -> SetHeader("Pt:");
  //leg1 -> AddEntry(h1, "WideJet1 MC", "L");
  //leg1 -> AddEntry(h2, "WideJet1 Mio ", "L");
  //leg1->Draw();

  // lower plot will be in pad
  c->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
    
  // Define the ratio plot
  TH1D *h3 = (TH1D*)h2->Clone("h3");
  h3->Divide(h1);

  h3->SetLineColor(kBlack);
  h3->SetMinimum(0.5);  // Define Y ..
  h3->SetMaximum(1.5); // .. range
  //   h3->Sumw2();
  h3->SetStats(0);      // No statistics on lower plot
  h3->SetTitle(""); // Remove the ratio title
  h3->SetMarkerStyle(21);
   
  // X axis ratio plot settings
  h3 -> SetXTitle(XTitle);
  h3->GetXaxis()->SetTitleSize(20);
  h3->GetXaxis()->SetTitleFont(43);
  h3->GetXaxis()->SetTitleOffset(4.);
  h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetXaxis()->SetLabelSize(15);
  h3->GetXaxis()->SetRangeUser(xmin,xmax);
  
  // Y axis ratio plot settings
  h3->GetYaxis()->SetTitle("ratio Smeared/MC ");
  h3->GetYaxis()->SetNdivisions(505);
  h3->GetYaxis()->SetTitleSize(20);
  h3->GetYaxis()->SetTitleFont(43);
  h3->GetYaxis()->SetTitleOffset(1.55);
  h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h3->GetYaxis()->SetLabelSize(15);
  h3->Draw("ep");       // Draw the ratio plot
  
  TF1 *f1 = new TF1("f1","1",xmin,xmax);
   f1->Draw("same");
  
  c -> SaveAs(namefile);
  c->Destructor();
}


////////////////
void Normalizer(TH1D* h1){

  int factor=1;

  h1 -> Scale(factor/h1->Integral());

  cout<<"Integrali normalized : "<<h1->Integral() << endl;
}
////////////////////////


////////////////////////////////////////////

int main(int argc, char* argv[]) {
   
  //Cambia file  
  TFile file1("rootfile_RSGravitonToQuarkQuark_kMpl01_M_3000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM.root");

  TH1D *H_JetSmeared1_Pt     = (TH1D*)file1.Get("H_JetSmeared1_Pt");
  TH1D *H_JetSmeared1_Eta   = (TH1D*)file1.Get("H_JetSmeared1_Eta");
  TH1D *H_JetSmeared1_Phi   = (TH1D*)file1.Get("H_JetSmeared1_Phi");
  TH1D *H_JetSmeared1_M     = (TH1D*)file1.Get("H_JetSmeared1_M");
  H_JetSmeared1_Pt ->Rebin(5); 
  H_JetSmeared1_Eta ->Rebin(5); 
  H_JetSmeared1_Phi ->Rebin(5);  
  H_JetSmeared1_M ->Rebin(5); 
  
  TH1D *H_JetSmeared2_Pt     = (TH1D*)file1.Get("H_JetSmeared2_Pt");
  TH1D *H_JetSmeared2_Eta   = (TH1D*)file1.Get("H_JetSmeared2_Eta");
  TH1D *H_JetSmeared2_Phi   = (TH1D*)file1.Get("H_JetSmeared2_Phi");
  TH1D *H_JetSmeared2_M     = (TH1D*)file1.Get("H_JetSmeared2_M");
  H_JetSmeared2_Pt ->Rebin(5); 
  H_JetSmeared2_Eta ->Rebin(5); 
  H_JetSmeared2_Phi ->Rebin(5); 
  H_JetSmeared2_M ->Rebin(5); 

  TH1D *H_WideJet1_Pt     = (TH1D*)file1.Get("H_WideJet1_Pt");
  TH1D *H_WideJet1_Eta   = (TH1D*)file1.Get("H_WideJet1_Eta");
  TH1D *H_WideJet1_Phi   = (TH1D*)file1.Get("H_WideJet1_Phi");
  TH1D *H_WideJet1_M     = (TH1D*)file1.Get("H_WideJet1_M");
  H_WideJet1_Pt ->Rebin(5); 
  H_WideJet1_Eta ->Rebin(5); 
  H_WideJet1_Phi ->Rebin(5); 
  H_WideJet1_M ->Rebin(5); 

  TH1D *H_WideJet2_Pt     = (TH1D*)file1.Get("H_WideJet2_Pt");
  TH1D *H_WideJet2_Eta   = (TH1D*)file1.Get("H_WideJet2_Eta");
  TH1D *H_WideJet2_Phi   = (TH1D*)file1.Get("H_WideJet2_Phi");
  TH1D *H_WideJet2_M     = (TH1D*)file1.Get("H_WideJet2_M");
  H_WideJet2_Pt ->Rebin(5); 
  H_WideJet2_Eta ->Rebin(5); 
  H_WideJet2_Phi ->Rebin(5); 
  H_WideJet2_M ->Rebin(5); 

  TH1D *H_DijetSmeared_Pt   = (TH1D*)file1.Get("H_DijetSmeared_Pt");
  TH1D *H_DijetSmeared_Eta = (TH1D*)file1.Get("H_DijetSmeared_Eta");
  TH1D *H_DijetSmeared_Phi = (TH1D*)file1.Get("H_DijetSmeared_Phi");
  TH1D *H_DijetSmeared_M   = (TH1D*)file1.Get("H_DijetSmeared_M");
  H_DijetSmeared_Pt->Rebin(5);
  H_DijetSmeared_Eta->Rebin(5);
  H_DijetSmeared_Phi->Rebin(5);
  H_DijetSmeared_M->Rebin(5);

  TH1D *H_dijetWide_Pt   = (TH1D*)file1.Get("H_dijetWide_Pt");
  TH1D *H_dijetWide_Eta = (TH1D*)file1.Get("H_dijetWide_Eta");
  TH1D *H_dijetWide_Phi = (TH1D*)file1.Get("H_dijetWide_Phi");
  TH1D *H_dijetWide_M   = (TH1D*)file1.Get("H_dijetWide_M");
  H_dijetWide_Pt->Rebin(5);
  H_dijetWide_Eta->Rebin(5);
  H_dijetWide_Phi->Rebin(5);
  H_dijetWide_M->Rebin(5);

  cout<<"File exist"<<endl;


  ////////////////////////////////////////////////////////////////////////

    int size_h1 = H_dijetWide_M      -> GetSize();
    int size_h2 = H_DijetSmeared_M-> GetSize();

    cout<<"size_h1: "<<size_h1<<endl;
    cout<<"size_h2: "<<size_h2<<endl;


    TH1D *H_Pull = new TH1D("H_Pull"," ", 1000, -5, 5);
    TH2D *H_Mass_Pull = new TH2D("H_Mass_Pull"," ",1000,0,10000,100, -5, 5);

  //calcolo il pull

  for(int nbin = 0; nbin<size_h1; nbin++){

    double binCenter          = H_dijetWide_M->GetXaxis()->GetBinCenter(nbin);
    double binContent_h1  = H_dijetWide_M -> GetBinContent(nbin);
    double binError_h1       = H_dijetWide_M -> GetBinError(nbin);
    double binContent_h2  = H_DijetSmeared_M -> GetBinContent(nbin);
    double binError_h2       = H_DijetSmeared_M -> GetBinError(nbin);
    
    if(binError_h1 ==0 && binError_h2 == 0) continue;

    cout<<"nbin "<<nbin<<endl;
    cout<<"binCenter " <<binCenter<<endl;
    cout<<"binContent_h1 " <<binContent_h1<<endl;
    cout<<"binError_h1 " <<binError_h1<<endl;
    cout<<"binContent_h2 " <<binContent_h2<<endl;
    cout<<"binError_h2 " <<binError_h2<<endl;
   
    double error_h1_square = binError_h1 * binError_h1 ;
    double error_h2_square = binError_h2 * binError_h2 ;

    double error_diff = sqrt( binError_h1 * binError_h1 + binError_h2* binError_h2 );

    double pull = (binContent_h2 - binContent_h1 ) / error_diff ;
    
    cout<<error_h1_square<<endl;
    cout<<error_h2_square<<endl;
    cout<<error_diff<<endl;
    cout<<"Pull "<<pull<<endl;
    
    H_Pull -> Fill(pull);
    H_Mass_Pull -> Fill(binCenter, pull);
    
  }
  
  TCanvas *cpull = new TCanvas("cpull","",800,800);
  
  //H_Mass_Pull ->GetXaxis() -> SetRangeUser(1500,4000);
  //H_Mass_Pull ->GetYaxis() -> SetRangeUser(-3,3);
  //H_Mass_Pull ->SetMarkerSize(1);
  //H_Mass_Pull ->SetMarkerStyle(21);
  //H_Mass_Pull->Draw("p");

  H_Pull ->SetMarkerSize(1);
  H_Pull ->SetMarkerStyle(21);
  H_Pull->Draw("ep");
  
  cpull -> SaveAs("H_Pull.png");
 
  TCanvas *Canv = new TCanvas("Canv","",800,800);
  // Upper plot will be in pad1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  
  H_dijetWide_M->SetStats(0);          // No statistics on upper plot
  H_dijetWide_M -> SetYTitle("Events");
  H_dijetWide_M -> SetTitle(" ");
  H_dijetWide_M -> GetYaxis()->SetTitleOffset(1.55);
  H_dijetWide_M -> SetLineColor(kBlue);
  H_dijetWide_M -> SetMinimum(-0.01);
  H_dijetWide_M -> GetXaxis()->SetRangeUser(1500,4000);
  H_DijetSmeared_M ->SetLineColor(kRed);
  H_dijetWide_M ->Draw();               // Draw H_dijetWide_M
  H_DijetSmeared_M ->Draw("same");         // Draw H_DijetSmeared_M on top of H_dijetWide_M
  
  //  TLegend *leg1 = new TLegend( 0.9, 0.9, 0.6, 0.6);
  //leg1 -> SetHeader("Pt:");
  //leg1 -> AddEntry(H_dijetWide_M, "WideJet1 MC", "L");
  //leg1 -> AddEntry(H_DijetSmeared_M, "WideJet1 Mio ", "L");
  //leg1->Draw();

  // lower plot will be in pad
  Canv->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  H_Mass_Pull->SetLineColor(kBlack);
  H_Mass_Pull ->GetXaxis() -> SetRangeUser(1500,4000);
  H_Mass_Pull ->GetYaxis() -> SetRangeUser(-4,4);
  H_Mass_Pull ->SetMarkerSize(1);
  H_Mass_Pull ->SetMarkerStyle(21);
  H_Mass_Pull->SetStats(0);      // No statistics on lower plot
  H_Mass_Pull->SetTitle(""); // Remove the ratio title 
  
  // X axis ratio plot settings
 H_Mass_Pull -> SetXTitle("M [GeV]");
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
 
 TF1 *f1 = new TF1("f1","0",1500,4000);
 f1->Draw("same");
 
 Canv -> SaveAs("Histo_WithPull.png");



  ///////////////////////////////////////////////////////////////////////////////////
  // Normalizer(H_JetSmeared1_Pt);
  // Normalizer(H_WideJet1_Pt);
  // Normalizer(H_JetSmeared2_Pt);
  // Normalizer(H_WideJet2_Pt);
  // Normalizer(H_DijetSmeared_M);
  // Normalizer(H_dijetWide_M);

  bool Comparison = true;

  if(Comparison == true){
    DrawRatioAndSave("Compare_Pt_Jet1.png",    "Pt [GeV]",   H_WideJet1_Pt ,     H_JetSmeared1_Pt,    0, 3000) ;
    DrawRatioAndSave("Compare_Eta_Jet1.png",  "#eta",          H_WideJet1_Eta ,   H_JetSmeared1_Eta, -4, 4) ;
    DrawRatioAndSave("Compare_Phi_Jet1.png",  "#phi",         H_WideJet1_Phi ,    H_JetSmeared1_Phi, -4, 4) ;
    DrawRatioAndSave("Compare_M_Jet1.png",    "M [GeV]",   H_WideJet1_M ,      H_JetSmeared1_M,    0, 1000) ;
    
    DrawRatioAndSave("Compare_Pt_Jet2.png",     "Pt [GeV]",   H_WideJet2_Pt ,     H_JetSmeared2_Pt,    0, 3000) ;
    DrawRatioAndSave("Compare_Eta_Jet2.png",   "#eta",          H_WideJet2_Eta ,   H_JetSmeared2_Eta, -4, 4) ;
    DrawRatioAndSave("Compare_Phi_Jet2.png",   "#phi",          H_WideJet2_Phi ,   H_JetSmeared2_Phi, -4, 4) ;
    DrawRatioAndSave("Compare_M_Jet2.png",     "M [GeV]",    H_WideJet2_M ,     H_JetSmeared2_M,    0, 1000) ;
    
    DrawRatioAndSave("Compare_Dijet_Pt.png",     "Pt [GeV]",      H_dijetWide_Pt,     H_DijetSmeared_Pt,    0, 3000) ;
    DrawRatioAndSave("Compare_Dijet_Eta.png",    "#eta",            H_dijetWide_Eta,   H_DijetSmeared_Eta, -4, 4) ;
    DrawRatioAndSave("Compare_Dijet_Phi.png",    "#phi",            H_dijetWide_Phi,   H_DijetSmeared_Phi, -4, 4) ;
    DrawRatioAndSave("Compare_Dijet_M.png",      "M [GeV]",      H_dijetWide_M,     H_DijetSmeared_M,    1500, 4000) ;
  }else{
    cout<<"Non sto facendo confronti"<<endl;
  }

   ///////////////////////////////////////////////////////

  bool SmearingFunctions = false;

  if(SmearingFunctions == true){
  TFile file10("rootFile_QuarkQuark_Bin500.root");
  TH1D *H_step_pt = (TH1D*)file10.Get("H_step_pt");
  int step_pt = H_step_pt->GetMean();
  //  cout<< step_pt<<endl;

  int n_bin = 4500 / step_pt ;   
  char HistoName[200];	 
  TH1D *h[100][100];


  for(int ii=0; ii<5; ii++){
    for(int jj=0; jj<n_bin ; jj++){	 	 

      int pt_bin_max = (jj*step_pt)+step_pt;
      
      sprintf(HistoName,"Histo_RRecowide_Parton_%d_%d",ii,pt_bin_max);

      h[ii][jj]= (TH1D*)file10.Get(HistoName);   
      
      h[ii][jj]->Integral();
      cout<< "ii: "<< ii<< "   jj: "<<jj<<endl;
      cout<<"      h[ii][jj]->Integral(): "<<       h[ii][jj]->Integral() <<endl;

     char extension[10] = ".png";
     char output_Dir[50] = "SmearingFunctions_Bin500/";

      strcat( HistoName, extension);

      strcat(output_Dir,HistoName);

      cout<<output_Dir<<endl;
      DrawAndSave(output_Dir, h[ii][jj], 0., 2.);
      
    }
  }
  }else{
    cout<<"Non sto graficando le funzioni di smearing"<<endl;
  }
  
}
