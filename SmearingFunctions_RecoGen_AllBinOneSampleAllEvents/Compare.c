#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TStyle.h>

using namespace std;

///////////////////////////

void DrawAndSave(const char *namefile ,TH1D* h1, int xmin, int xmax, int kColor, const char *XTitle, const char *YTitle){
   
  //  cout<<"working"<<endl;
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1    -> SetStats(0);
  h1    -> SetTitle(" ");
  h1     ->GetXaxis()->SetRangeUser(xmin,xmax);
  h1    -> SetLineColor(kColor);
  h1    -> SetXTitle(XTitle);
  h1    -> SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1    -> GetXaxis()->SetRangeUser(xmin, xmax);
  h1    -> Draw();
  
  canvas ->SaveAs(namefile);
  canvas ->Destructor();
}


//////////////////////////
void DrawAndSave(const char *namefile ,TH1D* h1, TH1D* h2, int xmin, int xmax, const char *XTitle, const char *YTitle, TLegend* leg){
   
  //  cout<<"working"<<endl;
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->SetStats(0);
  h1->SetTitle(" ");
  h1->SetXTitle(XTitle);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->Draw();
  h2->Draw("same");
  
  leg->Draw();

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
  h1 -> SetMinimum(-0.001);
  //h1 -> SetMaximum(0.12);
  h1 -> GetXaxis()->SetRangeUser(xmin,xmax);
  h1 -> SetMinimum(-0.001);
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
  h3->SetMinimum(0);  // Define Y ..
  h3->SetMaximum(5); // .. range
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
   
  
  TFile file1("rootfile_RSGravitonToQuarkQuark_kMpl01_M_3000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM.root");

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

  TH1D *H_dijetWide_Pt   = (TH1D*)file1.Get("H_dijetWide_Pt");
  TH1D *H_dijetWide_Eta = (TH1D*)file1.Get("H_dijetWide_Eta");
  TH1D *H_dijetWide_Phi = (TH1D*)file1.Get("H_dijetWide_Phi");
  TH1D *H_dijetWide_M   = (TH1D*)file1.Get("H_dijetWide_M");
  H_dijetWide_Pt   ->Rebin(5);
  H_dijetWide_Eta ->Rebin(5);
  H_dijetWide_Phi ->Rebin(5);
  H_dijetWide_M   ->Rebin(5);

  TH1D *H_GenWideJet1_Pt       = (TH1D*)file1.Get("H_GenWideJet1_Pt");
  TH1D *H_GenWideJet1_Eta     = (TH1D*)file1.Get("H_GenWideJet1_Eta");
  H_GenWideJet1_Pt   -> Rebin(5);
  H_GenWideJet1_Eta -> Rebin(5);

  TH1D *H_GenWideJet2_Pt       = (TH1D*)file1.Get("H_GenWideJet2_Pt");
  TH1D *H_GenWideJet2_Eta     = (TH1D*)file1.Get("H_GenWideJet2_Eta");
  H_GenWideJet2_Pt   -> Rebin(5);
  H_GenWideJet2_Eta -> Rebin(5);

  TH1D *H_invariant_mass_AllGenerator                          = (TH1D*)file1.Get("H_invariant_mass_AllGenerator");
  TH1D *H_invariant_mass_PartonEnergy_RecoDirection  = (TH1D*)file1.Get("H_invariant_mass_PartonEnergy_RecoDirection");
  TH1D *H_invariant_mass_PartonDirection_RecoEnergy  = (TH1D*)file1.Get("H_invariant_mass_PartonDirection_RecoEnergy");
  H_invariant_mass_AllGenerator                           -> Rebin(5);
  H_invariant_mass_PartonEnergy_RecoDirection  -> Rebin(5);
  H_invariant_mass_PartonDirection_RecoEnergy  -> Rebin(5);

  TH1D *H_invariant_mass_GenWideJets                          = (TH1D*)file1.Get("H_invariant_mass_GenWideJets");
  TH1D *H_invariant_mass_GenEnergy_RecoDirection      = (TH1D*)file1.Get("H_invariant_mass_GenEnergy_RecoDirection");
  TH1D *H_invariant_mass_GenDirection_RecoEnergy      = (TH1D*)file1.Get("H_invariant_mass_GenDirection_RecoEnergy");
  TH1D *H_invariant_mass_AllReco                                   = (TH1D*)file1.Get("H_invariant_mass_AllReco");
  H_invariant_mass_GenWideJets                          -> Rebin(5);
  H_invariant_mass_GenEnergy_RecoDirection      -> Rebin(5);
  H_invariant_mass_GenDirection_RecoEnergy      -> Rebin(5);
  H_invariant_mass_AllReco                                  -> Rebin(5);

  TH1D *H_invariant_mass_partons_Complete          = (TH1D*)file1.Get("H_invariant_mass_partons_Complete");
  TH1D *H_invariant_mass_partons_Approx              = (TH1D*)file1.Get("H_invariant_mass_partons_Approx");
  TH1D *H_invariant_mass_GenWideJets_Complete  =  (TH1D*)file1.Get("H_invariant_mass_GenWideJets_Complete");
  TH1D *H_invariant_mass_GenWideJets_Approx      = (TH1D*)file1.Get("H_invariant_mass_GenWideJets_Approx");
  TH1D *H_invariant_mass_WideJets_Complete        = (TH1D*)file1.Get("H_invariant_mass_WideJets_Complete");
  TH1D *H_invariant_mass_WideJets_Approx            = (TH1D*)file1.Get("H_invariant_mass_WideJets_Approx");
  H_invariant_mass_partons_Complete         -> Rebin(5);
  H_invariant_mass_partons_Approx             -> Rebin(5);
  H_invariant_mass_GenWideJets_Complete -> Rebin(5);
  H_invariant_mass_GenWideJets_Approx     -> Rebin(5);
  H_invariant_mass_WideJets_Complete        -> Rebin(5);
  H_invariant_mass_WideJets_Approx            -> Rebin(5);

  TH1D *H_diparton_M          = (TH1D*)file1.Get("H_diparton_M");
  TH1D *H_Genwdijet_M_R11          = (TH1D*)file1.Get("H_GendijetWide_M_R11");
  TH1D *H_Genwdijet_M_R16          = (TH1D*)file1.Get("H_GendijetWide_M_R16");
  TH1D *H_Genwdijet_M_R21          = (TH1D*)file1.Get("H_GendijetWide_M_R21");
  TH1D *H_Genwdijet_M_R26          = (TH1D*)file1.Get("H_GendijetWide_M_R26");
  TH1D *H_Genwdijet_M_R31          = (TH1D*)file1.Get("H_GendijetWide_M_R31");
  H_diparton_M -> Rebin(5);
  H_Genwdijet_M_R11 -> Rebin(5);
  H_Genwdijet_M_R16 -> Rebin(5);
  H_Genwdijet_M_R21 -> Rebin(5);
  H_Genwdijet_M_R26 -> Rebin(5);
  H_Genwdijet_M_R31 -> Rebin(5);

  TH1D *H_DeltaR_WideJet1_parton1    = (TH1D*)file1.Get("H_DeltaR_WideJet1_parton1");
  TH1D *H_DeltaR_WideJet2_parton2    = (TH1D*)file1.Get("H_DeltaR_WideJet2_parton2");
  TH1D *H_DeltaR_WideJet1_GenWideJet1    = (TH1D*)file1.Get("H_DeltaR_WideJet1_GenWideJet1");
  TH1D *H_DeltaR_WideJet2_GenWideJet2    = (TH1D*)file1.Get("H_DeltaR_WideJet2_GenWideJet2");

  cout<<"File exist"<<endl;

  //cosa deve fare il programma
  bool test_RCono = true;
  bool test_formula = false;
  bool test_parton_VS_Reco = false;
  bool test_Gen_VS_Reco = false;
  bool SmearingFunctions = false;

  ////////////////////////////////////////////////////////////////////
  // non normalizzo perche gi eventi si parlano 1:1
  // Normalizer(H_WideJet1_Pt);
  // Normalizer(H_WideJet2_Pt);
  // Normalizer(H_dijetWide_M);

  ////////////////////////////////////////////////////
  
  if( test_RCono == true){
  TGaxis::SetMaxDigits(4);
  TCanvas *Canvas = new TCanvas("Canvas","",800,800);
  
  H_diparton_M    -> SetStats(0);
  H_diparton_M    -> SetXTitle("M [GeV]");
  H_diparton_M    -> SetYTitle("Events");
  H_diparton_M    -> GetYaxis()->SetTitleOffset(1.40);
  H_diparton_M    -> SetTitle(" ");
  H_diparton_M    -> GetXaxis()->SetRangeUser(1500,4000);
  
  
  H_diparton_M           -> SetMarkerColor(kBlack);
  H_Genwdijet_M_R11           -> SetMarkerColor(kBlue);
  H_Genwdijet_M_R16           -> SetMarkerColor(kRed);
  H_Genwdijet_M_R21           -> SetMarkerColor(kGreen+1);
  H_Genwdijet_M_R26           -> SetMarkerColor(kViolet);
  H_Genwdijet_M_R31           -> SetMarkerColor(kYellow+1);
  
  H_diparton_M           -> SetMarkerStyle(2);
  H_Genwdijet_M_R11           -> SetMarkerStyle(3);
  H_Genwdijet_M_R16           -> SetMarkerStyle(4);
  H_Genwdijet_M_R21           -> SetMarkerStyle(5);
  H_Genwdijet_M_R26           -> SetMarkerStyle(6);
  H_Genwdijet_M_R31           -> SetMarkerStyle(7);
  
  H_diparton_M           -> SetLineColor(kBlack);
  H_Genwdijet_M_R11           -> SetLineColor(kBlue);
  H_Genwdijet_M_R16           -> SetLineColor(kRed);
  H_Genwdijet_M_R21           -> SetLineColor(kGreen+1);
  H_Genwdijet_M_R26           -> SetLineColor(kViolet);
  H_Genwdijet_M_R31           -> SetLineColor(kYellow+1);
  
  H_diparton_M          ->Draw();
  H_Genwdijet_M_R11          ->Draw("same");
  H_Genwdijet_M_R16          ->Draw("same");
  H_Genwdijet_M_R21          ->Draw("same");
  H_Genwdijet_M_R26          ->Draw("same");
  H_Genwdijet_M_R31          ->Draw("same");
  
  TLegend *Leg1 = new TLegend(0.1, 0.9, 0.5, 0.7);
  Leg1 -> SetHeader("Massa invariante:");
  Leg1  ->AddEntry(H_diparton_M,"Partoni","PL");
  Leg1  ->AddEntry(H_Genwdijet_M_R11,"Gen widejets con R=1.1","PL");
  Leg1  ->AddEntry(H_Genwdijet_M_R16,"Gen widejets con R=1.6","PL");
  Leg1  ->AddEntry(H_Genwdijet_M_R21,"Gen widejets con R=2.1","PL");
  Leg1  ->AddEntry(H_Genwdijet_M_R26,"Gen widejets con R=2.6","PL");
  Leg1  ->AddEntry(H_Genwdijet_M_R31,"Gen widejets con R=3.1","PL");
  
  Leg1->Draw();  
  
  Canvas ->SaveAs("H_GenWideJets_M_RCono.png");
  
  TLegend *Leg2 = new TLegend(0.9, 0.9, 0.6, 0.7);
  Leg2 -> SetHeader("dR:");
  Leg2  ->AddEntry(H_DeltaR_WideJet1_parton1,"WideJet - Partone","L");
  Leg2  ->AddEntry(H_DeltaR_WideJet1_GenWideJet1,"WideJet - GenWideJet","L");
  
  DrawAndSave("DeltaR_1.root", H_DeltaR_WideJet1_GenWideJet1, H_DeltaR_WideJet1_parton1, 0, 1, "dR", "Events",Leg2);
  DrawAndSave("DeltaR_2.root", H_DeltaR_WideJet2_GenWideJet2, H_DeltaR_WideJet2_parton2, 0, 1, "dR", "Events",Leg2);

  DrawAndSave("H_GenWideJet1_Pt.png", H_GenWideJet1_Pt, 0, 2500, kBlue, "Pt [GeV]", "Events");
  DrawAndSave("H_GenWideJet2_Pt.png", H_GenWideJet2_Pt, 0, 2500, kBlue, "Pt [GeV]", "Events");
  DrawAndSave("H_GenWideJet1_Eta.png", H_GenWideJet1_Eta, -4, 4, kBlue, "Pt [GeV]", "Events");
  DrawAndSave("H_GenWideJet2_Eta.png", H_GenWideJet2_Eta, -4, 4, kBlue, "Pt [GeV]", "Events");

  }















  ///////////////////////////////////////////////////////////////////////
  if(test_formula == true){

  TLegend *Leg = new TLegend(0.1, 0.9, 0.5, 0.7);
  Leg  -> SetHeader("Massa invariante calcolata con:");
  Leg  ->AddEntry(H_invariant_mass_partons_Complete,"Formula completa (con masse)","L");
  Leg  ->AddEntry(H_invariant_mass_partons_Approx,"Formula approssimata (senza masse)","L");

  //confronti
  DrawAndSave("H_invariant_mass_partons.png" , H_invariant_mass_partons_Complete, H_invariant_mass_partons_Approx, 1500, 4000, "M [GeV]", "Events",Leg);
  DrawAndSave("H_invariant_mass_GenWideJets.png" , H_invariant_mass_GenWideJets_Complete, H_invariant_mass_GenWideJets_Approx, 1500, 4000, "M [GeV]", "Events",Leg);
  DrawAndSave("H_invariant_mass_WideJets.png" , H_invariant_mass_WideJets_Complete, H_invariant_mass_WideJets_Approx, 1500, 4000, "M [GeV]", "Events",Leg);

  //formula completa
  TCanvas *Canvas = new TCanvas("Canvas","",800,800);
  H_invariant_mass_partons_Complete    -> SetStats(0);
  H_invariant_mass_partons_Complete    -> SetXTitle("M [GeV]");
  H_invariant_mass_partons_Complete    -> SetYTitle("Events");
  H_invariant_mass_partons_Complete    -> GetYaxis()->SetTitleOffset(1.55);
  H_invariant_mass_partons_Complete    -> SetTitle(" ");
  H_invariant_mass_partons_Complete    -> GetXaxis()->SetRangeUser(1500,4000);
  H_invariant_mass_partons_Complete         -> SetLineColor(kBlack);
  H_invariant_mass_GenWideJets_Complete -> SetLineColor(kRed);
  H_invariant_mass_WideJets_Complete       -> SetLineColor(kGreen+1);
  H_invariant_mass_partons_Complete          ->Draw();
  H_invariant_mass_GenWideJets_Complete ->Draw("same");
  H_invariant_mass_WideJets_Complete        ->Draw("same");
  TLegend *Leg1 = new TLegend(0.1, 0.9, 0.5, 0.7);
  Leg1 -> SetHeader("Massa invariante con formula completa:");
  Leg1  ->AddEntry(H_invariant_mass_partons_Complete,"Partoni","L");
  Leg1  ->AddEntry(H_invariant_mass_GenWideJets_Complete,"Gen widejets","L");
  Leg1  ->AddEntry(H_invariant_mass_WideJets_Complete,"Reco widejets","L");
  Leg1->Draw();  
  Canvas ->SaveAs("Invariant_mass_Complete.png");
  
  //formula approssimata
  TCanvas *Canvas1 = new TCanvas("Canvas1","",800,800);
  H_invariant_mass_partons_Approx    -> SetStats(0);
  H_invariant_mass_partons_Approx    -> SetXTitle("M [GeV]");
  H_invariant_mass_partons_Approx    -> SetYTitle("Events");
  H_invariant_mass_partons_Approx    -> GetYaxis()->SetTitleOffset(1.55);
  H_invariant_mass_partons_Approx    -> SetTitle(" ");
  H_invariant_mass_partons_Approx    -> GetXaxis()->SetRangeUser(1500,4000);
  H_invariant_mass_partons_Approx         -> SetLineColor(kBlack);
  H_invariant_mass_GenWideJets_Approx -> SetLineColor(kRed);
  H_invariant_mass_WideJets_Approx       -> SetLineColor(kGreen+1);
  H_invariant_mass_partons_Approx          ->Draw();
  H_invariant_mass_GenWideJets_Approx ->Draw("same");
  H_invariant_mass_WideJets_Approx        ->Draw("same");
  TLegend *Leg2 = new TLegend(0.1, 0.9, 0.5, 0.7);
  Leg2   -> SetHeader("Massa invariante con formula approssimata:");
  Leg2  ->AddEntry(H_invariant_mass_partons_Approx,"Partoni","L");
  Leg2  ->AddEntry(H_invariant_mass_GenWideJets_Approx,"Gen widejets","L");
  Leg2  ->AddEntry(H_invariant_mass_WideJets_Approx,"Reco widejets","L");
  Leg2  ->Draw();  
  Canvas1 ->SaveAs("Invariant_mass_Approx.png");
  }

  ////////////////////////////////////////////////////////////////////////////////

  if(test_parton_VS_Reco == true){
    //confronto
  DrawAndSave("H_invariant_mass_AllGenerator.png" , H_invariant_mass_AllGenerator, 1500, 4000, kBlack, "M [GeV]", "Events");
  DrawAndSave("H_invariant_mass_PartonEnergy_RecoDirection.png" , H_invariant_mass_PartonEnergy_RecoDirection, 1500, 4000, kRed, "M [GeV]", "Events");
  DrawAndSave("H_invariant_mass_PartonDirection_RecoEnergy.png" , H_invariant_mass_PartonDirection_RecoEnergy, 1500, 4000, kGreen+1, "M [GeV]", "Events");
  DrawAndSave("H_invariant_mass_AllReco.png" , H_invariant_mass_AllReco, 1500, 4000, kViolet, "M [GeV]", "Events");

  TCanvas *canvas = new TCanvas("canvas","",800,800);
  H_invariant_mass_AllGenerator                           -> SetStats(0);
  H_invariant_mass_PartonEnergy_RecoDirection  -> SetStats(0);
  H_invariant_mass_PartonDirection_RecoEnergy  -> SetStats(0);
  H_invariant_mass_AllReco                                   -> SetStats(0);
  H_invariant_mass_AllGenerator    -> SetXTitle("M [GeV]");
  H_invariant_mass_AllGenerator    -> SetYTitle("Events");
  H_invariant_mass_AllGenerator    -> GetYaxis()->SetTitleOffset(1.55);
  H_invariant_mass_AllGenerator    -> SetTitle(" ");
  H_invariant_mass_AllGenerator    -> GetXaxis()->SetRangeUser(1500,4000);
  H_invariant_mass_AllGenerator                           -> SetLineColor(kBlack);
  H_invariant_mass_PartonEnergy_RecoDirection  -> SetLineColor(kRed);
  H_invariant_mass_PartonDirection_RecoEnergy  -> SetLineColor(kGreen+1);
  H_invariant_mass_AllReco            -> SetLineColor(kViolet);
  H_invariant_mass_AllGenerator                           ->Draw();
  H_invariant_mass_PartonEnergy_RecoDirection  ->Draw("same");
  H_invariant_mass_PartonDirection_RecoEnergy  ->Draw("same");
  H_invariant_mass_AllReco                                  ->Draw("same");
  TLegend *leg = new TLegend(0.1, 0.9, 0.5, 0.7);
  leg  -> SetHeader("Invariant Mass with:");
  leg  ->AddEntry(H_invariant_mass_AllGenerator,"All generator","L");
  leg  ->AddEntry(H_invariant_mass_PartonEnergy_RecoDirection,"Energy Gen, Direction Reco","L");
  leg  ->AddEntry(H_invariant_mass_PartonDirection_RecoEnergy,"Energy Reco, Direction Gen","L");
  leg  ->AddEntry(H_invariant_mass_AllReco,"All Reco","L");
  leg->Draw();  
  canvas ->SaveAs("Invariant_mass_PartonVSReco_Cases.png");
  }
  ///////////////////////////////////////////////////////


  if(test_Gen_VS_Reco == true){
    //confronto
  DrawAndSave("H_invariant_mass_GenWideJets.png" , H_invariant_mass_GenWideJets, 1500, 4000, kBlack, "M [GeV]", "Events");
  DrawAndSave("H_invariant_mass_GenEnergy_RecoDirection.png" , H_invariant_mass_GenEnergy_RecoDirection, 1500, 4000, kRed, "M [GeV]", "Events");
  DrawAndSave("H_invariant_mass_GenDirection_RecoEnergy.png" , H_invariant_mass_GenDirection_RecoEnergy, 1500, 4000, kGreen+1, "M [GeV]", "Events");
  DrawAndSave("H_invariant_mass_AllReco.png" , H_invariant_mass_AllReco, 1500, 4000, kViolet, "M [GeV]", "Events");

  TCanvas *canvas = new TCanvas("canvas","",800,800);
  H_invariant_mass_GenWideJets                      -> SetStats(0);
  H_invariant_mass_GenEnergy_RecoDirection  -> SetStats(0);
  H_invariant_mass_GenDirection_RecoEnergy  -> SetStats(0);
  H_invariant_mass_AllReco                               -> SetStats(0);

  H_invariant_mass_GenWideJets    -> SetXTitle("M [GeV]");
  H_invariant_mass_GenWideJets    -> SetYTitle("Events");
  H_invariant_mass_GenWideJets    -> GetYaxis()->SetTitleOffset(1.55);
  H_invariant_mass_GenWideJets    -> SetTitle(" ");
  H_invariant_mass_GenWideJets    -> GetXaxis()->SetRangeUser(1500,4000);

  H_invariant_mass_GenWideJets                      -> SetLineColor(kBlack);
  H_invariant_mass_GenEnergy_RecoDirection  -> SetLineColor(kRed);
  H_invariant_mass_GenDirection_RecoEnergy  -> SetLineColor(kGreen+1);
  H_invariant_mass_AllReco                               -> SetLineColor(kViolet);

  H_invariant_mass_GenWideJets                      ->Draw();
  H_invariant_mass_GenEnergy_RecoDirection  ->Draw("same");
  H_invariant_mass_GenDirection_RecoEnergy  ->Draw("same");
  H_invariant_mass_AllReco                               ->Draw("same");
  TLegend *leg = new TLegend(0.1, 0.9, 0.5, 0.7);
  leg  -> SetHeader("Invariant Mass with:");
  leg  ->AddEntry(H_invariant_mass_GenWideJets,"Gen WideJets","L");
  leg  ->AddEntry(H_invariant_mass_GenEnergy_RecoDirection,"Energy Gen, Direction Reco","L");
  leg  ->AddEntry(H_invariant_mass_GenDirection_RecoEnergy,"Energy Reco, Direction Gen","L");
  leg  ->AddEntry(H_invariant_mass_AllReco,"Reco WideJets","L");
  leg->Draw();  
  canvas ->SaveAs("Invariant_mass_GenVSReco_Cases.png");
  }






































   ///////////////////////////////////////////////////////

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
      DrawAndSave(output_Dir, h[ii][jj], 0., 2.,1,"R","Events");
      
    }
  }
  }

}
