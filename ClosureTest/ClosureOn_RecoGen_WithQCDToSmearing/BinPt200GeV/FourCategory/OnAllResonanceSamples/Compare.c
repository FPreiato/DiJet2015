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

#include "/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/Compare.h"

using namespace std;

int main(int argc, char* argv[]) {

  string * inputFile = new string(argv[1]);
  string * output_dir = new string(argv[2]);

  char * input_file = new char [inputFile->length()+1];
  strcpy (input_file, inputFile->c_str());

  char * directory_output = new char [output_dir->length()+1];
  strcpy (directory_output, output_dir->c_str());

  //  TFile file1("ClosureOn_RecoGen_QCDSamples/rootfile_RSGravitonToQuarkQuark_kMpl01_M_3000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM.root");

  TFile file1(input_file);

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
  //  H_DijetSmeared_M->Rebin(5);
   H_DijetSmeared_M->Rebin(10);

  TH1D *H_dijetWide_Pt   = (TH1D*)file1.Get("H_dijetWide_Pt");
  TH1D *H_dijetWide_Eta = (TH1D*)file1.Get("H_dijetWide_Eta");
  TH1D *H_dijetWide_Phi = (TH1D*)file1.Get("H_dijetWide_Phi");
  TH1D *H_dijetWide_M   = (TH1D*)file1.Get("H_dijetWide_M");
  H_dijetWide_Pt->Rebin(5);
  H_dijetWide_Eta->Rebin(5);
  H_dijetWide_Phi->Rebin(5);
  // H_dijetWide_M->Rebin(5);
  H_dijetWide_M->Rebin(10);

  cout<<"File exist"<<endl;

  // Normalizer(H_JetSmeared1_Pt);
  // Normalizer(H_WideJet1_Pt);
  // Normalizer(H_JetSmeared2_Pt);
  // Normalizer(H_WideJet2_Pt);
  // Normalizer(H_DijetSmeared_M);
  // Normalizer(H_dijetWide_M);
    
 

  ///////////////////////////////////////////////////////////////////////////////////

  bool Comparison = true;

  if(Comparison == true){

    // TLegend *leg2 = new TLegend( 0.15, 0.9, 0.45, 0.6);
    // TLegend *leg2 = new TLegend( 0.9, 0.9, 0.6, 0.6);
    // leg2 -> SetBorderSize( 0);
    // leg2 -> SetFillColor ( 0);
    // leg2 -> SetFillStyle ( 0);
    // leg2 -> SetTextFont ( 42);
    // leg2 -> SetTextSize (0.035);
    // leg2 -> SetHeader("  ");
    // leg2 -> AddEntry(H_dijetWide_M, "Reco WideJet", "L");
    // leg2 -> AddEntry(H_DijetSmeared_M, "Jet Smeared ", "L");

  int xmin = 0;
  int xmax = 10000;


  DrawRatioAndSave(directory_output, "Compare_Dijet_M.png",    H_dijetWide_M,     H_DijetSmeared_M,    xmin, xmax, "M [GeV]", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Dijet_Pt.png",    H_dijetWide_Pt,     H_DijetSmeared_Pt,    0, 3000, "Pt [GeV]", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Dijet_Eta.png",  H_dijetWide_Eta,   H_DijetSmeared_Eta, -4, 4, "#eta", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Dijet_Phi.png",  H_dijetWide_Phi,   H_DijetSmeared_Phi, -4, 4, "#phi", "Events", true) ;

  DrawRatioAndSave(directory_output, "Compare_Pt_Jet1.png",   H_WideJet1_Pt,   H_JetSmeared1_Pt,    0, 3000, "Pt [GeV]", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Eta_Jet1.png", H_WideJet1_Eta , H_JetSmeared1_Eta, -4, 4, "#eta", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Phi_Jet1.png", H_WideJet1_Phi , H_JetSmeared1_Phi, -4, 4, "#phi", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_M_Jet1.png",   H_WideJet1_M ,   H_JetSmeared1_M,    0, 1000, "M [GeV]", "Events", true) ;
    
  DrawRatioAndSave(directory_output, "Compare_Pt_Jet2.png",   H_WideJet2_Pt,   H_JetSmeared2_Pt,    0, 3000, "Pt [GeV]", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Eta_Jet2.png", H_WideJet2_Eta , H_JetSmeared2_Eta, -4, 4, "#eta", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_Phi_Jet2.png", H_WideJet2_Phi , H_JetSmeared2_Phi, -4, 4, "#phi", "Events", true) ;
  DrawRatioAndSave(directory_output, "Compare_M_Jet2.png",   H_WideJet2_M ,   H_JetSmeared2_M,    0, 1000, "M [GeV]", "Events", true) ;
    
  DrawPullAndSave(directory_output, "Histo_WithPull.png", H_dijetWide_M, H_DijetSmeared_M, xmin, xmax, "M [GeV]", "Events") ;

  }else{
    cout<<"Non sto facendo confronti"<<endl;
  }

}






