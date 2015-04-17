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
#include <TF1.h>

using namespace std;

///////////////////////////

void DrawAndSave(const char *namefile ,TH1D* h1, int xmin, int xmax){
   
  //  cout<<"working"<<endl;
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetLineColor(kBlue);
  h1->SetXTitle("R");
  h1->Draw();
  
  canvas ->SaveAs(namefile);
  canvas ->Destructor();
}


//////////////////////////
void DrawAndSave(const char *namefile ,TH1D* h1, TH1D* h2, int xmin, int xmax){
   
  //  cout<<"working"<<endl;
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->SetStats(0);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->Draw();
  h2->Draw("same");

  TLegend *leg1 = new TLegend( 0.9, 0.9, 0.65, 0.7);
  //leg1 -> SetHeader("Pt:");
  leg1 -> AddEntry(h1, "q^{*} -> g q (pdgId =1)", "L");
  leg1 -> AddEntry(h2, "q^{*} -> g q (pdgId =2)", "L");
  leg1->Draw();



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
  h1 -> SetYTitle("Normalized to unit");
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
   
  
  TFile file1("/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/TestF/QStar_M3000_pdgId_1.root");

  TH1D *H_Qstar_M3000_pdgId_1_WideJet1_Eta   = (TH1D*)file1.Get("H_WideJet1_Eta");
  H_Qstar_M3000_pdgId_1_WideJet1_Eta ->Rebin(5); 
  TH1D *H_Qstar_M3000_pdgId_1_WideJet2_Eta   = (TH1D*)file1.Get("H_WideJet2_Eta");
  H_Qstar_M3000_pdgId_1_WideJet2_Eta ->Rebin(5); 
  TH1D *H_Qstar_M3000_pdgId_1_DEtaJJ   = (TH1D*)file1.Get("H_DeltaEtaJJWide");
  H_Qstar_M3000_pdgId_1_DEtaJJ ->Rebin(5); 
  TH1D *H_Qstar_M3000_pdgId_1_dijetWide_M = (TH1D*)file1.Get("H_dijetWide_M");
  H_Qstar_M3000_pdgId_1_dijetWide_M->Rebin(5);

  TFile file2("/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/TestF/QStar_M5000_pdgId_1.root");

  TH1D *H_Qstar_M5000_pdgId_1_WideJet1_Eta   = (TH1D*)file2.Get("H_WideJet1_Eta");
  H_Qstar_M5000_pdgId_1_WideJet1_Eta ->Rebin(5); 
  TH1D *H_Qstar_M5000_pdgId_1_WideJet2_Eta   = (TH1D*)file2.Get("H_WideJet2_Eta");
  H_Qstar_M5000_pdgId_1_WideJet2_Eta ->Rebin(5); 
  TH1D *H_Qstar_M5000_pdgId_1_DEtaJJ   = (TH1D*)file2.Get("H_DeltaEtaJJWide");
  H_Qstar_M5000_pdgId_1_DEtaJJ ->Rebin(5); 
  TH1D *H_Qstar_M5000_pdgId_1_dijetWide_M = (TH1D*)file2.Get("H_dijetWide_M");
  H_Qstar_M5000_pdgId_1_dijetWide_M->Rebin(5);


  TFile file3("/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/TestF/QStar_M3000_pdgId_2.root");

  TH1D *H_Qstar_M3000_pdgId_2_WideJet1_Eta   = (TH1D*)file3.Get("H_WideJet1_Eta");
  H_Qstar_M3000_pdgId_2_WideJet1_Eta ->Rebin(5); 
  TH1D *H_Qstar_M3000_pdgId_2_WideJet2_Eta   = (TH1D*)file3.Get("H_WideJet2_Eta");
  H_Qstar_M3000_pdgId_2_WideJet2_Eta ->Rebin(5); 
  TH1D *H_Qstar_M3000_pdgId_2_DEtaJJ   = (TH1D*)file3.Get("H_DeltaEtaJJWide");
  H_Qstar_M3000_pdgId_2_DEtaJJ ->Rebin(5); 
  TH1D *H_Qstar_M3000_pdgId_2_dijetWide_M = (TH1D*)file3.Get("H_dijetWide_M");
  H_Qstar_M3000_pdgId_2_dijetWide_M->Rebin(5);

  TFile file4("/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/TestF/QStar_M5000_pdgId_2.root");

  TH1D *H_Qstar_M5000_pdgId_2_WideJet1_Eta   = (TH1D*)file4.Get("H_WideJet1_Eta");
  H_Qstar_M5000_pdgId_2_WideJet1_Eta ->Rebin(5); 
  TH1D *H_Qstar_M5000_pdgId_2_WideJet2_Eta   = (TH1D*)file4.Get("H_WideJet2_Eta");
  H_Qstar_M5000_pdgId_2_WideJet2_Eta ->Rebin(5); 
  TH1D *H_Qstar_M5000_pdgId_2_DEtaJJ   = (TH1D*)file4.Get("H_DeltaEtaJJWide");
  H_Qstar_M5000_pdgId_2_DEtaJJ ->Rebin(5); 
  TH1D *H_Qstar_M5000_pdgId_2_dijetWide_M = (TH1D*)file4.Get("H_dijetWide_M");
  H_Qstar_M5000_pdgId_2_dijetWide_M->Rebin(5);


  ///////////////////////////////////////////////////////////////////////////////////////

  cout<<"File exist"<<endl;

  Normalizer(H_Qstar_M3000_pdgId_1_WideJet1_Eta);
  Normalizer(H_Qstar_M3000_pdgId_1_WideJet2_Eta);
  Normalizer(H_Qstar_M3000_pdgId_1_DEtaJJ);
  Normalizer(H_Qstar_M3000_pdgId_1_dijetWide_M);

  Normalizer(H_Qstar_M3000_pdgId_2_WideJet1_Eta);
  Normalizer(H_Qstar_M3000_pdgId_2_WideJet2_Eta);
  Normalizer(H_Qstar_M3000_pdgId_2_DEtaJJ);
  Normalizer(H_Qstar_M3000_pdgId_2_dijetWide_M);

  Normalizer(H_Qstar_M5000_pdgId_1_WideJet1_Eta);
  Normalizer(H_Qstar_M5000_pdgId_1_WideJet2_Eta);
  Normalizer(H_Qstar_M5000_pdgId_1_DEtaJJ);
  Normalizer(H_Qstar_M5000_pdgId_1_dijetWide_M);

  Normalizer(H_Qstar_M5000_pdgId_2_WideJet1_Eta);
  Normalizer(H_Qstar_M5000_pdgId_2_WideJet2_Eta);
  Normalizer(H_Qstar_M5000_pdgId_2_DEtaJJ);
  Normalizer(H_Qstar_M5000_pdgId_2_dijetWide_M);

  bool Comparison = true;

  if(Comparison == true){

    DrawAndSave("Qstar_M3000_Jet1_Eta.png", H_Qstar_M3000_pdgId_1_WideJet1_Eta, H_Qstar_M3000_pdgId_2_WideJet1_Eta, -4., 4.);
    DrawAndSave("Qstar_M3000_Jet2_Eta.png", H_Qstar_M3000_pdgId_1_WideJet2_Eta, H_Qstar_M3000_pdgId_2_WideJet2_Eta, -4., 4.);
    DrawAndSave("Qstar_M3000_DEtaJJ.png",    H_Qstar_M3000_pdgId_1_DEtaJJ,            H_Qstar_M3000_pdgId_2_DEtaJJ, 0., 9.);
    DrawAndSave("Qstar_M3000_MJJ.png",        H_Qstar_M3000_pdgId_1_dijetWide_M,   H_Qstar_M3000_pdgId_2_dijetWide_M, 0,5000);    

    DrawAndSave("Qstar_M5000_Jet1_Eta.png", H_Qstar_M5000_pdgId_1_WideJet1_Eta, H_Qstar_M5000_pdgId_2_WideJet1_Eta, -4., 4.);
    DrawAndSave("Qstar_M5000_Jet2_Eta.png", H_Qstar_M5000_pdgId_1_WideJet2_Eta, H_Qstar_M5000_pdgId_2_WideJet2_Eta, -4., 4.);
    DrawAndSave("Qstar_M5000_DEtaJJ.png",    H_Qstar_M5000_pdgId_1_DEtaJJ,            H_Qstar_M5000_pdgId_2_DEtaJJ, 0., 9.);
    DrawAndSave("Qstar_M5000_MJJ.png",        H_Qstar_M5000_pdgId_1_dijetWide_M,   H_Qstar_M5000_pdgId_2_dijetWide_M, 0,8000);    



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
