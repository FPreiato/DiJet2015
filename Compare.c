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
  //  h1    -> SetStats(0);
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

  cout<<"Integral normalized to: "<<h1->Integral() << endl;
}
////////////////////////


////////////////////////////////////////////

int main(int argc, char* argv[]) {
   
  
  bool Draw_SmearingFunctions;


  TFile file1("SmearingFunctions_AllBinAllSampleQuarkQuarkOddEvents.root");
  TH1D *H_R_QuarkQuark_0_1600 = (TH1D*)file1.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_QuarkQuark_1_2000 = (TH1D*)file1.Get("Histo_R_WideJet_GenWideJet_1_2000");


  TFile file2("SmearingFunctions_AllBinAllSampleGluonGluonOddEvents.root");
  TH1D *H_R_GluonGluon_0_1600 = (TH1D*)file2.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_GluonGluon_1_2000 = (TH1D*)file2.Get("Histo_R_WideJet_GenWideJet_1_2000");

  TFile file3("SmearingFunctions_AllBinAllSampleQstarOddEvents.root");
  TH1D *H_R_Qstar_0_1600 = (TH1D*)file3.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_Qstar_1_2000 = (TH1D*)file3.Get("Histo_R_WideJet_GenWideJet_1_2000");

  TFile file4("SmearingFunctions_AllBinAllSampleSampleQCDOddEvents.root");
  TH1D *H_R_QCD_0_1600 = (TH1D*)file4.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_QCD_1_2000 = (TH1D*)file4.Get("Histo_R_WideJet_GenWideJet_1_2000");
  TH1D *H_QCD_GenWideJet1_Pt = (TH1D*)file4.Get("H_GenWideJet2_Pt");
  TH1D *H_QCD_GenWideJet1_Eta = (TH1D*)file4.Get("H_GenWideJet2_Eta");


  TFile file5("SmearingFunctions_Signal.root");
  TH1D *H_Signal_GenWideJet1_Pt = (TH1D*)file5.Get("H_GenWideJet2_Pt");
  TH1D *H_Signal_GenWideJet1_Eta = (TH1D*)file5.Get("H_GenWideJet2_Eta");



  cout<<"Integral qq: "<<H_R_QuarkQuark_0_1600->Integral()<<endl;
  cout<<"Integral gg: "<<H_R_GluonGluon_0_1600->Integral()<<endl;
  cout<<"Integral gq: "<<H_R_Qstar_0_1600->Integral()<<endl;
  cout<<"Integral QCD: "<<H_R_QCD_0_1600->Integral()<<endl;

  Normalizer(H_R_QuarkQuark_0_1600);
  Normalizer(H_R_GluonGluon_0_1600);
  Normalizer(H_R_QCD_0_1600);
  Normalizer(H_R_Qstar_0_1600);


  cout<<"Integral qq: "<<H_R_QuarkQuark_1_2000->Integral()<<endl;
  cout<<"Integral gg: "<<H_R_GluonGluon_1_2000->Integral()<<endl;
  cout<<"Integral gq: "<<H_R_Qstar_1_2000->Integral()<<endl;
  cout<<"Integral QCD: "<<H_R_QCD_1_2000->Integral()<<endl;

  Normalizer(H_R_QuarkQuark_1_2000);
  Normalizer(H_R_GluonGluon_1_2000);
  Normalizer(H_R_QCD_1_2000);
  Normalizer(H_R_Qstar_1_2000);


  /////////////////////////////////////////////////////////
  //  TGaxis::SetMaxDigits(4);
  TCanvas *Canvas = new TCanvas("Canvas","",800,800);
  
  H_R_GluonGluon_0_1600    -> SetStats(0);
  H_R_GluonGluon_0_1600    -> SetXTitle("M [GeV]");
  H_R_GluonGluon_0_1600    -> SetYTitle("Events");
  H_R_GluonGluon_0_1600    -> GetYaxis()->SetTitleOffset(1.40);
  H_R_GluonGluon_0_1600    -> SetTitle(" ");
  H_R_GluonGluon_0_1600    -> GetXaxis()->SetRangeUser(0.4, 1.3);

  H_R_QCD_0_1600                 -> SetLineColor(kBlue);    
  H_R_QCD_0_1600                 -> SetMarkerColor(kBlue);
  H_R_QCD_0_1600                 -> SetMarkerStyle(2);
 
  H_R_GluonGluon_0_1600     -> SetLineColor(kRed);
  H_R_GluonGluon_0_1600     -> SetMarkerColor(kRed);
  H_R_GluonGluon_0_1600     -> SetMarkerStyle(3);
  
  H_R_Qstar_0_1600                -> SetLineColor(kGreen+1);
  H_R_Qstar_0_1600                -> SetMarkerColor(kGreen+1);
  H_R_Qstar_0_1600                -> SetMarkerStyle(4);
 
  H_R_QuarkQuark_0_1600     -> SetLineColor(kViolet);
  H_R_QuarkQuark_0_1600     -> SetMarkerColor(kViolet);
  H_R_QuarkQuark_0_1600     -> SetMarkerStyle(5);
 

   H_R_GluonGluon_0_1600     ->Draw();   
   H_R_QCD_0_1600                 ->Draw("same");
   H_R_QuarkQuark_0_1600      ->Draw("same");
   H_R_Qstar_0_1600                ->Draw("same");
   
  
   TLegend *Leg1 = new TLegend(0.1, 0.9, 0.4, 0.7);
   Leg1 -> SetHeader("Sample:");
   Leg1  ->AddEntry(H_R_QCD_0_1600,"QCD","PL");
   Leg1  ->AddEntry(H_R_QuarkQuark_0_1600,"QuarkQuark","PL");
   Leg1  ->AddEntry(H_R_Qstar_0_1600,"QuarkGluon","PL");
   Leg1  ->AddEntry(H_R_GluonGluon_0_1600,"GluonGluon","PL");
  
   Leg1->Draw();  
   
   Canvas ->SaveAs("H_R_Confronto_0_1600.png");
   ////////////////////////////////////////////////////////////////

  TCanvas *Canvas_2 = new TCanvas("Canvas_2","",800,800);
  
  H_R_GluonGluon_1_2000    -> SetStats(0);
  H_R_GluonGluon_1_2000    -> SetXTitle("M [GeV]");
  H_R_GluonGluon_1_2000    -> SetYTitle("Events");
  H_R_GluonGluon_1_2000    -> GetYaxis()->SetTitleOffset(1.40);
  H_R_GluonGluon_1_2000    -> SetTitle(" ");
  H_R_GluonGluon_1_2000    -> GetXaxis()->SetRangeUser(0.4, 1.3);

  H_R_QCD_1_2000                 -> SetLineColor(kBlue);    
  H_R_QCD_1_2000                 -> SetMarkerColor(kBlue);
  H_R_QCD_1_2000                 -> SetMarkerStyle(2);

  H_R_GluonGluon_1_2000     -> SetLineColor(kRed);
  H_R_GluonGluon_1_2000     -> SetMarkerColor(kRed);
  H_R_GluonGluon_1_2000     -> SetMarkerStyle(3);
  
  H_R_Qstar_1_2000                -> SetLineColor(kGreen+1);
  H_R_Qstar_1_2000                -> SetMarkerColor(kGreen+1);
  H_R_Qstar_1_2000                -> SetMarkerStyle(4);
 
  H_R_QuarkQuark_1_2000      -> SetLineColor(kViolet);
  H_R_QuarkQuark_1_2000      -> SetMarkerColor(kViolet);
  H_R_QuarkQuark_1_2000      -> SetMarkerStyle(5);
 

   H_R_GluonGluon_1_2000     ->Draw();   
   H_R_QCD_1_2000                 ->Draw("same");
   H_R_QuarkQuark_1_2000      ->Draw("same");
   H_R_Qstar_1_2000                ->Draw("same");
   
  
   TLegend *Leg2 = new TLegend(0.1, 0.9, 0.4, 0.7);
   Leg2 -> SetHeader("Sample:");
   Leg2  ->AddEntry(H_R_QCD_1_2000,"QCD","PL");
   Leg2  ->AddEntry(H_R_QuarkQuark_1_2000,"QuarkQuark","PL");
   Leg2  ->AddEntry(H_R_Qstar_1_2000,"QuarkGluon","PL");
   Leg2  ->AddEntry(H_R_GluonGluon_1_2000,"GluonGluon","PL");
  
   Leg2->Draw();  
   
   Canvas_2 ->SaveAs("H_R_Confronto_1_2000.png");

   /////////////////////////////////////////////////////////////////

  TCanvas *Canvas_3 = new TCanvas("Canvas_3","",800,800);
  H_QCD_GenWideJet1_Pt    -> SetStats(0);
  H_QCD_GenWideJet1_Pt    -> SetXTitle("Pt [GeV]");
  H_QCD_GenWideJet1_Pt    -> SetYTitle("Normalized");
  H_QCD_GenWideJet1_Pt    -> GetYaxis()->SetTitleOffset(1.55);
  H_QCD_GenWideJet1_Pt    -> SetTitle(" ");
  H_QCD_GenWideJet1_Pt    -> GetXaxis()->SetRangeUser(0, 5000);
  H_QCD_GenWideJet1_Pt                 -> SetLineColor(kBlue);    
  H_QCD_GenWideJet1_Pt                 -> SetMarkerColor(kBlue);
  H_QCD_GenWideJet1_Pt                 -> SetMarkerStyle(2);
  H_Signal_GenWideJet1_Pt     -> SetLineColor(kRed);
  H_Signal_GenWideJet1_Pt     -> SetMarkerColor(kRed);
  H_Signal_GenWideJet1_Pt     -> SetMarkerStyle(3);
   H_QCD_GenWideJet1_Pt       ->DrawNormalized();
   H_Signal_GenWideJet1_Pt     ->DrawNormalized("same");   
   TLegend *Leg3 = new TLegend(0.9, 0.9, 0.7, 0.7);
   Leg3 -> SetHeader("Sample:");
   Leg3  ->AddEntry(H_QCD_GenWideJet1_Pt,"QCD","PL");
   Leg3  ->AddEntry(H_Signal_GenWideJet1_Pt,"Signal","PL");
   Leg3->Draw();  
   Canvas_3 ->SaveAs("H_SignalvsQCD_GenWideJet2_Pt.png");
   ////////////////////////////////////////////////////////////////

  TCanvas *Canvas_4 = new TCanvas("Canvas_4","",800,800);
  H_QCD_GenWideJet1_Eta    -> SetStats(0);
  H_QCD_GenWideJet1_Eta    -> SetXTitle("#eta");
  H_QCD_GenWideJet1_Eta    -> SetYTitle("Normalized");
  H_QCD_GenWideJet1_Eta    -> GetYaxis()->SetTitleOffset(1.55);
  H_QCD_GenWideJet1_Eta    -> SetTitle(" ");
  H_QCD_GenWideJet1_Eta    -> GetXaxis()->SetRangeUser(-3, 3);
  H_QCD_GenWideJet1_Eta                 -> SetLineColor(kBlue);    
  H_QCD_GenWideJet1_Eta                 -> SetMarkerColor(kBlue);
  H_QCD_GenWideJet1_Eta                 -> SetMarkerStyle(2);
  H_Signal_GenWideJet1_Eta     -> SetLineColor(kRed);
  H_Signal_GenWideJet1_Eta     -> SetMarkerColor(kRed);
  H_Signal_GenWideJet1_Eta     -> SetMarkerStyle(3);
   H_QCD_GenWideJet1_Eta       ->DrawNormalized();
   H_Signal_GenWideJet1_Eta     ->DrawNormalized("same");   
   //   TLegend *Leg3 = new TLegend(0.9, 0.9, 0.7, 0.7);
   //Leg3 -> SetHeader("Sample:");
   //Leg3  ->AddEntry(H_QCD_GenWideJet1_Eta,"QCD","PL");
   //Leg3  ->AddEntry(H_Signal_GenWideJet1_Eta,"Signal","PL");
   Leg3->Draw();  
   Canvas_4 ->SaveAs("H_SignalvsQCD_GenWideJet2_Eta.png");
   ////////////////////////////////////////////////////////////////


















   ///////////////////////////////////////////////////////

  if(Draw_SmearingFunctions == true){
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
