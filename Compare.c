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

#include "/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/Compare.h"

using namespace std;

///////////////////////////


int main(int argc, char* argv[]) {
  
  bool Draw_SmearingFunctions = false;
  
  char directory_output[300] ="SmearingFunctions_RecoGen_QCDSample/" ;
  
  TFile file1("SmearingFunctions_RecoGen_QCDSample/SmearingFunctions_AllBinAllSampleQuarkQuarkOddEvents.root");
  TH1D *H_R_QuarkQuark_0_1600 = (TH1D*)file1.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_QuarkQuark_1_2000 = (TH1D*)file1.Get("Histo_R_WideJet_GenWideJet_1_2000");

  TFile file2("SmearingFunctions_RecoGen_QCDSample/SmearingFunctions_AllBinAllSampleGluonGluonOddEvents.root");
  TH1D *H_R_GluonGluon_0_1600 = (TH1D*)file2.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_GluonGluon_1_2000 = (TH1D*)file2.Get("Histo_R_WideJet_GenWideJet_1_2000");

  TFile file3("SmearingFunctions_RecoGen_QCDSample/SmearingFunctions_AllBinAllSampleQstarOddEvents.root");
  TH1D *H_R_Qstar_0_1600 = (TH1D*)file3.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_Qstar_1_2000 = (TH1D*)file3.Get("Histo_R_WideJet_GenWideJet_1_2000");

  TFile file4("SmearingFunctions_RecoGen_QCDSample/SmearingFunctions_AllBinAllSampleSampleQCDOddEvents.root");
  TH1D *H_R_QCD_0_1600 = (TH1D*)file4.Get("Histo_R_WideJet_GenWideJet_0_1600");
  TH1D *H_R_QCD_1_2000 = (TH1D*)file4.Get("Histo_R_WideJet_GenWideJet_1_2000");
  TH1D *H_QCD_GenWideJet1_Pt = (TH1D*)file4.Get("H_GenWideJet1_Pt");
  TH1D *H_QCD_GenWideJet1_Eta = (TH1D*)file4.Get("H_GenWideJet1_Eta");

  TFile file5("SmearingFunctions_RecoGen_QCDSample/SmearingFunctions_Signal.root");
  TH1D *H_Signal_GenWideJet1_Pt = (TH1D*)file5.Get("H_GenWideJet1_Pt");
  TH1D *H_Signal_GenWideJet1_Eta = (TH1D*)file5.Get("H_GenWideJet1_Eta");

  // cout<<"Integral qq: "<<H_R_QuarkQuark_0_1600->Integral()<<endl;
  // cout<<"Integral gg: "<<H_R_GluonGluon_0_1600->Integral()<<endl;
  // cout<<"Integral gq: "<<H_R_Qstar_0_1600->Integral()<<endl;
  // cout<<"Integral QCD: "<<H_R_QCD_0_1600->Integral()<<endl;

  cout<<"Normalizing.."<<endl;
  Normalizer(H_R_QuarkQuark_0_1600);
  Normalizer(H_R_GluonGluon_0_1600);
  Normalizer(H_R_QCD_0_1600);
  Normalizer(H_R_Qstar_0_1600);
  Normalizer(H_R_QuarkQuark_1_2000);
  Normalizer(H_R_GluonGluon_1_2000);
  Normalizer(H_R_QCD_1_2000);
  Normalizer(H_R_Qstar_1_2000);
  Normalizer(H_Signal_GenWideJet1_Pt);
  Normalizer(H_Signal_GenWideJet1_Eta);


  TLegend *Leg1 = new TLegend(0.1, 0.9, 0.4, 0.7);
  Leg1 -> SetHeader("Sample:");
  Leg1  ->AddEntry(H_R_QCD_0_1600,"QCD","PL");
  Leg1  ->AddEntry(H_R_QuarkQuark_0_1600,"QuarkQuark","PL");
  Leg1  ->AddEntry(H_R_Qstar_0_1600,"QuarkGluon","PL");
  Leg1  ->AddEntry(H_R_GluonGluon_0_1600,"GluonGluon","PL");
 
  DrawAndSave(directory_output, "H_R_Confronto_0_1600.png", H_R_QCD_0_1600, H_R_QuarkQuark_0_1600, H_R_Qstar_0_1600, H_R_GluonGluon_0_1600, 0.4, 1.8, "M [GeV]", "Events", Leg1);
    DrawAndSave(directory_output, "H_R_Confronto_1_2000.png", H_R_QCD_1_2000, H_R_QuarkQuark_1_2000, H_R_Qstar_1_2000, H_R_GluonGluon_1_2000, 0.4, 1.8, "M [GeV]", "Events", Leg1);

  TLegend *Leg3 = new TLegend(0.9, 0.9, 0.7, 0.7);
  Leg3 -> SetHeader("Sample:");
  Leg3  ->AddEntry(H_QCD_GenWideJet1_Pt,"QCD","PL");
  Leg3  ->AddEntry(H_Signal_GenWideJet1_Pt,"Signal","PL");

  DrawAndSave(directory_output, "H_SignalvsQCD_GenWideJet1_Pt.png", H_QCD_GenWideJet1_Pt, H_Signal_GenWideJet1_Pt, 0, 5000, "Pt [GeV]", "Normalized", Leg3);
  DrawAndSave(directory_output, "H_SignalvsQCD_GenWideJet1_Eta.png", H_QCD_GenWideJet1_Eta, H_Signal_GenWideJet1_Eta, -3, 3, "#eta", "Normalized", Leg3);

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

	char output_Dir[200] = "Histogram_SmearingFunctions";
	
	char extension[10] = ".png";
	strcat( HistoName, extension);
	
	cout<<HistoName<<endl;
	DrawAndSave(output_Dir, HistoName, h[ii][jj], 0., 2.,"R","Events");
      
      }
    }
  }

}
