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
  
  bool Draw_SmearingFunctions = true;
  

  if(Draw_SmearingFunctions == true){

    TFile file1("SmearingFunctions_RecoGen_ResonanceSamples_Bin100GeV/SmearingFunctions_ResonanceSamplesToQuarkQuark_Bin100.root ");
    TH1D *H_step_pt = (TH1D*)file1.Get("H_step_pt");
    int step_pt = H_step_pt->GetMean();
    TFile file2("SmearingFunctions_RecoGen_QCDSamples_Bin100GeV/SmearingFunctions_QCDSamplesToQuarkQuark_Bin100.root ");    

    int n_bin = 4500 / step_pt ;   
    char HistoName[200];	 
    TH1D *h[100][100];
    TH1D *j[100][100];

    TLegend *leg = new TLegend( 0.9, 0.9, 0.6, 0.6);
    leg -> SetHeader("Smearing Functions: ");
    leg -> AddEntry(h[0][0], "Resonance Samples", "L");
    leg -> AddEntry(j[0][0], "QCD Samples", "L");


    for(int ii=0; ii<5; ii++){
      for(int jj=0; jj<n_bin ; jj++){	 	 
	
	int pt_bin_max = (jj*step_pt)+step_pt;

	sprintf(HistoName,"Histo_RQuarks_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	j[ii][jj]= (TH1D*)file2.Get(HistoName);   


	Normalizer(h[ii][jj]);
	Normalizer(j[ii][jj]);

	//h[ii][jj]->Integral();
	//cout<< "ii: "<< ii<< "   jj: "<<jj<<endl;
	cout<<"      h[ii][jj]->Integral(): "<<       h[ii][jj]->Integral() <<endl;
	cout<<"      j[ii][jj]->Integral(): "<<       j[ii][jj]->Integral() <<endl;

	char output_Dir[200] = "Histogram_SmearingFunctions/";
	char extension[10] = ".png";
	strcat( HistoName, extension);
	cout<<HistoName<<endl;

	if( j[ii][jj]->Integral() >= h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, j[ii][jj], h[ii][jj], 0., 2.,"R","Normalized", kBlue, kRed);
	if( j[ii][jj]->Integral() < h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, h[ii][jj], j[ii][jj], 0., 2.,"R","Normalized", kRed, kBlue);

      }
    }
  }

  /*
  TFile file10("output/rootFile.root");
  TH1D *H_step_pt = (TH1D*)file10.Get("H_step_pt");
  int step_pt = H_step_pt->GetMean();
  cout<< step_pt<<endl;
  
  int n_bin = 4500 / step_pt ;   
  char HistoName[200];	 
  TH1D *h;
  double sum_integral;  

  for(int ii=0; ii<5; ii++){
    for(int jj=0; jj<n_bin ; jj++){	 	 
      
      int pt_bin_max = (jj*step_pt)+step_pt;
      sprintf(HistoName,"Histo_RQuarks_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
      h= (TH1D*)file10.Get(HistoName);   
      h->Integral();
      //cout<< "ii: "<< ii<< "   jj: "<<jj<<endl;
      cout<<"      Integral(): "<< HistoName<<" : "<<h->Integral() <<endl;
      
      sum_integral += h->Integral();
    }
  }

  cout<<"Sum Integral: "<< sum_integral<<endl;

  double sum_integral_g;  
  for(int ii=0; ii<5; ii++){
    for(int jj=0; jj<n_bin ; jj++){	 	 
      
      int pt_bin_max = (jj*step_pt)+step_pt;
      sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
      h= (TH1D*)file10.Get(HistoName);   
      h->Integral();
      //cout<< "ii: "<< ii<< "   jj: "<<jj<<endl;
      cout<<"      Integral(): "<< HistoName<<" : "<<h->Integral() <<endl;
      
      sum_integral_g += h->Integral();
    }
  }

  cout<<"Sum Integral_g: "<< sum_integral_g<<endl;
*/



}
