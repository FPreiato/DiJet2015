#include <stdlib.h>
#include <sys/stat.h>  

#include "TParameter.h"
#include "TPaveText.h"
#include "TError.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "drawClosureTest.h"
//#include "Compare.h"


#include "etaBinning.h"
#include "ptBinning.h"

#include <TColor.h>

int main(int argc, char* argv[]) {

  if (argc != 3 && argc != 4) {

    std::cout << "USAGE: ./drawPhotonJet [dataset_StandardAnalysis] [dataset_ClosureTest]" << std::endl;
    exit(23);
  }

  std::string outputDir = "PlotComparison";                 
  mkdir(outputDir.c_str(), 0777); 

  //truth
  std::string data1_dataset(argv[1]);
  //smeared
  std::string data2_dataset(argv[2]);

  TString dataFileName1;
  dataFileName1 = TString::Format("%s", data1_dataset.c_str() );
  TString dataFileName2;
  dataFileName2 = TString::Format("%s", data2_dataset.c_str() );

  TFile* dataFile1 = TFile::Open(dataFileName1);
  TFile* dataFile2 = TFile::Open(dataFileName2);

  if (dataFile1) {
    std::cout << "Opened data file '" << dataFileName1 << "'." << std::endl;
  }
  if (dataFile2) {
    std::cout << "Opened data file '" << dataFileName2 << "'." << std::endl;
  }

  //  bool log = true;
  gErrorIgnoreLevel = kWarning;
  
  //  TString HistoName = TString::Format("smearing_Quark");
  // si puo fare un vector di stringhe?
  //cosi passo il vettore

  vector<TString> HistoName;
  HistoName.push_back(TString::Format("cutHisto_noCuts_____dijetWide_M"));
  HistoName.push_back(TString::Format("cutHisto_noCuts_____dijetWide_pT"));
  HistoName.push_back(TString::Format("cutHisto_noCuts_____WideJet1_pT"));
  HistoName.push_back(TString::Format("cutHisto_noCuts_____WideJet1_Eta"));
  HistoName.push_back(TString::Format("cutHisto_noCuts_____WideJet2_pT"));
  HistoName.push_back(TString::Format("cutHisto_noCuts_____WideJet2_Eta"));
  HistoName.push_back(TString::Format("cutHisto_noCuts_____nVtx"));
 

  vector<TString> XAxis;
  XAxis.push_back(TString::Format("M(jj) [GeV]"));
  XAxis.push_back(TString::Format("pT(jj) [GeV]"));
  XAxis.push_back(TString::Format("pT(WideJet1) [GeV]"));
  XAxis.push_back(TString::Format("#eta(WideJet1) [GeV]"));
  XAxis.push_back(TString::Format("pT(WideJet2) [GeV]"));
  XAxis.push_back(TString::Format("#eta(WideJet2) [GeV]"));
  XAxis.push_back(TString::Format("n vertex"));
  
  
  size_t size = HistoName.size();
  std::cout<< size << std::endl;

  for(size_t ii = 0; ii < size; ii++){
    std::cout<< HistoName.at(ii) << std::endl;
    
    TH1F *h1 = (TH1F*)dataFile1->Get(HistoName.at(ii) );
    TH1F *h2 = (TH1F*)dataFile2->Get(HistoName.at(ii) );
    double mean_h1 = h1->GetMean();
    double mean_h2 = h2->GetMean();

    std::cout<< mean_h1 << std::endl; 
    std::cout<< mean_h2 << std::endl;

    h1->Rebin(10);
    h2->Rebin(10);

    Normalizer(h1);
    Normalizer(h2);

    TLegend* legend = new TLegend(0.70, 0.80, 0.90, 0.90);
    legend->SetTextFont(42);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.036);
    legend->AddEntry(h1, "Standard An.", "PL");
    legend->AddEntry(h2, "Closure", "PL");
   	 
    DrawPullAndSave("PlotComparison/", HistoName.at(ii)+"_Pull.png", h1, h2, XAxis.at(ii) , "Events", legend);
    DrawRatioAndSave("PlotComparison/", HistoName.at(ii)+"_Ratio.png", h1, h2, XAxis.at(ii) , "Events", legend);
    
  }
  
  
  
  return 0;
  
}

