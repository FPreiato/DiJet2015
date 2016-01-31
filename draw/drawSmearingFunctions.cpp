#include <stdlib.h>
#include <sys/stat.h>

#include "TParameter.h"
#include "TPaveText.h"
#include "TError.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "drawSmearingFunctions.h"

#include "etaBinning.h"
#include "ptBinning.h"

#include <TColor.h>


int main(int argc, char* argv[]) {

  if (argc != 2 && argc != 3) {

    std::cout << "USAGE: ./drawPhotonJet [data_dataset]" << std::endl;
    exit(23);
  }

  std::string outputDir = "SmearingFunction";
  mkdir(outputDir.c_str(), 0777);

  std::string data_dataset(argv[1]);
  TString dataFileName;
  dataFileName = TString::Format("%s", data_dataset.c_str() );
  TFile* dataFile = TFile::Open(dataFileName);
  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
  }

  int n_cat;  
  TParameter<double>* cat = static_cast<TParameter<double>*>(dataFile->Get("n_categories"));                                                    
  n_cat = cat->GetVal();                                                                                    
  int n_files;  
  TParameter<double>* file = static_cast<TParameter<double>*>(dataFile->Get("n_files"));                              
  n_files = file->GetVal();

  int n_categories;  
  n_categories = n_cat / n_files;

  std::cout << "N categories choose "<< n_categories << std::endl; 
  
  gErrorIgnoreLevel = kWarning;
  
  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();
  PtBinning ptBinning;  
  size_t ptBinningSize = ptBinning.size();

  for(int kk =0; kk<n_categories; kk++){
    
    for (size_t i = 0; i < etaBinningSize; i++) {
      const std::string etaName = etaBinning.getBinName(i);
      const std::pair<float, float> eta_bin = etaBinning.getBinValue(i);
      
      for (size_t j = 0; j < ptBinningSize; j++) {
	
	const std::pair<float, float> pt_bin = ptBinning.getBinValue(j);
	std::stringstream ss;
	ss << "pt_" << (int) pt_bin.first << "_" << (int) pt_bin.second;
	
	TString responseName;
	if(n_categories == 2){
	  if(kk == 0)  responseName = TString::Format("smearing_Quark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 1)  responseName = TString::Format("smearing_Gluon_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	} else if(n_categories == 3){
	  if(kk == 0)  responseName = TString::Format("smearing_LightQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 1)  responseName = TString::Format("smearing_HeavyQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 2)  responseName = TString::Format("smearing_Gluon_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	}else if(n_categories == 4){
	  if(kk == 0)  responseName = TString::Format("smearing_LightQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 1)  responseName = TString::Format("smearing_CharmQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 2)  responseName = TString::Format("smearing_BottomQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 3)  responseName = TString::Format("smearing_Gluon_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	}else if (n_categories == 6){
	  if(kk == 0)  responseName = TString::Format("smearing_UpQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 1)  responseName = TString::Format("smearing_DownQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 2)  responseName = TString::Format("smearing_StrangeQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 3)  responseName = TString::Format("smearing_CharmQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 4)  responseName = TString::Format("smearing_BottomQuark_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	  if(kk == 5)  responseName = TString::Format("smearing_Gluon_%s_%s", etaBinning.getBinName(i).c_str(), ss.str().c_str() );
	}
	
	TH1F *h = (TH1F*)dataFile->Get("smearingFunction/"+responseName);	 
	std::cout<< responseName << std::endl;
	
	Float_t label_cuts_xMin  = 0.66;
	Float_t label_cuts_yMin  = 0.7;
	Float_t label_cuts_xMax =0.95;
	Float_t label_cuts_yMax = 0.9;
	
	TPaveText* range = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax, "nbNDC");
	range->SetFillColor(kWhite);
	range->SetTextSize(0.035);
	range->SetTextFont(42);
	char etaRange[100];
	
	sprintf(etaRange, " %.1f < |#eta| < %.1f", eta_bin.first, eta_bin.second);
	range->AddText(etaRange);
	char ptRange[70];
	sprintf(ptRange, " %.0f < p_{T}^{reco} < %.0f GeV", pt_bin.first , pt_bin.second);
	range->AddText(ptRange);
	   
	DrawAndSave("SmearingFunction/", responseName+".png", h, 0, 2, "Ratio" , "Events", range);	   
      }
    }   
  }
  return 0;
  
}
