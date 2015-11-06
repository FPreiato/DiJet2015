#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TMath.h>
#include "ptBinning.h"
#include "etaBinning.h"

bool verbose = false;

// Reading file with Smearing Functions
TFile file1("/cmshome/fpreiato/DiJet/test/CMSSW_7_4_14/src/CMSDIJET/DijetRootTreeAnalyzer/output/rootFile_smearing.root");
int n_categories;

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;
  
  string jetAlgo = getPreCutString1("jetAlgo");
  double rParam = getPreCutValue1("DeltaR");
  
  if( jetAlgo == "AntiKt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::antikt_algorithm, rParam) );
  else if( jetAlgo == "Kt" )
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::kt_algorithm, rParam) );
  else 
    fjJetDefinition = JetDefPtr( new fastjet::JetDefinition(fastjet::cambridge_algorithm, rParam) );
  
  TParameter<double>* cat = static_cast<TParameter<double>*>(file1.Get("n_categories"));
  n_categories = cat->GetVal();
  std::cout << "N categories choose "<< n_categories << std::endl;  
  
  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;
  
  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}


void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop() begins" <<std::endl;   
  
  if (fChain == 0) return;
  
  /////////initialize variables
  
  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   
  
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //for (Long64_t jentry=0; jentry<1000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code starts here ///////////////////////
    ///Stuff to be done for every event
    resetCuts();
    
    // Define partons from the resonance
    TLorentzVector p1, p2;
    double p1_pdgId;
    double p2_pdgId;
    double parton1_pdgId;
    double parton2_pdgId;
    TLorentzVector parton1, parton2;
    TLorentzVector diparton;
    int n_partons = gen_pt->size();
    
    if(n_partons < 2) continue;
    
    p1.SetPxPyPzE(gen_px->at(n_partons -2) , gen_py->at(n_partons -2) , gen_pz->at(n_partons -2), gen_energy->at(n_partons -2) );
    p1_pdgId = gen_pdgId -> at(n_partons-2);
    p2.SetPxPyPzE(gen_px->at(n_partons -1 ) , gen_py->at(n_partons -1) , gen_pz->at(n_partons -1), gen_energy->at(n_partons -1) );
    p2_pdgId = gen_pdgId -> at(n_partons -1);
    
    // Re-order the partons in pt
    if( p1.Pt() > p2.Pt() ){
      parton1 = p1;
      parton2 = p2;
      parton1_pdgId = p1_pdgId;
      parton2_pdgId = p2_pdgId;
    }else{
      parton1 = p2;
      parton2 = p1;
      parton1_pdgId = p2_pdgId;
      parton2_pdgId = p1_pdgId;
    }
    
    if( parton1.Pt()<0 || parton2.Pt()<0 ) continue ;
    
    if(verbose){
      cout<<"Parton 1"<<endl;
      cout<< "Pt: "<<parton1.Pt()<<"         Eta: "<< parton1.Eta()<<"         Phi: "<<parton1.Phi()<<"          M: "<<parton1.M()<<" pdgId: "<<parton1_pdgId<<endl;
      cout<<"Parton 2"<<endl;
      cout<< "Pt: "<<parton2.Pt()<<"         Eta: "<< parton2.Eta()<<"         Phi: "<<parton2.Phi()<<"          M: "<<parton2.M()<<" pdgId: "<<parton2_pdgId<<endl;
    }
    // Create dijet(reco) system
    diparton = parton1 + parton2;
    
    //+++++++++++++++++++++++++++
    // Gen-level: construct Gen_widejet from genjet_ak4     
    size_t no_Genjets_ak4 = jetPtGenAK4->size();
    TLorentzVector Genjet1, Genjet2;
    TLorentzVector Genwj1_tmp, Genwj2_tmp;
    TLorentzVector Genwj1, Genwj2;
    TLorentzVector Genwdijet;
    double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");
    
    if(no_Genjets_ak4 < 2) continue;
    
    if( fabs(jetEtaGenAK4->at(0)) > getPreCutValue1("jetFidRegion") || jetPtGenAK4->at(0) < getPreCutValue1("pt0Cut") ) continue;
    if( fabs(jetEtaGenAK4->at(1)) > getPreCutValue1("jetFidRegion") || jetPtGenAK4->at(1) < getPreCutValue1("pt1Cut") ) continue;
    
    //TLorentzVector leading Genjet 
    Genjet1.SetPtEtaPhiM(jetPtGenAK4->at(0), jetEtaGenAK4->at(0), jetPhiGenAK4->at(0), jetMassGenAK4->at(0) );
    Genjet2.SetPtEtaPhiM(jetPtGenAK4->at(1), jetEtaGenAK4->at(1), jetPhiGenAK4->at(1), jetMassGenAK4->at(1) );
    
    for( Long64_t ijet=0;  ijet<no_Genjets_ak4;  ijet++){//genjet loop for ak4
      TLorentzVector currentGenJet;
      
      if(fabs(jetEtaGenAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
	 //	  && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
	 && jetPtGenAK4->at(ijet) > getPreCutValue1("ptCut")){   
	
	TLorentzVector currentGenJet;
	currentGenJet.SetPtEtaPhiM(jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet) );   
	
	double DeltaR1 = currentGenJet.DeltaR(Genjet1);
	double DeltaR2 = currentGenJet.DeltaR(Genjet2);
	
	if(DeltaR1 < DeltaR2 && DeltaR1 < GenwideJetDeltaR_){
	  Genwj1_tmp += currentGenJet;
	} else if(DeltaR2 < GenwideJetDeltaR_){
	  Genwj2_tmp += currentGenJet;
	}
      }    
    }//end of ak4 loop
    
    if( fabs(Genwj1_tmp.Eta() ) > 2.5 || fabs(Genwj2_tmp.Eta() ) > 2.5) continue;
    if( Genwj1_tmp.Pt() <30 || Genwj2_tmp.Pt() <30) continue;
    
    double DeltaR_GenWidejet1_tmp_parton1 = Genwj1_tmp.DeltaR(parton1);     
    double DeltaR_GenWidejet2_tmp_parton1 = Genwj2_tmp.DeltaR(parton1);     
    
    if(DeltaR_GenWidejet1_tmp_parton1 < DeltaR_GenWidejet2_tmp_parton1){
      Genwj1 = Genwj1_tmp;
      Genwj2 = Genwj2_tmp;
    } else{
      Genwj1 = Genwj2_tmp; 
      Genwj2 = Genwj1_tmp; 
    }
    
    double DeltaR_GenWideJet1_parton1;   
    double DeltaR_GenWideJet2_parton2;   
    DeltaR_GenWideJet1_parton1 = Genwj1.DeltaR(parton1);
    DeltaR_GenWideJet2_parton2 = Genwj2.DeltaR(parton2);
    
    CreateAndFillUserTH1D("H_DeltaR_GenWideJet1_Parton1", 100, 0. , 2. , DeltaR_GenWideJet1_parton1);
    CreateAndFillUserTH1D("H_DeltaR_GenWideJet2_Parton2", 100, 0. , 2. , DeltaR_GenWideJet2_parton2);     
    
    if(verbose){
      cout<<"GenWideJet 1"<<endl;
      cout<< "Pt: "<<Genwj1.Pt()<<"         Eta: "<< Genwj1.Eta()<<"         Phi: "<<Genwj1.Phi()<<"         M: "<<Genwj1.M()<<endl;
      cout<<"GenWideJet 2"<<endl;
      cout<< "Pt: "<<Genwj2.Pt()<<"         Eta: "<< Genwj2.Eta()<<"         Phi: "<<Genwj2.Phi()<<"         M: "<<Genwj2.M()<<endl;
    }
    
    // Create dijet system
    Genwdijet = Genwj1 + Genwj2;
    
    //+++++++++++++++++++++++++++
    // Reco-level: construct the WideJet from recojet ak4     
    size_t no_jets_ak4=jetPtAK4->size();
    TLorentzVector jet1, jet2, dijet;
    TLorentzVector wj1_tmp, wj2_tmp;
    TLorentzVector wj1, wj2, wdijet;
    double wideJetDeltaR_ = getPreCutValue1("DeltaR");
    
    if(no_jets_ak4<2) continue;
    
    if(fabs(jetEtaAK4->at(0)) > getPreCutValue1("jetFidRegion") || jetPtAK4->at(0) < getPreCutValue1("pt0Cut")) continue;
    if(fabs(jetEtaAK4->at(1)) > getPreCutValue1("jetFidRegion") || jetPtAK4->at(1) < getPreCutValue1("pt1Cut")) continue;
    
    //    TLorentzVector jet1, jet2 leading reco jet
    jet1.SetPtEtaPhiM(jetPtAK4->at(0),jetEtaAK4->at(0),jetPhiAK4->at(0),jetMassAK4->at(0));
    jet2.SetPtEtaPhiM(jetPtAK4->at(1),jetEtaAK4->at(1),jetPhiAK4->at(1),jetMassAK4->at(1));
    
    for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
      TLorentzVector currentJet;
      
      if(fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
	 && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
	 && jetPtAK4->at(ijet) > getPreCutValue1("ptCut")) {       
	
	TLorentzVector currentJet;  
	currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   
	
	double DeltaR1 = currentJet.DeltaR(jet1);
	double DeltaR2 = currentJet.DeltaR(jet2);
	
	if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
	  wj1_tmp += currentJet;
	} else if(DeltaR2 < wideJetDeltaR_){
	  wj2_tmp += currentJet;
	}
      }			 
    } //end of ak4 jet loop		     
    
    if( wj1_tmp.Pt()<0 || wj2_tmp.Pt()<0 ) continue;
    
    double DeltaR_Widejet1_tmp_parton1 = wj1_tmp.DeltaR(parton1);     
    double DeltaR_Widejet2_tmp_parton1 = wj2_tmp.DeltaR(parton1);     
    
    if(DeltaR_Widejet1_tmp_parton1 < DeltaR_Widejet2_tmp_parton1){
      wj1 = wj1_tmp;
      wj2 = wj2_tmp;
    } else{
      wj1 = wj2_tmp; 
      wj2 = wj1_tmp; 
    }
    
    // check the spatial resolution Widejet-Parton
    double Eta_difference_Widejet1_parton1 = wj1.Eta() - parton1.Eta();
    double Phi_difference_Widejet1_parton1 = wj1.Phi() - parton1.Phi() ;
    double DeltaR_WideJet1_parton1            = wj1.DeltaR(parton1);
    double Eta_difference_Widejet2_parton2 = wj2.Eta() - parton2.Eta() ;
    double Phi_difference_Widejet2_parton2 = wj2.Phi() - parton2.Phi() ;
    double DeltaR_WideJet2_parton2            = wj2.DeltaR(parton2);
    double Angle_Widejet1_parton1 = wj1.Angle(parton1.Vect());
    double Angle_Widejet2_parton2= wj2.Angle(parton2.Vect());
    double Angle_Partons= parton1.Angle(parton2.Vect());
    
    CreateAndFillUserTH1D("H_Eta_difference_Widejet1_parton1", 100, -2. , 2. , Eta_difference_Widejet1_parton1);
    CreateAndFillUserTH1D("H_Phi_difference_Widejet1_parton1", 100, -2. , 2. , Phi_difference_Widejet1_parton1);
    CreateAndFillUserTH1D("H_Eta_difference_Widejet2_parton2", 100, -2. , 2. , Eta_difference_Widejet2_parton2);
    CreateAndFillUserTH1D("H_Phi_difference_Widejet2_parton2", 100, -2. , 2. , Phi_difference_Widejet2_parton2);       
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_parton1", 100, 0. , 3. , DeltaR_WideJet1_parton1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_parton2", 100, 0. , 3. , DeltaR_WideJet1_parton1);
    CreateAndFillUserTH1D("H_Angle_difference_WideJet1_parton1", 100, -3.15 , 3.15, Angle_Widejet1_parton1 );
    CreateAndFillUserTH1D("H_Angle_difference_WideJet2_parton2", 100, -3.15 , 3.15, Angle_Widejet2_parton2 );
    CreateAndFillUserTH1D("H_Angle_difference_Partons", 100, -3.15 , 3.15, Angle_Partons );     
    
    if(verbose){
      cout<<"WideJet 1"<<endl;
      cout<< "Pt: "<<wj1.Pt()<<"         Eta: "<< wj1.Eta()<<"         Phi: "<<wj1.Phi()<<"         M: "<<wj1.M()<<endl;
      cout<<"WideJet 2"<<endl;
      cout<< "Pt: "<<wj2.Pt()<<"         Eta: "<< wj2.Eta()<<"         Phi: "<<wj2.Phi()<<"         M: "<<wj2.M()<<endl;
    }
    // Create dijet system
    wdijet = wj1 + wj2;
    
    //++++++++ SMEARING ++++++++
    
    double DeltaR_WideJet1_GenWideJet1 = wj1.DeltaR(Genwj1);     
    double DeltaR_WideJet2_GenWideJet2 = wj2.DeltaR(Genwj2);     
    
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_GenWideJet1", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_GenWideJet2", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
    
    double Parton_pdgId[2];
    double GenWideJet_Pt[2];
    double GenWideJet_Eta[2];
    double GenWideJet_Phi[2];
    double GenWideJet_M[2];    
    Parton_pdgId[0]      = parton1_pdgId;
    GenWideJet_Pt[0]   = Genwj1.Pt();
    GenWideJet_Eta[0] = Genwj1.Eta();
    GenWideJet_Phi[0] = Genwj1.Phi();
    GenWideJet_M[0]   = Genwj1.M();
    Parton_pdgId[1]      = parton2_pdgId;
    GenWideJet_Pt[1]   = Genwj2.Pt();
    GenWideJet_Eta[1] = Genwj2.Eta();
    GenWideJet_Phi[1] = Genwj2.Phi();
    GenWideJet_M[1]   = Genwj2.M();    
    
    // Geometrical matching
    if(DeltaR_WideJet1_GenWideJet1 > 0.3 || DeltaR_WideJet2_GenWideJet2 > 0.3) continue;
    
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_GenWideJet1_Match", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_GenWideJet2_Match", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
    
    double factor;
    double JetSmeared_Pt[2]   ={-1001, -1001};
    double JetSmeared_Eta[2] ={-1001, -1001};
    double JetSmeared_Phi[2] ={-1001, -1001};
    double JetSmeared_M[2]   ={-1001, -1001};
    
    TLorentzVector JetSmeared1_tmp, JetSmeared2_tmp;
    TLorentzVector JetSmeared1, JetSmeared2;
    
    for(int kk = 0 ; kk < 2 ; kk++){ // loop on GenWidejets
      
      int etaBin = mEtaBinning.getBin( fabs(GenWideJet_Eta[kk]) );
      int ptBin = mPtBinning.getPtBin( GenWideJet_Pt[kk] );
      
      const std::pair<float, float> pt_bin = mPtBinning.getBinValue(ptBin);
      std::stringstream ss;
      ss << "pt_" << (int) pt_bin.first << "_" << (int) pt_bin.second;      
      
      if(verbose){
	const std::pair<float, float> etaBins = mEtaBinning.getBinValue(etaBin);
	cout<<"etaGen  "<< fabs(GenWideJet_Eta[kk]) << "    pTGen   "<<  GenWideJet_Pt[kk] <<endl;
	cout<<"etaBin  "<< etaBin << "    pTBin   "<<ptBin<<endl;
	cout<<"etaBin.first  "<< etaBins.first << "    etaBin.second   "<<etaBins.second<<endl;
	cout<<"ptBin.first  "<< pt_bin.first << "    ptBin.second   "<<pt_bin.second<<endl;
      }
      
      TString responseName;
      
      // Categories
      if(n_categories ==2){
	if(Parton_pdgId[kk] == 21){
	  responseName = TString::Format("smearing_Gluon_%s_%s",  mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else{
	  responseName = TString::Format("smearing_Quark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	} 
      } else if(n_categories == 3 ){
	if(Parton_pdgId[kk] == 21){
	  responseName = TString::Format("smearing_Gluon_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 1 || fabs(Parton_pdgId[kk]) == 2  || fabs(Parton_pdgId[kk]) == 3 ){ 
	  responseName = TString::Format("smearing_LightQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 4 || fabs(Parton_pdgId[kk]) == 5 ){ 
	  responseName = TString::Format("smearing_HeavyQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}
      }else if(n_categories == 4 ){
	if(Parton_pdgId[kk] == 21){
	  responseName = TString::Format("smearing_Gluon_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 1 || fabs(Parton_pdgId[kk]) == 2  || fabs(Parton_pdgId[kk]) == 3 ){ 
	  responseName = TString::Format("smearing_LightQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 4 ){ 
	  responseName = TString::Format("smearing_CharmQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 5 ){ 
	  responseName = TString::Format("smearing_BottomQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}
      }else if(n_categories == 6 ){
	if(Parton_pdgId[kk] == 21){
	  responseName = TString::Format("smearing_Gluon_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 1){ 
	  responseName = TString::Format("smearing_UpQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 2){ 
	  responseName = TString::Format("smearing_DownQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 3){ 
	  responseName = TString::Format("smearing_StrangeQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 4 ){ 
	  responseName = TString::Format("smearing_CharmQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}else if( fabs(Parton_pdgId[kk]) == 5 ){ 
	  responseName = TString::Format("smearing_BottomQuark_%s_%s", mEtaBinning.getBinName(etaBin).c_str(), ss.str().c_str() );
	}
      }
      
      TH1F *h = (TH1F*)file1.Get("smearingFunction/"+responseName);
      
	      // if( h ->Integral() == 0){
	      // factor = 1;
	      // }else{
	      // factor = h ->GetRandom();
	      // }
      
      factor = h ->GetRandom();
      
      CreateAndFillUserTH1D("H_SmearingFactor", 150, 0 , 2. ,  factor);

      // Define smeared jet	      
      JetSmeared_Pt[kk]   = GenWideJet_Pt[kk] * factor ;
      JetSmeared_Eta[kk] = GenWideJet_Eta[kk] ;
      JetSmeared_Phi[kk] = GenWideJet_Phi[kk] ;
      JetSmeared_M[kk]   = GenWideJet_M[kk] ;

    }
  
    // check but probably it doesn't necessary    
    if(JetSmeared_Pt[0] <0 || JetSmeared_Pt[1] <0) continue;
    
    JetSmeared1_tmp.SetPtEtaPhiM(JetSmeared_Pt[0], JetSmeared_Eta[0], JetSmeared_Phi[0], JetSmeared_M[0] );
    JetSmeared2_tmp.SetPtEtaPhiM(JetSmeared_Pt[1], JetSmeared_Eta[1], JetSmeared_Phi[1], JetSmeared_M[1] );
    
    double dR1 = JetSmeared1_tmp.DeltaR(wj1);
    double dR2 = JetSmeared2_tmp.DeltaR(wj1);
    
    if(dR1 < dR2){
      JetSmeared1 = JetSmeared1_tmp ;
      JetSmeared2 = JetSmeared2_tmp ;
    }else{
      JetSmeared1 = JetSmeared2_tmp ;
      JetSmeared2 = JetSmeared1_tmp ;
    }
    
    TLorentzVector DijetSmeared;           
    DijetSmeared = JetSmeared1 + JetSmeared2 ;
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    //== Fill Variables ==
    
    fillVariableWithValue("run",runNo);     
    fillVariableWithValue("event",evtNo);     
    fillVariableWithValue("lumi",lumi);     
    fillVariableWithValue("nVtx",nvtx);     
    fillVariableWithValue("Parton1_pT", parton1.Pt() );     
    fillVariableWithValue("Parton1_Eta", parton1.Eta() );     
    fillVariableWithValue("Parton1_Phi", parton1.Phi() );     
    fillVariableWithValue("Parton1_M", parton1.M() ); 
    fillVariableWithValue("Parton1_pdgId", parton1_pdgId ); 
    fillVariableWithValue("Parton2_pT", parton2.Pt() );     
    fillVariableWithValue("Parton2_Eta", parton2.Eta() );     
    fillVariableWithValue("Parton2_Phi", parton2.Phi() );     
    fillVariableWithValue("Parton2_M", parton2.M() );     
    fillVariableWithValue("Parton2_pdgId", parton2_pdgId ); 
    fillVariableWithValue("Diparton_pT", diparton.Pt() );     
    fillVariableWithValue("Diparton_Eta", diparton.Eta() );     
    fillVariableWithValue("Diparton_Phi", diparton.Phi() );     
    fillVariableWithValue("Diparton_M", diparton.M() );     
    fillVariableWithValue("GenWideJet1_pT", Genwj1.Pt() );     
    fillVariableWithValue("GenWideJet1_Eta", Genwj1.Eta() );     
    fillVariableWithValue("GenWideJet1_Phi", Genwj1.Phi() );     
    fillVariableWithValue("GenWideJet1_M", Genwj1.M() );    
    fillVariableWithValue("GenWideJet2_pT", Genwj2.Pt() );     
    fillVariableWithValue("GenWideJet2_Eta", Genwj2.Eta() );     
    fillVariableWithValue("GenWideJet2_Phi", Genwj2.Phi() );     
    fillVariableWithValue("GenWideJet2_M", Genwj2.M() );
    fillVariableWithValue("GendijetWide_pT", Genwdijet.Pt() );     
    fillVariableWithValue("GendijetWide_Eta", Genwdijet.Eta() );     
    fillVariableWithValue("GendijetWide_Phi", Genwdijet.Phi() );     
    fillVariableWithValue("GendijetWide_M", Genwdijet.M() );     
    fillVariableWithValue("WideJet1_pT", wj1.Pt() );     
    fillVariableWithValue("WideJet1_Eta", wj1.Eta() );     
    fillVariableWithValue("WideJet1_Phi", wj1.Phi() );     
    fillVariableWithValue("WideJet1_M", wj1.M() );    
    fillVariableWithValue("WideJet2_pT", wj2.Pt() );     
    fillVariableWithValue("WideJet2_Eta", wj2.Eta() );     
    fillVariableWithValue("WideJet2_Phi", wj2.Phi() );     
    fillVariableWithValue("WideJet2_M", wj2.M() );
    fillVariableWithValue("dijetWide_pT", wdijet.Pt() );     
    fillVariableWithValue("dijetWide_Eta", wdijet.Eta() );     
    fillVariableWithValue("dijetWide_Phi", wdijet.Phi() );     
    fillVariableWithValue("dijetWide_M", wdijet.M() );     
    fillVariableWithValue("JetSmeared1_Pt", JetSmeared1.Pt() );
    fillVariableWithValue("JetSmeared1_Eta", JetSmeared1.Eta() );
    fillVariableWithValue("JetSmeared1_Phi", JetSmeared1.Phi() );
    fillVariableWithValue("JetSmeared1_M", JetSmeared1.M() );
    fillVariableWithValue("JetSmeared2_Pt", JetSmeared2.Pt() );
    fillVariableWithValue("JetSmeared2_Eta", JetSmeared2.Eta() );
    fillVariableWithValue("JetSmeared2_Phi", JetSmeared2.Phi() );
    fillVariableWithValue("JetSmeared2_M", JetSmeared2.M() );
    fillVariableWithValue("DijetSmeared_Pt", DijetSmeared.Pt() );
    fillVariableWithValue("DijetSmeared_Eta", DijetSmeared.Eta() );
    fillVariableWithValue("DijetSmeared_Phi", DijetSmeared.Phi() );
    fillVariableWithValue("DijetSmeared_M", DijetSmeared.M() );
    
    // Evaluate cuts (but do not apply them)
    evaluateCuts();
    
    //see
     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     /*  
     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedAllPreviousCuts("mjj") && passedCut("mjj") ) 
     {
	 	 fillReducedSkimTree();
		 
		 // ===== Take a look at this =====
		 // //Example on how to investigate quickly the data
 	 // if(getVariableValue("mjj")>4000)
	 //   {
	 //     //fast creation and filling of histograms
	 //     CreateAndFillUserTH1D("h_dphijj_mjjgt4000", 100, 0, 3.15, getVariableValue("deltaPHIjj"));
	 //     CreateAndFillUserTH1D("h_htak4_mjjgt4000", 1000, 0, 10000, getVariableValue("HTAK4"));
	 //     CreateAndFillUserTH1D("h_nvtx_mjjgt4000", 31, -0.5, 30.5, getVariableValue("nVtx"));
	 //   }

            } //fine see
*/

     // ===== Example of mjj spectrum after HLT selection =====
     // if( passedAllPreviousCuts("mjj") )
     //   {
     // 	 if(getVariableValue("passHLT")>0)
     // 	   {
     // 	     //fast creation and filling of histograms
     // 	     CreateAndFillUserTH1D("h_mjj_passHLT", getHistoNBins("mjj"), getHistoMin("mjj"), getHistoMax("mjj"), getVariableValue("mjj"));
     // 	   }
     //   }

     // reject events that did not pass level 0 cuts
     //if( !passedCut("0") ) continue;
     // ......
     
     // reject events that did not pass level 1 cuts
     //if( !passedCut("1") ) continue;
     // ......

     // reject events that did not pass the full cut list
     //if( !passedCut("all") ) continue;
     // ......

     // if( widejets.size() >= 2) {
     //  h_nJetFinal->Fill(widejets.size());
     //  h_DijetMass->Fill(wdijet.M());
     //  h_pT1stJet->Fill(widejets[0].Pt());
     //  h_pT2ndJet->Fill(widejets[1].Pt());
     //  h_eta1stJet->Fill(widejets[0].Eta());
     //  h_eta2ndJet->Fill(widejets[1].Eta());
     // }
     ////////////////////// User's code ends here ///////////////////////

            
   } // End loop over events
   

   //////////write histos 

   // h_nVtx->Write();
   // h_trueVtx->Write();
   // h_nJetFinal->Write();
   // h_pT1stJet->Write();
   // h_pT2ndJet->Write();
   // h_DijetMass->Write();
   // h_eta1stJet->Write();
   // h_eta2ndJet->Write();

   // //pT of both jets, to be built using the histograms produced automatically by baseClass
   // TH1F * h_pTJets = new TH1F ("h_pTJets","", getHistoNBins("pT1stJet"), getHistoMin("pT1stJet"), getHistoMax("pT1stJet"));
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT1stJet") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT2ndJet") );
   // //one could also do:  *h_pTJets = getHisto_noCuts_or_skim("pT1stJet") + getHisto_noCuts_or_skim("pT2ndJet");
   // h_pTJets->Write();
   // //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

