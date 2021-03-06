#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TMath.h>

//cambia file per leggere le giuste funzioni di smearing

 TFile file1("/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/FourCategory/SmearingFunctions_QCDSamples_Bin200.root");


TH1D *H_step_pt = (TH1D*)file1.Get("H_step_pt");
int step_pt = H_step_pt->GetMean();
int n_bin;

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

   char HistoName[200];	 
   n_bin = 4500 / step_pt ;   

   /*
     TH1D *h[100][100];
     for(int ii=0; ii<5; ii++){
     for(int jj=0; jj<n_bin ; jj++){	 	 
     int pt_bin_max = (jj*step_pt)+step_pt;
     sprintf(HistoName,"Histo_R_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
     h[ii][jj] = new TH1D(HistoName, HistoName, 50,0.,2.);   
     sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
     h[ii][jj] = new TH1D(HistoName, HistoName, 50,0.,2.);   
      }
     }
   */
   cout <<"Binning in pT = "<< step_pt <<endl;
   
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

     // Uso i dispari per le funzioni di smearing, quindi qui uso i pari
     // Se uso le funzioni di smearing calcolate su samples di QCD posso usare tutti gli eventi
     // if( (int)jentry%2 ==! 0 ) continue;
     
     // cout << " " << endl;
     // cout <<"Evento numero "<< jentry << endl;

     
     // Definisco i partoni prodotti dalla risonanza     

     TLorentzVector p1, p2;
     double p1_pdgId;
     double p2_pdgId;
     double parton1_pdgId;
     double parton2_pdgId;
     TLorentzVector parton1, parton2;
     TLorentzVector diparton;

     if(gen_pt->size() == 3){
       p1.SetPxPyPzE(gen_px->at(1) , gen_py->at(1) , gen_pz->at(1), gen_energy->at(1) );
       p1_pdgId     = gen_pdgId->at(1);
       p2.SetPxPyPzE(gen_px->at(2) , gen_py->at(2) , gen_pz->at(2), gen_energy->at(2) );
       p2_pdgId     = gen_pdgId->at(2);
     }

     if(gen_pt->size() == 4){
       p1.SetPxPyPzE(gen_px->at(2) , gen_py->at(2) , gen_pz->at(2), gen_energy->at(2) );
       p1_pdgId     = gen_pdgId->at(2);
       p2.SetPxPyPzE(gen_px->at(3) , gen_py->at(3) , gen_pz->at(3), gen_energy->at(3) );
       p2_pdgId     = gen_pdgId->at(3);
     }

     if(gen_pt->size() == 5){
       p1.SetPxPyPzE(gen_px->at(3) , gen_py->at(3) , gen_pz->at(3), gen_energy->at(3) );
       p1_pdgId     = gen_pdgId->at(3);
       p2.SetPxPyPzE(gen_px->at(4) , gen_py->at(4) , gen_pz->at(4), gen_energy->at(4) );
       p2_pdgId     = gen_pdgId->at(4);
     }

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

     // cout<<"Partone 1"<<endl;
     // cout<< "Pt: "<<parton1.Pt()<<"         Eta: "<< parton1.Eta()<<"         Phi: "<<parton1.Phi()<<"          M: "<<parton1.M()<<" pdgId: "<<parton1_pdgId<<endl;
     // cout<<"Partone 2"<<endl;
     // cout<< "Pt: "<<parton2.Pt()<<"         Eta: "<< parton2.Eta()<<"         Phi: "<<parton2.Phi()<<"          M: "<<parton2.M()<<" pdgId: "<<parton2_pdgId<<endl;
 
     // Create dijet(reco) system
     diparton = parton1 + parton2;
     ///////////////////////////////////////////////////////////////////////////////////////
     //Livello gen
     //ricostruisce i Gen_widejet partendo dai genjet_ak4
     
     size_t no_Genjets_ak4 = jetPtGenAK4->size();
     TLorentzVector Genjet1, Genjet2;
     TLorentzVector Genwj1_tmp, Genwj2_tmp;
     TLorentzVector Genwj1, Genwj2;
     TLorentzVector Genwdijet;
     double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     // if(no_Genjets_ak4 < 2) cout<<"# Genjet AK4 <2. Butto evento"<<endl;
     if(no_Genjets_ak4 < 2) continue;
     
     // cout<<"Costruisco GenWidejet"<<endl;
     // cout<<"Ciclo sui recojet AK4"<<endl;
     
     //TLorentzVector Genjet1, Genjet2; 
     Genjet1.SetPtEtaPhiM(jetPtGenAK4->at(0), jetEtaGenAK4->at(0), jetPhiGenAK4->at(0), jetMassGenAK4->at(0) );
     Genjet2.SetPtEtaPhiM(jetPtGenAK4->at(1), jetEtaGenAK4->at(1), jetPhiGenAK4->at(1), jetMassGenAK4->at(1) );
     
     // costruisco Gen widejet attorno ai genjet a pT piu alto
     for( Long64_t ijet=0;  ijet<no_Genjets_ak4;  ijet++){//genjet loop for ak4
       
       TLorentzVector currentGenJet;
       currentGenJet.SetPtEtaPhiM(jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet) );   
       
       //cout<<"Genjet numero "<<ijet<<endl;
       //cout<<"Pt: "<<currentGenJet.Pt()<<"          Eta: "<<currentGenJet.Eta()<<"         Phi: "<<currentGenJet.Phi()<<"         M: "<<currentGenJet.M()<<endl;
       
       double DeltaR1 = currentGenJet.DeltaR(Genjet1);
       double DeltaR2 = currentGenJet.DeltaR(Genjet2);
       //cout<<"DeltaR tra il genjet numero "<<ijet<<" e il genjet 0: "<<DeltaR1<<endl;
       //cout<<"DeltaR tra il genjet numero "<<ijet<<" e il genjet 1: "<<DeltaR2<<endl;
       
       if(DeltaR1 < DeltaR2 && DeltaR1 < GenwideJetDeltaR_){
	 //cout<<"Sommo il jet numero "<<ijet<<" al genjet numero 0"<<endl;
	 Genwj1_tmp += currentGenJet;
       }
       else if(DeltaR2 < GenwideJetDeltaR_){
	 //cout<<"Sommo il jet numero "<<ijet<<" al genjet numero 1"<<endl;
	 Genwj2_tmp += currentGenJet;
       }
     }//end of ak4 loop
     
     // if( Genwj1_tmp.Pt() <0 || Genwj2_tmp.Pt() <0) cout<<"Pt GenWidejet <0"<<endl;
     if( Genwj1_tmp.Pt() <0 || Genwj2_tmp.Pt() <0) continue;
     
     //associo i Gen widejets ai partoni
     double DeltaR_GenWidejet1_tmp_parton1 = Genwj1_tmp.DeltaR(parton1);     
     double DeltaR_GenWidejet2_tmp_parton1 = Genwj2_tmp.DeltaR(parton1);     
     
     if(DeltaR_GenWidejet1_tmp_parton1 < DeltaR_GenWidejet2_tmp_parton1){
       //cout<<"GenWidejet_tmp 1 associato al partone 1 Allora Genwj1 = Genwj1_tmp"<<endl;
       //cout<<"GenWidejet_tmp 2 associato al partone 2 Allora Genwj2 = Genwj2_tmp"<<endl;
       Genwj1 = Genwj1_tmp;
       Genwj2 = Genwj2_tmp;
     }
     else{
       //cout<<"GenWidejet_tmp 2 associato al partone 1 Allora Genwj1 = Genwj2_tmp"<<endl;
       //cout<<"GenWidejet_tmp 1 associato al partone 2 Allora Genwj2 = Genwj1_tmp"<<endl;
       Genwj1 = Genwj2_tmp; 
       Genwj2 = Genwj1_tmp; 
     }
     
     // cout<<"GenWideJet 1"<<endl;
     // cout<< "Pt: "<<Genwj1.Pt()<<"         Eta: "<< Genwj1.Eta()<<"         Phi: "<<Genwj1.Phi()<<"         M: "<<Genwj1.M()<<endl;
     // cout<<"GenWideJet 2"<<endl;
     // cout<< "Pt: "<<Genwj2.Pt()<<"         Eta: "<< Genwj2.Eta()<<"         Phi: "<<Genwj2.Phi()<<"         M: "<<Genwj2.M()<<endl;
     
     //funzioni di smearing le calcolo fino a questo range   
     // if( fabs(Genwj1.Eta() ) > 2.5 || fabs(Genwj2.Eta() ) > 2.5) cout<<"Uno dei due GenWideJet e' fuori range in eta"<<endl;
     if( fabs(Genwj1.Eta() ) > 2.5 || fabs(Genwj2.Eta() ) > 2.5) continue; //see
     
     // Create dijet system
     Genwdijet = Genwj1 + Genwj2;
    
      ////////////////////////////////////////////////////////////////
     //Livello reco -> AK4 dentro la ntupla
     //ricostruisco i widejet partendo di recojet_ak4
     
     size_t no_jets_ak4=jetPtAK4->size();
     TLorentzVector jet1, jet2, dijet;
     TLorentzVector wj1_tmp, wj2_tmp;
     TLorentzVector wj1, wj2, wdijet;
     double wideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     // if(no_jets_ak4<2) cout<<"# recojet AK4 < 2. Butto evento"<<endl;;    
     if(no_jets_ak4<2) continue;
     
     // cout<<"Costruisco Widejet"<<endl;    
     // cout<<"Ciclo sui recojet AK4"<<endl;
     
     //    TLorentzVector jet1, jet2;
     jet1.SetPtEtaPhiM(jetPtAK4->at(0),jetEtaAK4->at(0),jetPhiAK4->at(0),jetMassAK4->at(0));
     jet2.SetPtEtaPhiM(jetPtAK4->at(1),jetEtaAK4->at(1),jetPhiAK4->at(1),jetMassAK4->at(1));
     
     //costruisco reco widejet attorno ai jet a pT piu alto
     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
       
       TLorentzVector currentJet;
       currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   
       
       // cout<<"Jet numero "<<ijet<<endl;
       // cout<<"Pt: "<<currentJet.Pt()<<"         Eta: "<<currentJet.Eta()<<"          Phi: "<<currentJet.Phi()<<"         M: "<<currentJet.M()<<endl;
       
       double DeltaR1 = currentJet.DeltaR(jet1);
       double DeltaR2 = currentJet.DeltaR(jet2);
       //cout<<"DeltaR tra il jet numero "<<ijet<<" e il jet 0 : "<<DeltaR1<<endl;
       //cout<<"DeltaR tra il jet numero "<<ijet<<" e il jet 1 : "<<DeltaR2<<endl; 
       
       if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
	 //cout<<"Sommo il jet numero "<<ijet<<" al jet numero 0"<<endl;
	 wj1_tmp += currentJet;
      }
       else if(DeltaR2 < wideJetDeltaR_){
	 //cout<<"Sommo il jet numero "<<ijet<<" al jet numero 1"<<endl;
	 wj2_tmp += currentJet;
       }			 
    } //end of ak4 jet loop		     
     
     // if( wj1_tmp.Pt()<0 || wj2_tmp.Pt()<0 ) cout<<"Pt Widejet <0"<<endl;;
     if( wj1_tmp.Pt()<0 || wj2_tmp.Pt()<0 ) continue;
     
     //associo i widejets ai partoni
     double DeltaR_Widejet1_tmp_parton1 = wj1_tmp.DeltaR(parton1);     
     double DeltaR_Widejet2_tmp_parton1 = wj2_tmp.DeltaR(parton1);     
     
     if(DeltaR_Widejet1_tmp_parton1 < DeltaR_Widejet2_tmp_parton1){
       //cout<<"Widejet_tmp 1 associato al partone 1 Allora wj1 = wj1_tmp"<<endl;
       //cout<<"Widejet_tmp 2 associato al partone 2 Allora wj2 = wj2_tmp"<<endl;
       wj1 = wj1_tmp;
       wj2 = wj2_tmp;
     }
     else{
       //cout<<"Widejet_tmp 2 associato al partone 1 Allora wj1 = wj2_tmp"<<endl;
       //cout<<"Widejet_tmp 1 associato al partone 2 Allora wj2 = wj1_tmp"<<endl;
       wj1 = wj2_tmp; 
       wj2 = wj1_tmp; 
     }
     
     // cout<<"WideJet 1"<<endl;
     // cout<< "Pt: "<<wj1.Pt()<<"         Eta: "<< wj1.Eta()<<"         Phi: "<<wj1.Phi()<<"         M: "<<wj1.M()<<endl;
     // cout<<"WideJet 2"<<endl;
     // cout<< "Pt: "<<wj2.Pt()<<"         Eta: "<< wj2.Eta()<<"         Phi: "<<wj2.Phi()<<"         M: "<<wj2.M()<<endl;
     
     // Create dijet system
    wdijet = wj1 + wj2;
    
    //+++++++++++++++++++++ FUNZIONI DI SMEARING ++++++++++++++++++++++++++++++++++++++
    
    double DeltaR_WideJet1_GenWideJet1 = wj1.DeltaR(Genwj1);     
    double DeltaR_WideJet2_GenWideJet2 = wj2.DeltaR(Genwj2);     
    
    // cout<< "DeltaR tra widejet 1 e Genwidejet 1 = "<<DeltaR_WideJet1_GenWideJet1<<endl;
    // cout<< "DeltaR tra widejet 2 e Genwidejet 2 = "<<DeltaR_WideJet2_GenWideJet2<<endl;
    
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

    // if(DeltaR_WideJet1_GenWideJet1 > 0.3 || DeltaR_WideJet2_GenWideJet2 > 0.3) cout<<"Uno dei due widejet non matchea il genwidejet . Butto evento"<<endl;
    if(DeltaR_WideJet1_GenWideJet1 > 0.3 || DeltaR_WideJet2_GenWideJet2 > 0.3) continue;
    
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_GenWideJet1_Match", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_GenWideJet2_Match", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
    
    // se commentato -> sto usando i jet in ogni bin     
    // see // considero i GenWideJets in un solo bin di |Eta| e pT
    // if( fabs(GenWideJet_Eta[0] ) >=0.5 || fabs( GenWideJet_Eta[1] ) >= 0.5) continue;
    // if( GenWideJet_Pt[0] <1400 || GenWideJet_Pt[0] >=1600) continue;
    // if( GenWideJet_Pt[1] <1400 || GenWideJet_Pt[1] >=1600) continue;
    
    // cout<<"CALCOLO I MIEI JET SMEARED"<<endl;
    
    double factor;
    double JetSmeared_Pt[2]   ={-1001, -1001};
    double JetSmeared_Eta[2] ={-1001, -1001};
    double JetSmeared_Phi[2] ={-1001, -1001};
    double JetSmeared_M[2]   ={-1001, -1001};
       
    TLorentzVector JetSmeared1_tmp, JetSmeared2_tmp;
    TLorentzVector JetSmeared1, JetSmeared2;
    
    // cout<<" " <<endl;
    // cout<<"Costruisco i miei jet smeared"<<endl; 
    
    for(int kk = 0 ; kk < 2 ; kk++){
      
      if(kk == 0) CreateAndFillUserTH1D("H_DeltaR_Totale", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
      if(kk == 1) CreateAndFillUserTH1D("H_DeltaR_Totale", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
      
      // cout<<"GenWideJet "<< kk+1<<endl;
      // cout<<"Pt: "<<GenWideJet_Pt[kk]<<"         |Eta|: "<< fabs(GenWideJet_Eta[kk])<<endl;
      
      for(int ii=0; ii<5; ii++){
	double eta_bin_min = ii/2.;       
	double eta_bin_max = ii/2. +0.5;   
	
	if( fabs(GenWideJet_Eta[kk])>=eta_bin_min && fabs(GenWideJet_Eta[kk]) < eta_bin_max){
	  // cout<<"Passed bin in eta"<<endl;	
	  // cout << "eta_bin: "<< eta_bin_min <<" - "<<eta_bin_max << endl;
	  
	  for(int jj=0; jj<n_bin ; jj++){	 
	    int pt_bin_min = ( (jj-1)*step_pt)+step_pt;
	    int pt_bin_max = (jj*step_pt)+step_pt;
	    
	    if(GenWideJet_Pt[kk] >=pt_bin_min && GenWideJet_Pt[kk] < pt_bin_max){ 
	      // cout<<"Passed bin in pt"<<endl;
	      // cout << "pt_bin: "<< pt_bin_min <<" - "<<pt_bin_max<< endl;
	      
	      if( Parton_pdgId[kk] == 21){
		sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	      }
	      if( fabs(Parton_pdgId[kk]) == 1 || fabs(Parton_pdgId[kk]) == 2 || fabs(Parton_pdgId[kk]) == 3 ){ 
		sprintf(HistoName,"Histo_RQuarkLIGHT_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	      }
	      if( fabs(Parton_pdgId[kk]) == 4){ 
		sprintf(HistoName,"Histo_RQuarkCHARM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	      }
	      if( fabs(Parton_pdgId[kk]) == 5){ 
		sprintf(HistoName,"Histo_RQuarkBOTTOM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	      }
	      
	      // cout<<"Estraggo il fattore dall'istogramma: "<< HistoName << endl;	      
	      
	      TH1D *h = (TH1D*)file1.Get(HistoName);		
	      //cout<<"h[ii][jj]-> Integral(): "<<h[ii][jj]->Integral()<<endl;
	      
	      if( h ->Integral() == 0){
		// cout<<HistoName<<" ha integrale nullo"<<endl;
		factor = 1;
	      }else{
		factor = h ->GetRandom();
	      }
	      
	      // cout<<"Fattore estratto: "<<factor<<endl;
	      
	      JetSmeared_Pt[kk]   = GenWideJet_Pt[kk] * factor ;
	      JetSmeared_Eta[kk] = GenWideJet_Eta[kk] ;
	      JetSmeared_Phi[kk] = GenWideJet_Phi[kk] ;
	      JetSmeared_M[kk]   = GenWideJet_M[kk] ;
	      
	      // cout<<"Jet Smeared "<<kk+1<<endl;
	      // cout<<"Pt: "<<JetSmeared_Pt[kk]<<"          Eta: "<<JetSmeared_Eta[kk]<<"          Phi: "<<JetSmeared_Phi[kk]<<"          M: "<<JetSmeared_M[kk]<<endl;
	      
	    }//if range di pt
	  }//for sui bin in pt
	}//if range di eta
      }//for sui bin in eta
    }//for sui due GenWideJets
    
    //se non trova il match -> quindi non ci sono due jet smeared    
    // if(JetSmeared_Pt[0] <0 || JetSmeared_Pt[1] <0) cout<<"Non ci sono due jet smeared. Butto evento"<<endl;;
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
    
    // cout<<"Massa invariante dei partoni: "<<diparton.M()<<endl;
    // cout<<"Massa invariante dei widejets: "<<wdijet.M()<<endl;
    // cout<<"Massa invariante dei Genwidejets: "<<Genwdijet.M()<<endl;
    // cout<<"Massa invariante dei jets smeared: "<<DijetSmeared.M()<<endl;
    ///////////////////////////////////////////////////////////////////     
    // Riempio qui, cosi hanno gli stessi eventi
    //Partoni
    CreateAndFillUserTH1D("H_parton1_Pt", 1000, 0. , 10000. , parton1.Pt() );
    CreateAndFillUserTH1D("H_parton1_Eta", 100, -4. , 4. , parton1.Eta() );
    CreateAndFillUserTH1D("H_parton1_Phi", 100, -4. , 4. , parton1.Phi() );
    CreateAndFillUserTH1D("H_parton1_M", 50000, 0. , 10000. , parton1.M() );
    CreateAndFillUserTH1D("H_parton1_pdgId", 80, -10. , 10. , parton1_pdgId );
    CreateAndFillUserTH1D("H_parton2_Pt", 1000, 0. , 10000. , parton2.Pt() );
    CreateAndFillUserTH1D("H_parton2_Eta", 100, -4. , 4. , parton2.Eta() );  
    CreateAndFillUserTH1D("H_parton2_Phi", 100, -4. , 4. , parton2.Phi() );
    CreateAndFillUserTH1D("H_parton2_M", 50000, 0. , 10000. , parton2.M() );
    CreateAndFillUserTH1D("H_parton2_pdgId", 80, -10. , 10. , parton2_pdgId );
    //Variabili del sistema diparton
    CreateAndFillUserTH1D("H_diparton_Pt", 1000, 0. , 10000. , diparton.Pt() );
    CreateAndFillUserTH1D("H_diparton_Eta", 100, -4. , 4. , diparton.Eta() );
    CreateAndFillUserTH1D("H_diparton_Phi", 100, -4. , 4. , diparton.Phi() );
    CreateAndFillUserTH1D("H_diparton_M", 1000, 0. , 10000. , diparton.M() );
    //Gen WideJet
    CreateAndFillUserTH1D("H_GenWideJet1_Pt", 1000, 0. , 10000. , Genwj1.Pt()); 
    CreateAndFillUserTH1D("H_GenWideJet1_Eta", 100, -4. , 4. , Genwj1.Eta()); 
    CreateAndFillUserTH1D("H_GenWideJet1_Phi", 100, -4. , 4. , Genwj1.Phi());
    CreateAndFillUserTH1D("H_GenWideJet1_M", 1000, 0. , 10000. , Genwj1.M());  
    CreateAndFillUserTH1D("H_GenWideJet2_Pt", 1000, 0. , 10000. , Genwj2.Pt()); 
    CreateAndFillUserTH1D("H_GenWideJet2_Eta", 100, -4. , 4. , Genwj2.Eta()); 
    CreateAndFillUserTH1D("H_GenWideJet2_Phi", 100, -4. , 4. , Genwj2.Phi());
    CreateAndFillUserTH1D("H_GenWideJet2_M", 1000, 0. , 10000. , Genwj2.M());  
    //Variabili del sistema Gen Dijet Wide	   
    CreateAndFillUserTH1D("H_GendijetWide_Pt", 1000, 0. , 10000. , Genwdijet.Pt() );
    CreateAndFillUserTH1D("H_GendijetWide_Eta", 100, -4. , 4. , Genwdijet.Eta() );
    CreateAndFillUserTH1D("H_GendijetWide_Phi", 100, -4. , 4. , Genwdijet.Phi() );
    CreateAndFillUserTH1D("H_GendijetWide_M", 1000, 0. , 10000. , Genwdijet.M() );
    // Reco WideJet
    CreateAndFillUserTH1D("H_WideJet1_Pt", 1000, 0. , 10000. , wj1.Pt()); 
    CreateAndFillUserTH1D("H_WideJet1_Eta", 100, -4. , 4. , wj1.Eta()); 
    CreateAndFillUserTH1D("H_WideJet1_Phi", 100, -4. , 4. , wj1.Phi());
    CreateAndFillUserTH1D("H_WideJet1_M", 1000, 0. , 10000. , wj1.M());  
    CreateAndFillUserTH1D("H_WideJet2_Pt", 1000, 0. , 10000. , wj2.Pt()); 
    CreateAndFillUserTH1D("H_WideJet2_Eta", 100, -4. , 4. , wj2.Eta()); 
    CreateAndFillUserTH1D("H_WideJet2_Phi", 100, -4. , 4. , wj2.Phi());
    CreateAndFillUserTH1D("H_WideJet2_M", 1000, 0. , 10000. , wj2.M());  
    CreateAndFillUserTH1D("H_dijetWide_Pt", 1000, 0. , 10000. , wdijet.Pt() );
    CreateAndFillUserTH1D("H_dijetWide_Eta", 100, -4. , 4. , wdijet.Eta() );
    CreateAndFillUserTH1D("H_dijetWide_Phi", 100, -4. , 4. , wdijet.Phi() );
    CreateAndFillUserTH1D("H_dijetWide_M", 1000, 0. , 10000. , wdijet.M() );       
    // Jet Smeared
    CreateAndFillUserTH1D("H_JetSmeared1_Pt", 1000, 0. , 10000. ,JetSmeared1.Pt()); 
    CreateAndFillUserTH1D("H_JetSmeared1_Eta", 100, -4. , 4. ,JetSmeared1.Eta()); 
    CreateAndFillUserTH1D("H_JetSmeared1_Phi", 100, -4. , 4. ,JetSmeared1.Phi()); 
    CreateAndFillUserTH1D("H_JetSmeared1_M", 1000, 0. , 10000. ,JetSmeared1.M()); 
    CreateAndFillUserTH1D("H_JetSmeared2_Pt", 1000, 0. , 10000. ,JetSmeared2.Pt()); 
    CreateAndFillUserTH1D("H_JetSmeared2_Eta", 100, -4. , 4. ,JetSmeared2.Eta()); 
    CreateAndFillUserTH1D("H_JetSmeared2_Phi", 100, -4. , 4. ,JetSmeared2.Phi()); 
    CreateAndFillUserTH1D("H_JetSmeared2_M", 1000, 0. , 10000. ,JetSmeared2.M()); 
    CreateAndFillUserTH1D("H_DijetSmeared_Pt", 1000, 0. , 10000. , DijetSmeared.Pt()); 
    CreateAndFillUserTH1D("H_DijetSmeared_Eta", 100, -4. , 4. , DijetSmeared.Eta()); 
    CreateAndFillUserTH1D("H_DijetSmeared_Phi", 100, -4. , 4. , DijetSmeared.Phi()); 
    CreateAndFillUserTH1D("H_DijetSmeared_M", 1000, 0. , 10000. , DijetSmeared.M());   
    
    //////////////////////////////////////////////////////////////////////////////////////     
    //== Fill Variables ==
    
    //fillVariableWithValue("run",runNo);     
    //fillVariableWithValue("event",evtNo);     
    //fillVariableWithValue("lumi",lumi);     
    //fillVariableWithValue("nVtx",nvtx);     
    //fillVariableWithValue("nJet",widejets.size());
    
    // Trigger
    int NtriggerBits = triggerResult->size();
    if( NtriggerBits > 0)
      //fillVariableWithValue("passHLT",triggerResult->at(0));// HLT_PFHT900_v*    

      //     if( no_jets_ak4 >=1 )
       //fillVariableWithValue("IdTight_j1",idTAK4->at(0));

       //     if( no_jets_ak4 >=2 )
       //fillVariableWithValue("IdTight_j2",idTAK4->at(1));

       //     if( widejets.size() >= 1 )
       //       {
         //fillVariableWithValue( "pT_j1", widejets[0].Pt() );
         //fillVariableWithValue( "eta_j1", widejets[0].Eta());

	 //no cuts on these variables, just to store in output
         //fillVariableWithValue( "phi_j1", widejets[0].Phi());
         //fillVariableWithValue( "neutrHadEnFrac_j1", jetNhfAK4->at(0));
         //fillVariableWithValue( "chargedHadEnFrac_j1", jetChfAK4->at(0));
         //fillVariableWithValue( "photonEnFrac_j1", jetPhfAK4->at(0));
         //fillVariableWithValue( "eleEnFract_j1", jetElfAK4->at(0));
         //fillVariableWithValue( "muEnFract_j1", jetMufAK4->at(0));
       //}

     //     if( widejets.size() >= 2 )
     // {
         //fillVariableWithValue( "pT_j2", widejets[1].Pt() );
         //fillVariableWithValue( "eta_j2", widejets[1].Eta());
	 //fillVariableWithValue( "GendeltaETAjj_hp", diparton_DeltaEtaJJ ) ;
         //fillVariableWithValue( "Genmjj_hp", diparton_MJJ ) ;
	 //fillVariableWithValue( "deltaETAjj", DeltaEtaJJ ) ;
         //fillVariableWithValue( "mjj", MJJ ) ; 
	 //fillVariableWithValue( "deltaETAjj_wide", DeltaEtaJJWide ) ;
         //fillVariableWithValue( "mjj_wide", MJJWide ) ;


	 //no cuts on these variables, just to store in output
         //fillVariableWithValue( "phi_j2", widejets[1].Phi());	
         //fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(1));
         //fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(1));
         //fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(1));
         //fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(1));
         //fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(1));
	 //fillVariableWithValue( "deltaPHIjj", DeltaPhiJJWide ) ;

     //       }

     //no cuts on these variables, just to store in output
     //fillVariableWithValue("trueVtx",PileupInteractions->at(12));
     //fillVariableWithValue("MET",met);
     double METoverHTAK4=double(met/htAK4);
     //fillVariableWithValue("METoverHTAK4",METoverHTAK4);
     //fillVariableWithValue("HTAK4",htAK4);
     //fillVariableWithValue("ptHat",ptHat);

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



////////////////////////////////////////////////////////////////////////

// rapporto con i genjet AK4/AK8 non servono piu perche abbiamo visto che hanno una risoluzione peggiore rispetto ai widejet
// sia generatore che reco non servono

    /*    
    double deltaR_currentGenjet_parton1 = 1000;
    double deltaRmin_Genjet_parton1 = 1000;
    int Id1_min = -1000;
    
    // cout<< " " << endl;
    // cout<< "COMINCIA QUI"<<endl;
    
    for(Long64_t ijet=0; ijet<no_Genjets_ak4; ijet++){ 
      
      TLorentzVector currentGenJet;
      currentGenJet.SetPtEtaPhiM( jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet));   
      
      deltaR_currentGenjet_parton1 = currentGenJet.DeltaR(parton1);
      
      //	   cout<< "deltaR_currentGenjet_parton1: "<< deltaR_currentGenjet_parton1 << endl;
      
      if( deltaR_currentGenjet_parton1 < deltaRmin_Genjet_parton1){
	
	     deltaRmin_Genjet_parton1 = deltaR_currentGenjet_parton1;
	     Id1_min = ijet;
	     //cout<< "DeltaR1 minimo: "<< deltaRmin_Genjet_parton1<<endl;
	     //cout<< "ID1_minimo: "<< Id1_min<<endl;
      }
    }
    
    //cout<< "DeltaR1 minimo: "<< deltaRmin_Genjet_parton1<<endl;
    //cout<< "ID1_minimo: "<< Id1_min<<endl;
      
    if(Id1_min < -1) continue;;
    
    CreateAndFillUserTH1D("H_DeltaR_Gen_hp1", 100, 0. , 10. , deltaRmin_Genjet_parton1 );
    CreateAndFillUserTH1D("H_Id1_min_Gen", 16, -0.5 , 15.5 , Id1_min );
    
    if(deltaRmin_Genjet_parton1 <= 0.3){
      //Calcolo ratio con un taglio sensato in DeltaR
      
      double     R_PtGenjet_Ptjet1hp;
      R_PtGenjet_Ptjet1hp = jetPtGenAK4->at(Id1_min) / parton1.Pt() ;  	
      
      CreateAndFillUserTH1D("H_R_PtGenjet_Ptjet1hp", 100, 0. , 10. , R_PtGenjet_Ptjet1hp );
      CreateAndFillUserTProfile("Profile_PtGen", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_PtGenjet_Ptjet1hp );  
      
      if(	fabs( parton1.Eta() ) >= 0. && fabs( parton1.Eta() ) < 0.5){ 
	CreateAndFillUserTProfile("Profile_PtGen_05", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_PtGenjet_Ptjet1hp ); 	 
      }
      if(	fabs( parton1.Eta() ) >= 0.5 && fabs( parton1.Eta() ) < 1.){ 
	CreateAndFillUserTProfile("Profile_PtGen_1", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_PtGenjet_Ptjet1hp ); 
      }
      if(	fabs( parton1.Eta() ) >= 1. && fabs( parton1.Eta() ) < 1.5){ 
	CreateAndFillUserTProfile("Profile_PtGen_15", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_PtGenjet_Ptjet1hp );
      }
      if(	fabs( parton1.Eta() ) >= 1.5 && fabs( parton1.Eta() ) < 2.){ 
	CreateAndFillUserTProfile("Profile_PtGen_2", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_PtGenjet_Ptjet1hp ); 
      }
      if(	fabs( parton1.Eta() ) >= 2. && fabs( parton1.Eta() ) < 2.5){ 
	CreateAndFillUserTProfile("Profile_PtGen_25", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_PtGenjet_Ptjet1hp );
      }
    }
    */


     /*     
     double deltaRmin_Jet_parton1 = 1000;
     double deltaR_currentJet_parton1 = 1000;
     int IdJet_min = -1000;
     
     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ 
       
       //if(  fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
       //   && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
       //   && jetPtAK4->at(ijet) > getPreCutValue1("ptCut") ){
	       
	   TLorentzVector currentJet;
	   currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   

	   deltaR_currentJet_parton1 = currentJet.DeltaR(parton1);
	   
	   if( deltaR_currentJet_parton1 < deltaRmin_Jet_parton1){
	     
	     deltaRmin_Jet_parton1 = deltaR_currentJet_parton1;
	     IdJet_min = ijet;
	   }
	   //}
     }	   
     
     if(IdJet_min < -1) continue; // qui ci sono i tagli in cinematica
     
     CreateAndFillUserTH1D("H_DeltaRmin_Jet_parton1", 100, 0. , 10. , deltaRmin_Jet_parton1 );
     CreateAndFillUserTH1D("H_Id_min_Jet", 16, -0.5 , 15.5 , IdJet_min );
     
     if(deltaRmin_Jet_parton1 <= 0.3){
       //Calcolo ratio con un taglio sensato in DeltaR
       
       double     R_Ptjet_PtParton1;
       
       R_Ptjet_PtParton1 = jetPtAK4->at(IdJet_min) / parton1.Pt() ;  	
       
       CreateAndFillUserTH1D("H_R_Ptjet_PtParton1", 100, 0. , 10. , R_Ptjet_PtParton1 );
       CreateAndFillUserTProfile("Profile_PtReco", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_Ptjet_PtParton1 );  //, Double_t weight)	 
       
       if(	fabs( parton1.Eta() ) >= 0. && fabs( parton1.Eta() ) < 0.5 ){ 
	 CreateAndFillUserTProfile("Profile_PtReco_05", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_Ptjet_PtParton1 );  //, Double_t weight)	 
       }
       if(	fabs( parton1.Eta() ) >= 0.5 && fabs( parton1.Eta() ) < 1. ){ 
	 CreateAndFillUserTProfile("Profile_PtReco_1", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_Ptjet_PtParton1 );  //, Double_t weight)	 
       }
       if(	fabs( parton1.Eta() ) >= 1. && fabs( parton1.Eta() ) < 1.5 ){ 
	 CreateAndFillUserTProfile("Profile_PtReco_15", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_Ptjet_PtParton1 );  //, Double_t weight)	 
       }
       if(	fabs( parton1.Eta() ) >= 1.5 && fabs( parton1.Eta() ) < 2. ){ 
	 CreateAndFillUserTProfile("Profile_PtReco_2", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_Ptjet_PtParton1 );  //, Double_t weight)	 
       }
       if(	fabs( parton1.Eta() ) >= 2. && fabs( parton1.Eta() ) < 2.5 ){ 
	 CreateAndFillUserTProfile("Profile_PtReco_25", 1000, 0. , 10000. , 0., 3., parton1.Pt(), R_Ptjet_PtParton1 );  //, Double_t weight)	 
       }
     }		   
     */	



/*      // mi calcolo il widejet partendo dai partoni con le funzioni di smearing
    
    double parton_Pt[2];
    double parton_Eta[2];
    double parton_Phi[2];
    double parton_M[2];    
    parton_Pt[0]   = parton1.Pt();
    parton_Eta[0] = parton1.Eta();
    parton_Phi[0] = parton1.Phi();
    parton_M[0]   = parton1.M();

    parton_Pt[1]   = parton2.Pt();
    parton_Eta[1] = parton2.Eta();
    parton_Phi[1] = parton2.Phi();
    parton_M[1]   = parton2.M();    

    if(DeltaR[0] > 0.3 || DeltaR[1] > 0.3) cout<<"Uno dei due jet non matchea il partone. Butto evento"<<endl;
    if(DeltaR[0] > 0.3 || DeltaR[1] > 0.3) continue;
    
    CreateAndFillUserTH1D("H_DeltaR1_Match", 100, 0. , 3. , DeltaR[0]);
    CreateAndFillUserTH1D("H_DeltaR2_Match", 100, 0. , 3. , DeltaR[1]);

    double factor;
    double JetSmeared_Pt[2]   ={-1001, -1001};
    double JetSmeared_Eta[2] ={-1001, -1001};
    double JetSmeared_Phi[2] ={-1001, -1001};
    double JetSmeared_M[2]   ={-1001, -1001};
    double Massjet = 50; //GeV

    TLorentzVector JetSmeared1_tmp, JetSmeared2_tmp;
    TLorentzVector JetSmeared1, JetSmeared2;

    cout<<" " <<endl;
    cout<<"Costruisco i miei jet smeared"<<endl; 

    for(int kk = 0 ; kk < 2 ; kk++){

      cout<<"Parton "<< kk+1<<endl;
      cout<<"Pt: "<<parton_Pt[kk]<<"         |Eta|: "<< fabs(parton_Eta[kk])<<endl;

      CreateAndFillUserTH1D("H_DeltaR_Totale", 100, 0. , 3. , DeltaR[kk]);

	for(int ii=0; ii<5; ii++){
	  double eta_bin_min = ii/2.;       
	  double eta_bin_max = ii/2. +0.5;

	  if( fabs(parton_Eta[kk])>=eta_bin_min && fabs(parton_Eta[kk]) < eta_bin_max){
	    cout<<"Passed bin in eta"<<endl;	
	    cout << "eta_bin: "<< eta_bin_min <<" - "<<eta_bin_max << endl;
	    
	    for(int jj=0; jj<n_bin ; jj++){	 
	      int pt_bin_min = ( (jj-1)*step_pt)+step_pt;
	      int pt_bin_max = (jj*step_pt)+step_pt;
	      
	      if(parton_Pt[kk] >=pt_bin_min && parton_Pt[kk] < pt_bin_max){ 
		cout<<"Passed bin in pt"<<endl;
		cout << "pt_bin: "<< pt_bin_min <<" - "<<pt_bin_max<< endl;

		sprintf(HistoName,"Histo_RRecowide_Parton_%d_%d",ii,pt_bin_max);	    			
		cout<<"Estraggo il fattore dall'istogramma: "<< HistoName << endl;	      
		
		h[ii][jj]= (TH1D*)file1.Get(HistoName);		
		//cout<<"h[ii][jj]-> Integral(): "<<h[ii][jj]->Integral()<<endl;
		
		if( h[ii][jj]->Integral() == 0){
		  cout<<HistoName<<" ha integrale nullo"<<endl;
		  factor = 1;
		}else{
		  factor = h[ii][jj]->GetRandom();
		}

		cout<<"Fattore estratto: "<<factor<<endl;
		
		JetSmeared_Pt[kk] = parton_Pt[kk] * factor ;
		JetSmeared_Eta[kk] = parton_Eta[kk] ;
		JetSmeared_Phi[kk] = parton_Phi[kk] ;
		JetSmeared_M[kk] = Massjet;
		
		cout<<"Jet Smeared "<<kk+1<<endl;
		cout<<"Pt: "<<JetSmeared_Pt[kk]<<"          Eta: "<<JetSmeared_Eta[kk]<<"          Phi: "<<JetSmeared_Phi[kk]<<"          M: "<<JetSmeared_M[kk]<<endl;
		
	      }//if range di pt
	    }//for sui bin in pt
	  }//if range di eta
	}//for sui bin in eta
    }//for sui due partoni
    
    //se non trova il match    
    if(JetSmeared_Pt[0] <0 || JetSmeared_Pt[1] <0) cout<<"Non ci sono due jet smeared. Butto evento"<<endl;;
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
*/
