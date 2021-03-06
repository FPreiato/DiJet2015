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
#include <TMath.h>


// Variables to change
int step_pt=200;       // GeV for pt binning
int pT_max = 4800;  // maximum value of pT range [GeV]
int n_categories= 4;  // number of categories for the smearing functions: 2 || 4 || 6

// Histogram parameters : range and binning
int xmin = 0.0 ;
int xmax = 2.0 ;
int histo_bin = 100; 
int n_bin;
char HistoName[200];

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

   // Check variables
   if( step_pt <0){
     cout<<" "<<endl;
     cout<<"Wrong pT binning"<<endl;
     exit(1);
   }   
   if( n_categories != 2 && n_categories != 4 && n_categories != 6){
     cout<<" "<<endl;
     cout<<"Wrong number of categories"<<endl;
     cout<<"Choose between 2, 4 or 6"<<endl;
     exit(1);
   }   

   //create the histograms for the smearing functions
   for(int ii=0; ii<5; ii++){
     n_bin = pT_max/ step_pt;
     for(int jj=0; jj<n_bin ; jj++){	 
       int pt_bin_max = (jj*step_pt)+step_pt;
       sprintf(HistoName,"Histo_RGluons_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       if(n_categories == 2){// all quarks together
       sprintf(HistoName,"Histo_RQuarks_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       }
       if(n_categories == 4){// quark light, charm e bottom
       sprintf(HistoName,"Histo_RQuarkLIGHT_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       sprintf(HistoName,"Histo_RQuarkCHARM_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       sprintf(HistoName,"Histo_RQuarkBOTTOM_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       }
       if(n_categories == 6){// up, down, strange, charm, bottom
       sprintf(HistoName,"Histo_RQuarkDOWN_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       sprintf(HistoName,"Histo_RQuarkUP_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       sprintf(HistoName,"Histo_RQuarkSTRANGE_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       sprintf(HistoName,"Histo_RQuarkCHARM_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       sprintf(HistoName,"Histo_RQuarkBOTTOM_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
       }
     }
   }

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
     
     CreateAndFillUserTH1D("H_step_pt", 1000, 0. , 1000. , step_pt );
     CreateAndFillUserTH1D("H_n_categories", 7, -0.5 , 6.5 , n_categories );
     CreateAndFillUserTH1D("H_pT_max", 5000, 0. , 5000 , pT_max );

     // Parton Level: Define partons from the resonance
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
     
     // cout<<"PARTON 1"<<endl;
     // cout<<"Pt: "<<parton1.Pt()<< " Eta: "<<parton1.Eta()<<" Phi: "<<parton1.Phi()<<" M: "<<parton1.M()<<" pdgId: "<<parton1_pdgId<<endl;
     // cout<<"PARTON 2"<<endl;
     // cout<<"Pt: "<<parton2.Pt()<< " Eta: "<<parton2.Eta()<<" Phi: "<<parton2.Phi()<<" M: "<<parton2.M()<<" pdgId: "<<parton2_pdgId<<endl;
     
     // Create dijet(reco) system
     diparton = parton1 + parton2;
     
     //+++++++++++++++++++++++++++++++++++++++++++++++++++ 
     // Gen-level : Construct the GenWidejet from genjet_ak4
     size_t no_Genjets_ak4 = jetPtGenAK4->size();
     TLorentzVector Genjet1, Genjet2;
     TLorentzVector Genwj1_tmp, Genwj2_tmp;
     TLorentzVector Genwj1, Genwj2;
     TLorentzVector Genwdijet;
     double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     if(no_Genjets_ak4 < 2) continue;
     
     //TLorentzVector leading Genjet
     Genjet1.SetPtEtaPhiM(jetPtGenAK4->at(0), jetEtaGenAK4->at(0), jetPhiGenAK4->at(0), jetMassGenAK4->at(0) );
     Genjet2.SetPtEtaPhiM(jetPtGenAK4->at(1), jetEtaGenAK4->at(1), jetPhiGenAK4->at(1), jetMassGenAK4->at(1) );
     
     double n_genjet = 0;
     
     for( Long64_t ijet=0;  ijet<no_Genjets_ak4;  ijet++){//genjet loop for ak4
       TLorentzVector currentGenJet;
       currentGenJet.SetPtEtaPhiM(jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet) );   
       
       if( currentGenJet.Pt() > 30 && fabs(currentGenJet.Eta())< 2.5){
	 n_genjet = n_genjet +1;
       }
       
       double DeltaR1 = currentGenJet.DeltaR(Genjet1);
       double DeltaR2 = currentGenJet.DeltaR(Genjet2);
       
       if(DeltaR1 < DeltaR2 && DeltaR1 < GenwideJetDeltaR_){
	 Genwj1_tmp += currentGenJet;
       }
       else if(DeltaR2 < GenwideJetDeltaR_){
	 Genwj2_tmp += currentGenJet;
       }
     }//end of ak4 loop
     
     CreateAndFillUserTH1D("H_Number_Genjets", 31, -0.5 , 30.5 , n_genjet);
     
     if( fabs(Genwj1_tmp.Eta() ) > 2.5 || fabs(Genwj2_tmp.Eta() ) > 2.5) continue;
     if( Genwj1_tmp.Pt() <30 || Genwj2_tmp.Pt() <30) continue;
     
     double DeltaR_GenWidejet1_tmp_parton1 = Genwj1_tmp.DeltaR(parton1);     
     double DeltaR_GenWidejet2_tmp_parton1 = Genwj2_tmp.DeltaR(parton1);     
     
     if(DeltaR_GenWidejet1_tmp_parton1 < DeltaR_GenWidejet2_tmp_parton1){
       Genwj1 = Genwj1_tmp;
       Genwj2 = Genwj2_tmp;
     }
     else{
       Genwj1 = Genwj2_tmp; 
       Genwj2 = Genwj1_tmp; 
     }
     
     double DeltaR_GenWideJet1_parton1;   
     double DeltaR_GenWideJet2_parton2;   
     DeltaR_GenWideJet1_parton1 = Genwj1.DeltaR(parton1);
     DeltaR_GenWideJet2_parton2 = Genwj2.DeltaR(parton2);
     
     CreateAndFillUserTH1D("H_DeltaR_GenWideJet1_Parton1", 100, 0. , 2. , DeltaR_GenWideJet1_parton1);
     CreateAndFillUserTH1D("H_DeltaR_GenWideJet2_Parton2", 100, 0. , 2. , DeltaR_GenWideJet2_parton2);     
     
     // cout<<"GenWideJet 1"<<endl;
     // cout<< "Pt: "<<Genwj1.Pt()<<"         Eta: "<< Genwj1.Eta()<<"         Phi: "<<Genwj1.Phi()<<"         M: "<<Genwj1.M()<<endl;
     // cout<<"GenWideJet 2"<<endl;
     // cout<< "Pt: "<<Genwj2.Pt()<<"         Eta: "<< Genwj2.Eta()<<"         Phi: "<<Genwj2.Phi()<<"         M: "<<Genwj2.M()<<endl;
          
     // Create dijet system
     Genwdijet = Genwj1 + Genwj2;
     
     /////////////////////////////////////////////////////////
     // Reco level : Construct  widejet from recojet_ak4
     size_t no_jets_ak4=jetPtAK4->size();    
     TLorentzVector jet1, jet2;
     TLorentzVector wj1_tmp, wj2_tmp;
     TLorentzVector wj1, wj2; 
     TLorentzVector wdijet;
     double wideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     if(no_jets_ak4<2) continue;
     
     // TLorentzVector jet1, jet2 leading reco jet
     jet1.SetPtEtaPhiM(jetPtAK4->at(0), jetEtaAK4->at(0), jetPhiAK4->at(0), jetMassAK4->at(0));
     jet2.SetPtEtaPhiM(jetPtAK4->at(1), jetEtaAK4->at(1), jetPhiAK4->at(1), jetMassAK4->at(1));
     
     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
       TLorentzVector currentJet;
       currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet), jetEtaAK4->at(ijet), jetPhiAK4->at(ijet), jetMassAK4->at(ijet));   
       
       double DeltaR1 = currentJet.DeltaR(jet1);
       double DeltaR2 = currentJet.DeltaR(jet2);
       
       if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
	 wj1_tmp += currentJet;
       }
       else if(DeltaR2 < wideJetDeltaR_){
	 wj2_tmp += currentJet;
       }			 
     } //end of ak4 jet loop		     
     
     if( wj1_tmp.Pt()<0 || wj2_tmp.Pt()<0 ) continue;
     
     double DeltaR_Widejet1_tmp_parton1 = wj1_tmp.DeltaR(parton1);     
     double DeltaR_Widejet2_tmp_parton1 = wj2_tmp.DeltaR(parton1);     
     
     if(DeltaR_Widejet1_tmp_parton1 < DeltaR_Widejet2_tmp_parton1){
       wj1 = wj1_tmp;
       wj2 = wj2_tmp;
     }
     else{
       wj1 = wj2_tmp; 
       wj2 = wj1_tmp; 
     }
     
     // check the spatial resolution Widejet-Parton
     // double Eta_difference_Widejet1_parton1 = wj1.Eta() - parton1.Eta();
     // double Phi_difference_Widejet1_parton1 = wj1.Phi() - parton1.Phi() ;
     // double DeltaR_WideJet1_parton1            = wj1.DeltaR(parton1);
     // double Eta_difference_Widejet2_parton2 = wj2.Eta() - parton2.Eta() ;
     // double Phi_difference_Widejet2_parton2 = wj2.Phi() - parton2.Phi() ;
     // double DeltaR_WideJet2_parton2            = wj2.DeltaR(parton2);
     // double Angle_difference_Widejet1_parton1 = wj1.Angle(parton1.Vect());
     // double Angle_difference_Widejet2_parton2= wj2.Angle(parton2.Vect());
     // double Angle_difference_Partons= parton1.Angle(parton2.Vect());
     
     // CreateAndFillUserTH1D("H_Eta_difference_Widejet1_parton1", 100, -2. , 2. , Eta_difference_Widejet1_parton1);
     // CreateAndFillUserTH1D("H_Phi_difference_Widejet1_parton1", 100, -2. , 2. , Phi_difference_Widejet1_parton1);
     // CreateAndFillUserTH1D("H_Eta_difference_Widejet2_parton2", 100, -2. , 2. , Eta_difference_Widejet2_parton2);
     // CreateAndFillUserTH1D("H_Phi_difference_Widejet2_parton2", 100, -2. , 2. , Phi_difference_Widejet2_parton2);       
     // CreateAndFillUserTH1D("H_DeltaR_WideJet1_parton1", 100, 0. , 3. , DeltaR_WideJet1_parton1);
     // CreateAndFillUserTH1D("H_DeltaR_WideJet2_parton2", 100, 0. , 3. , DeltaR_WideJet1_parton1);
     // CreateAndFillUserTH1D("H_Angle_difference_WideJet1_parton1", 100, -3.15 , 3.15, Angle_difference_Widejet1_parton1 );
     // CreateAndFillUserTH1D("H_Angle_difference_WideJet2_parton2", 100, -3.15 , 3.15, Angle_difference_Widejet2_parton2 );
     // CreateAndFillUserTH1D("H_Angle_difference_Partons", 100, -3.15 , 3.15, Angle_difference_Partons );
     
     // cout<<"WideJet 1"<<endl;
     // cout<< "Pt: "<<wj1.Pt()<<"         Eta: "<< wj1.Eta()<<"         Phi: "<<wj1.Phi()<<"         M: "<<wj1.M()<<endl;
     // cout<<"WideJet 2"<<endl;
     // cout<< "Pt: "<<wj2.Pt()<<"         Eta: "<< wj2.Eta()<<"         Phi: "<<wj2.Phi()<<"         M: "<<wj2.M()<<endl;
     
     // Create dijet system
     wdijet = wj1 + wj2;

     //+++++++++++++++++++++++++++++++++++++++++++++++
     // TOOLS For SMEARING FUNCTIONS
     double DeltaR_WideJet1_GenWideJet1 = wj1.DeltaR(Genwj1);     
     double DeltaR_WideJet2_GenWideJet2 = wj2.DeltaR(Genwj2);         
     CreateAndFillUserTH1D("H_DeltaR_WideJet1_GenWideJet1", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
     CreateAndFillUserTH1D("H_DeltaR_WideJet2_GenWideJet2", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
     
     // Calculate response and resolution
     double R_PtWidejet_PtGenWideJet[2];
     double Diff_EtaWidejet_EtaGenWideJet[2];
     double Diff_PhiWidejet_PhiGenWideJet[2];

     R_PtWidejet_PtGenWideJet[0]        = wj1.Pt() / Genwj1.Pt() ;
     R_PtWidejet_PtGenWideJet[1]        = wj2.Pt() / Genwj2.Pt() ;
     Diff_EtaWidejet_EtaGenWideJet[0] = wj1.Eta() - Genwj1.Eta() ;
     Diff_EtaWidejet_EtaGenWideJet[1] = wj2.Eta() - Genwj2.Eta() ;
     Diff_PhiWidejet_PhiGenWideJet[0] = wj1.Phi() - Genwj1.Phi() ;
     Diff_PhiWidejet_PhiGenWideJet[1] = wj2.Phi() - Genwj2.Phi() ;

     double Parton_pdgId[2];
     double GenWideJet_Pt[2];
     double GenWideJet_Eta[2];    
     GenWideJet_Pt[0]   = Genwj1.Pt();
     GenWideJet_Eta[0] = Genwj1.Eta();
     Parton_pdgId[0]      = parton1_pdgId;        
     GenWideJet_Pt[1]   = Genwj2.Pt();
     GenWideJet_Eta[1] = Genwj2.Eta();    
     Parton_pdgId[1]      = parton2_pdgId;    
     
     // geometrical matching
     if(DeltaR_WideJet1_GenWideJet1 > 0.3 || DeltaR_WideJet2_GenWideJet2 > 0.3) continue;
     
     for(int kk=0; kk<2; kk++){// loop on GenWideJets
          
       CreateAndFillUserTH1D("H_Diff_EtaWidejet_EtaGenWideJet", 100, -1, 1, Diff_EtaWidejet_EtaGenWideJet[kk] );
       CreateAndFillUserTH1D("H_Diff_PhiWidejet_PhiGenWideJet", 100, -1, 1, Diff_PhiWidejet_PhiGenWideJet[kk] );
       
       for(int ii=0; ii<5; ii++){ // loop on eta bins	  
	 double eta_bin_min = ii/2.;       
	 double eta_bin_max = ii/2. +0.5;	 
	 if( fabs(GenWideJet_Eta[kk]) >=eta_bin_min && fabs(GenWideJet_Eta[kk]) < eta_bin_max){	    

	   for(int jj=0; jj<n_bin ; jj++){ // loop on pT bins	 	      
	     int pt_bin_min = ( (jj-1)*step_pt)+step_pt;
	     int pt_bin_max = (jj*step_pt)+step_pt;
	     if(GenWideJet_Pt[kk] >=pt_bin_min && GenWideJet_Pt[kk] < pt_bin_max){ 
	       
	       // Categories 
	       if(n_categories == 2){ // Quarks // Gluons
		 if(Parton_pdgId[kk] == 21){ 
		   sprintf(HistoName,"Histo_RGluons_%d_%d",ii,pt_bin_max);
		 }
		 else{
		   sprintf(HistoName,"Histo_RQuarks_%d_%d",ii,pt_bin_max);
		 }
	       }
	       if(n_categories == 4){ // Light Quarks, Charm, Bottom // Gluons	       
		 if(Parton_pdgId[kk] == 21){ 
		   sprintf(HistoName,"Histo_RGluons_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 1 || fabs(Parton_pdgId[kk]) == 2  || fabs(Parton_pdgId[kk]) == 3 ){ 
		   sprintf(HistoName,"Histo_RQuarkLIGHT_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 4){ 
		   sprintf(HistoName,"Histo_RQuarkCHARM_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 5){ 
		   sprintf(HistoName,"Histo_RQuarkBOTTOM_%d_%d",ii,pt_bin_max);
		 }
	       }
	       if(n_categories == 6){ // Up, Down, Strange, Charm, Bottom // Gluons
		 if(Parton_pdgId[kk] == 21){ 
		   sprintf(HistoName,"Histo_RGluons_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 1){ 
		   sprintf(HistoName,"Histo_RQuarkDOWN_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 2){ 
		   sprintf(HistoName,"Histo_RQuarkUP_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 3){ 
		   sprintf(HistoName,"Histo_RQuarkSTRANGE_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 4){ 
		   sprintf(HistoName,"Histo_RQuarkCHARM_%d_%d",ii,pt_bin_max);
		 }
		 if( fabs(Parton_pdgId[kk]) == 5){ 
		   sprintf(HistoName,"Histo_RQuarkBOTTOM_%d_%d",ii,pt_bin_max);
		 }
	       }

	       FillUserTH1D(HistoName, R_PtWidejet_PtGenWideJet[kk] );
	       
	     }//end if on parton_Pt
	   }//end loop on pT bins
	 }//end if on parton_Eta
       }//end loop on eta bins
     }//end loop on GenWidejet
     
     ///////////////////////////////////////////////////////////////////////////////////////////
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

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedAllPreviousCuts("dijetWide_M") && passedCut("dijetWide_M") ) 
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
	 
       }
   
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
   
     ////////////////////// User's code ends here ///////////////////////
   
   } // End loop over events   
   

   // //pT of both jets, to be built using the histograms produced automatically by baseClass
   // TH1F * h_pTJets = new TH1F ("h_pTJets","", getHistoNBins("pT1stJet"), getHistoMin("pT1stJet"), getHistoMax("pT1stJet"));
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT1stJet") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   // h_pTJets->Add( & getHisto_noCuts_or_skim("pT2ndJet") );
   // //one could also do:  *h_pTJets = getHisto_noCuts_or_skim("pT1stJet") + getHisto_noCuts_or_skim("pT2ndJet");
   // h_pTJets->Write();
   // //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
