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


inline double delta_phi(double phi1, double phi2) {
  double dphi = TMath::Abs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
  
}

inline double delta_eta(double eta1, double eta2) {
  return (eta2 >= 0 ? eta1 - eta2 : eta2 - eta1);
}


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

     // Generated jet from hard process

     TLorentzVector p1, p2;
     TLorentzVector parton1_tmp, parton2_tmp;
     TLorentzVector parton1, parton2;
     TLorentzVector diparton;

     double diparton_MJJ;
     double diparton_DeltaEtaJJ;
     double diparton_DeltaPhiJJ;

     double parton1_pt;
     double parton1_eta;
     double parton1_phi;
     double parton1_mass;
     double parton1_pdgId;

     double parton2_pt;
     double parton2_eta;
     double parton2_phi;
     double parton2_mass;
     double parton2_pdgId;

     //++++++++++++++++++++++++++++++++++++++++++++++++

     //       cout<< "size: " << gen_pt->size() << endl;
     //for(int ii= 0 ; ii < gen_pt->size(); ii++){
     //cout<<"pT "<< ii<<"     "<<gen_pt->at(ii) <<endl;
     //}

     if(gen_pt->size() == 3){

       p1.SetPxPyPzE(gen_px->at(1) , gen_py->at(1) , gen_pz->at(1), gen_energy->at(1) );
       parton1_pt       = gen_pt->at(1);
       parton1_eta     = gen_eta->at(1);
       parton1_phi     = gen_phi->at(1);
       parton1_mass  = p1.M();
       parton1_pdgId = gen_pdgId -> at(1);

       p2.SetPxPyPzE(gen_px->at(2) , gen_py->at(2) , gen_pz->at(2), gen_energy->at(2) );
       parton2_pt       =  gen_pt->at(2);
       parton2_eta      = gen_eta->at(2);
       parton2_phi     = gen_phi->at(2);
       parton2_mass   = p2.M();
       parton2_pdgId = gen_pdgId -> at(2);
     }

     if(gen_pt->size() == 4){

       p1.SetPxPyPzE(gen_px->at(2) , gen_py->at(2) , gen_pz->at(2), gen_energy->at(2) );
       parton1_pt       = gen_pt->at(2);
       parton1_eta      = gen_eta->at(2);
       parton1_phi      = gen_phi->at(2);
       parton1_mass   = p1.M();
       parton1_pdgId = gen_pdgId -> at(2);

       p2.SetPxPyPzE(gen_px->at(3) , gen_py->at(3) , gen_pz->at(3), gen_energy->at(3) );
       parton2_pt       = gen_pt->at(3);
       parton2_eta     = gen_eta->at(3);
       parton2_phi     = gen_phi->at(3);
       parton2_mass  = p2.M();
       parton2_pdgId = gen_pdgId -> at(3);     
     }

     if(gen_pt->size() == 5){

       p1.SetPxPyPzE(gen_px->at(3) , gen_py->at(3) , gen_pz->at(3), gen_energy->at(3) );
       parton1_pt       = gen_pt->at(3);
       parton1_eta     = gen_eta->at(3);
       parton1_phi     = gen_phi->at(3);
       parton1_mass  = p1.M();
       parton1_pdgId = gen_pdgId -> at(3);

       p2.SetPxPyPzE(gen_px->at(4) , gen_py->at(4) , gen_pz->at(4), gen_energy->at(4) );
       parton2_pt       = gen_pt->at(4);
       parton2_eta     = gen_eta->at(4);
       parton2_phi     = gen_phi->at(4);
       parton2_mass  = p2.M();
       parton2_pdgId = gen_pdgId -> at(4);   
     }

     //     cout << "parton1_pdgId: " << parton1_pdgId << endl;
     //     cout << "parton2_pdgId: " << parton2_pdgId << endl;

     //     cout << "parton1_mass: " << parton1_mass << endl;
     //     cout << "parton2_mass: " << parton2_mass << endl;

     parton1_tmp.SetPtEtaPhiM(parton1_pt, parton1_eta, parton1_phi, parton1_mass );
     parton2_tmp.SetPtEtaPhiM(parton2_pt, parton2_eta, parton2_phi, parton2_mass );

     //     if(parton1_mass <0 || parton2_mass <0) cout << "++++++++++++++++++ WARNING!!!" << endl;

     // Re-order the genjets in pt
     if( parton1_tmp.Pt() > parton2_tmp.Pt() ){
       parton1 = parton1_tmp;
       parton2 = parton2_tmp;
     }else{
       parton1 = parton2_tmp;
       parton2 = parton1_tmp;
     }
     
     if( parton1.Pt()>0 && parton2.Pt()>0 ){ 
       if( fabs(parton1.Eta() ) < 2.5 && fabs(parton2.Eta() ) < 2.5){
       
       // Create dijet(reco) system
       diparton = parton1 + parton2;
       diparton_MJJ = diparton.M();
       diparton_DeltaEtaJJ = fabs( parton1.Eta() - parton2.Eta() );
       diparton_DeltaPhiJJ = fabs(parton1.DeltaPhi(parton2));
       
       //Oggetti che sommo
       CreateAndFillUserTH1D("H_parton1_Pt", 1000, 0. , 10000. , parton1.Pt() );
       CreateAndFillUserTH1D("H_parton1_Eta", 100, -4. , 4. , parton1.Eta() );
       CreateAndFillUserTH1D("H_parton1_Phi", 100, -4. , 4. , parton1.Phi() );
       CreateAndFillUserTH1D("H_parton1_M", 50000, 0. , 10000. , parton1.M() );
       CreateAndFillUserTH1D("H_parton2_Pt", 1000, 0. , 10000. , parton2.Pt() );
       CreateAndFillUserTH1D("H_parton2_Eta", 100, -4. , 4. , parton2.Eta() );  
       CreateAndFillUserTH1D("H_parton2_Phi", 100, -4. , 4. , parton2.Phi() );
       CreateAndFillUserTH1D("H_parton2_M", 50000, 0. , 10000. , parton2.M() );
       
       //Variabili del sistema diparton
       CreateAndFillUserTH1D("H_diparton_Pt", 1000, 0. , 10000. , diparton.Pt() );
       CreateAndFillUserTH1D("H_diparton_Eta", 100, -4. , 4. , diparton.Eta() );
       CreateAndFillUserTH1D("H_diparton_Phi", 100, -4. , 4. , diparton.Phi() );
       CreateAndFillUserTH1D("H_diparton_M", 1000, 0. , 10000. , diparton_MJJ );

       CreateAndFillUserTH1D("H_diparton_DeltaEtaJJ", 100, 0. , 8. , diparton_DeltaEtaJJ );
       
       CreateAndFillUserTH2D("H_diparton_Pt_Eta", 1000, 0., 10000., 120, -10. , 10. , diparton.Pt(), diparton.Eta() );
       
       // Put widejets in the container -> see
       //Genjets.push_back( Genjet1 );
       //Genjets.push_back( Genjet2 );
       }
     }  
     //     cout<< "diparton_MJJ:  " << diparton_MJJ<< endl;
     
     //+++++++++++++++++++++++++++++++++++++++++++++++++++

     // Generated jet from parton shower
     
     size_t no_Genjets_ak4 = jetPtGenAK4->size();
     
     double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");

     double GenMJJ; 
     double GenDeltaEtaJJ ;
     double GenDeltaPhiJJ ;     
     
     double GenMJJWide; 
     double GenDeltaEtaJJWide;
     double GenDeltaPhiJJWide;
     
     TLorentzVector Genjet1, Genjet2, Gendijet;
     vector<TLorentzVector> Genjets;
     TLorentzVector Genwj1_tmp, Genwj2_tmp;
     TLorentzVector Genwj1, Genwj2, Genwdijet;
    vector<TLorentzVector> Genwidejets;
 
    
     if(no_Genjets_ak4 < 2) cout << "++++++++++++++++ Pochi jets: "<< no_Genjets_ak4 << endl;
     
     if(no_Genjets_ak4>=2){
       
       Genjet1.SetPtEtaPhiM(jetPtGenAK4->at(0), jetEtaGenAK4->at(0), jetPhiGenAK4->at(0), jetMassGenAK4->at(0) );
       Genjet2.SetPtEtaPhiM(jetPtGenAK4->at(1), jetEtaGenAK4->at(1), jetPhiGenAK4->at(1), jetMassGenAK4->at(1) );

       //       cout<< "GenJet AK4_0 Mass: " << Genjet1.M()<<endl;
       //       cout<< "GenJet AK4_1 Mass: " << Genjet2.M()<<endl;
       
       CreateAndFillUserTH1D("H_Genjet1_Pt", 1000, 0. , 10000. , Genjet1.Pt() );
       CreateAndFillUserTH1D("H_Genjet1_Eta", 100, -4. , 4. , Genjet1.Eta() );
       CreateAndFillUserTH1D("H_Genjet1_Phi", 100, -4. , 4. , Genjet1.Phi() );
       CreateAndFillUserTH1D("H_Genjet1_M", 1000, 0. , 10000. , Genjet1.M() );

       CreateAndFillUserTH1D("H_Genjet2_Pt", 1000, 0. , 10000. , Genjet2.Pt() );
       CreateAndFillUserTH1D("H_Genjet2_Eta", 100, -4. , 4. , Genjet2.Eta() );
       CreateAndFillUserTH1D("H_Genjet2_Phi", 100, -4. , 4. , Genjet2.Phi() );
       CreateAndFillUserTH1D("H_Genjet2_M", 1000, 0. , 10000. , Genjet2.M() );


       if( Genjet1.Pt()>0 && Genjet2.Pt()>0 ){
	 // Create dijet(reco) system
	 Gendijet = Genjet1 + Genjet2;
	 GenMJJ = Gendijet.M();
	 GenDeltaEtaJJ = fabs(Genjet1.Eta()-Genjet2.Eta());
	 GenDeltaPhiJJ = fabs(Genjet1.DeltaPhi(Genjet2));
	 
	 CreateAndFillUserTH1D("H_Gendijet_Pt", 1000, 0. , 10000. , Gendijet.Pt() );
	 CreateAndFillUserTH1D("H_Gendijet_Eta", 100, -4. , 4. , Gendijet.Eta() );
	 CreateAndFillUserTH1D("H_Gendijet_Phi", 100, -4. , 4. , Gendijet.Phi() );
	 CreateAndFillUserTH1D("H_Gendijet_M", 1000, 0. , 10000. , GenMJJ );
	 
	 CreateAndFillUserTH1D("H_GenDeltaEtaJJ", 100, 0. , 8. , GenDeltaEtaJJ );
	 
	 // Put widejets in the container
	 //Genjets.push_back( Genjet1 );
	 //Genjets.push_back( Genjet2 );
       }
       
       //       cout << "GenMJJ:  " << GenMJJ << endl;
       
       //       cout<< "no_Genjets_ak4: " << no_Genjets_ak4 << endl;
       
       // Generated jet clustering -> wide genjet
       for( Long64_t ijet=0;  ijet<no_Genjets_ak4;  ijet++){
	 
	 TLorentzVector currentGenJet;
	 currentGenJet.SetPtEtaPhiM(jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet) );   

	 CreateAndFillUserTH1D("H_AllGenjet_Pt", 10000, 0. , 10000. , currentGenJet.Pt() );
	 CreateAndFillUserTH1D("H_AllGenjet_Eta", 2000, -10. , 10. , currentGenJet.Eta() );

	 
	 double DeltaR1 = currentGenJet.DeltaR(Genjet1);
	 double DeltaR2 = currentGenJet.DeltaR(Genjet2);
	 
	 if(DeltaR1 < DeltaR2){
	   CreateAndFillUserTH1D("H_DeltaR1_pre", 10, 0. , 10. , DeltaR1 );
	   if( DeltaR1 < GenwideJetDeltaR_){
	     
	     CreateAndFillUserTH1D("H_DeltaR1_post", 10, 0. , 10. , DeltaR1 );
	     CreateAndFillUserTH1D("H_PtGenJetAK4_enter1", 1000, 0. , 10000. , currentGenJet.Pt()); 
	     
	     Genwj1_tmp += currentGenJet;
	   }
	 }else{
	   CreateAndFillUserTH1D("H_DeltaR2_pre", 10, 0. , 10. , DeltaR2 );
	   if(DeltaR2 < GenwideJetDeltaR_){
	     CreateAndFillUserTH1D("H_DeltaR2_post", 10, 0. , 10. , DeltaR2 );
	     CreateAndFillUserTH1D("H_PtGenJetAK4_enter2", 1000, 0. , 10000. , currentGenJet.Pt()); 
	     
	     Genwj2_tmp += currentGenJet;
	   }			 
	 }
       }//end of ak4 jet loop		     
     }// end of two jets.
     
     // Re-order the wide genjets in pt
     if( Genwj1_tmp.Pt() > Genwj2_tmp.Pt()){
       Genwj1 = Genwj1_tmp;
       Genwj2 = Genwj2_tmp;
     }else{
       Genwj1 = Genwj2_tmp;
       Genwj2 = Genwj1_tmp;
     }
     
     //     cout<< "Genwj1.M() :  " << Genwj1.M()<< endl;
     //     cout<< "Genwj2.M() :  " << Genwj2.M()<< endl;
     
     if( Genwj1.Pt()>0 && Genwj2.Pt()>0 ){
       // Create dijet system
       Genwdijet = Genwj1 + Genwj2;
       GenMJJWide = Genwdijet.M();
       GenDeltaEtaJJWide = fabs(Genwj1.Eta() - Genwj2.Eta());
       GenDeltaPhiJJWide = fabs(Genwj1.DeltaPhi(Genwj2));
       
       //Oggetti che sommo
       CreateAndFillUserTH1D("H_GenWideJet1_Pt", 1000, 0. , 10000. , Genwj1.Pt()); 
       CreateAndFillUserTH1D("H_GenWideJet1_Eta", 100, -4. , 4. , Genwj1.Eta()); 
       CreateAndFillUserTH1D("H_GenWideJet1_Phi", 100, -4. , 4. , Genwj1.Phi());
       CreateAndFillUserTH1D("H_GenWideJet1_M", 1000, 0. , 10000. , Genwj1.M());  
       CreateAndFillUserTH1D("H_GenWideJet2_Pt", 1000, 0. , 10000. , Genwj2.Pt()); 
       CreateAndFillUserTH1D("H_GenWideJet2_Eta", 100, -4. , 4. , Genwj2.Eta()); 
       CreateAndFillUserTH1D("H_GenWideJet2_Phi", 100, -4. , 4. , Genwj2.Phi());
       CreateAndFillUserTH1D("H_GenWideJet2_M", 1000, 0. , 10000. , Genwj2.M());  
       
       //Variabili del sistema Dijet Wide	   
       CreateAndFillUserTH1D("H_GendijetWide_Pt", 1000, 0. , 10000. , Genwdijet.Pt() );
       CreateAndFillUserTH1D("H_GendijetWide_Eta", 100, -4. , 4. , Genwdijet.Eta() );
       CreateAndFillUserTH1D("H_GendijetWide_Phi", 100, -4. , 4. , Genwdijet.Phi() );
       CreateAndFillUserTH1D("H_GendijetWide_M", 1000, 0. , 10000. , GenMJJWide );

       CreateAndFillUserTH1D("H_GenDeltaEtaJJWide", 100, 0. , 8. , GenDeltaEtaJJWide );
       
       CreateAndFillUserTH2D("H_PtGenWideJJ_EtaGenWideJJ", 1000, 0., 10000., 120, -10. , 10. , Genwdijet.Pt(), Genwdijet.Eta() );
       
       //Confronto MJJ wide in funzione del Pt e Eta dei jet piu energetici
       CreateAndFillUserTH2D("H_GenMJJWide_PtGenjet1", 1000, 0., 10000., 1000, 0. , 10000. , Genjet1.Pt(), GenMJJWide );
       CreateAndFillUserTH2D("H_GenMJJWide_PtGenjet2", 1000, 0., 10000., 1000, 0. , 10000. , Genjet2.Pt(), GenMJJWide );
       CreateAndFillUserTH2D("H_GenMJJWide_EtaGenjet1", 120, -20., 20., 1000, 0. , 10000. , Genjet1.Eta(), GenMJJWide );
       CreateAndFillUserTH2D("H_GenMJJWide_EtaGenjet2", 120, -20., 20., 1000, 0. , 10000. , Genjet2.Eta(), GenMJJWide );

       // Put widejets in the container
       //       Genwidejets.push_back( Genwj1 );
       //Genwidejets.push_back( Genwj2 );

       //     cout<< "GenMJJWide :  " << GenMJJWide<< endl;       

     }
   
     //////////////////////////////////////////////////////////////////////////////////////////

     /*
     double deltaR_currentGenjet_parton1 = 1000;
     double deltaRmin_Genjet_parton1 = 1000;
     int Id1_min = -1000;

     double deltaR_currentGenjet_parton2 = 1000;     
     double deltaRmin_Genjet_parton2 = 1000;
     int Id2_min = -1000;

     if(no_Genjets_ak4>=2){

       cout<< " " << endl;
       cout<< "COMINCIA QUI"<<endl;

       for(Long64_t ijet=0; ijet<no_Genjets_ak4; ijet++)
	 { 
	   
	   TLorentzVector currentGenJet;
	   currentGenJet.SetPtEtaPhiM( jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet));   
	   
	   deltaR_currentGenjet_parton1 = currentGenJet.DeltaR(parton1);
	   deltaR_currentGenjet_parton2 = currentGenJet.DeltaR(parton2);
	  
	   cout<< "deltaR_currentGenjet_parton1: "<< deltaR_currentGenjet_parton1 << endl;
	   cout<< "deltaR_currentGenjet_parton2: "<< deltaR_currentGenjet_parton2 << endl;

 
	   if(deltaR_currentGenjet_parton1 < deltaR_currentGenjet_parton2){
	     if( deltaR_currentGenjet_parton1 < deltaRmin_Genjet_parton1){
	     
	     deltaRmin_Genjet_parton1 = deltaR_currentGenjet_parton1;
	     Id1_min = ijet;

	     deltaRmin_Genjet_parton2 = deltaR_currentGenjet_parton2;
	     Id2_min = ijet;

	     cout<< "DeltaR1 minimo: "<< deltaRmin_Genjet_parton1<<endl;
	     cout<< "ID1_minimo: "<< Id1_min<<endl;

	     }
	   }else{
	     if(deltaR_currentGenjet_parton2< deltaRmin_Genjet_parton2){
	       deltaRmin_Genjet_parton2 = deltaR_currentGenjet_parton2;
	       Id2_min = ijet;

	       cout<< "DeltaR2 minimo: "<< deltaRmin_Genjet_parton2<<endl;
	       cout<< "ID2_minimo: "<< Id2_min<<endl;

	     }
	   }
	 }   
       
       cout<< "DeltaR1 minimo: "<< deltaRmin_Genjet_parton1<<endl;
       cout<< "ID1_minimo: "<< Id1_min<<endl;
       cout<< "DeltaR2 minimo: "<< deltaRmin_Genjet_parton2<<endl;
       cout<< "ID2_minimo: "<< Id2_min<<endl;


       if(Id1_min < -1 || Id2_min < -1){

	 cout << "+++++++++++++++++++++++                         WARNING -1000                            ++++++++++++++++++++++++++++++++++"<<endl;
	 for(Long64_t ijet=0; ijet<no_Genjets_ak4; ijet++)
	   { 
	     
	     TLorentzVector currentGenJet;
	     currentGenJet.SetPtEtaPhiM( jetPtGenAK4->at(ijet), jetEtaGenAK4->at(ijet), jetPhiGenAK4->at(ijet), jetMassGenAK4->at(ijet));   
	     
	     deltaR_currentGenjet_parton1 = currentGenJet.DeltaR(parton1);
	     deltaR_currentGenjet_parton2 = currentGenJet.DeltaR(parton2);
	     
	     cout<< "deltaR_currentGenjet_parton1: "<< deltaR_currentGenjet_parton1 << endl;
	     cout<< "deltaR_currentGenjet_parton2: "<< deltaR_currentGenjet_parton2 << endl;
	     cout<< "DeltaR1 minimo: "<< deltaRmin_Genjet_parton1<<endl;
	     cout<< "ID1_minimo: "<< Id1_min<<endl;
	     cout<< "DeltaR2 minimo: "<< deltaRmin_Genjet_parton2<<endl;
	     cout<< "ID2_minimo: "<< Id2_min<<endl;
	   }
       }

*/
       
       /*
	   if(Id1_min > -1){
	     
	     CreateAndFillUserTH1D("H_DeltaR_Gen_hp1", 100, 0. , 10. , deltaR_Gen_hp1_min );
	     CreateAndFillUserTH1D("H_Id1_min_Gen", 16, -0.5 , 15.5 , Id1_min );
	  	     
	       if(deltaR_Gen_hp_min <= 0.3){
	       //Calcolo ratio con un taglio sensato in DeltaR
	       
	       double     R_PtGenjet_Ptjet1hp;
	       
	       R_PtGenjet_Ptjet1hp = jetPtGenAK4->at(Id_min) / Genjet1_hp.Pt() ;  	
	   
	   CreateAndFillUserTH1D("H_R_PtGenjet_Ptjet1hp", 100, 0. , 10. , R_PtGenjet_Ptjet1hp );
	   CreateAndFillUserTProfile("Profile_PtGen", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenjet_Ptjet1hp );  //, Double_t weight)	 
	   
	   if(	fabs( Genjet1_hp.Eta() ) >= 0. && fabs( Genjet1_hp.Eta() ) < 0.5){ 
	     CreateAndFillUserTProfile("Profile_PtGen_05", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 0.5 && fabs( Genjet1_hp.Eta() ) < 1.){ 
	     CreateAndFillUserTProfile("Profile_PtGen_1", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1. && fabs( Genjet1_hp.Eta() ) < 1.5){ 
	     CreateAndFillUserTProfile("Profile_PtGen_15", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1.5 && fabs( Genjet1_hp.Eta() ) < 2.){ 
	     CreateAndFillUserTProfile("Profile_PtGen_2", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 2. && fabs( Genjet1_hp.Eta() ) < 2.5){ 
	     CreateAndFillUserTProfile("Profile_PtGen_25", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   
	 }
	 		   
	   
       
	   if(Id2_min > -1){
	     
	     CreateAndFillUserTH1D("H_DeltaR_Gen_hp2", 100, 0. , 10. , deltaR_Gen_hp2_min );
	     CreateAndFillUserTH1D("H_Id2_min_Gen", 16, -0.5 , 15.5 , Id2_min );
	   }
       

       if( Genwj1.Pt()>0 && Genwj2.Pt()>0 ){
	 
	 double DeltaR_GenWide_hp = 1000 ;
	 double DeltaR1_GenWide_hp = Genwj1.DeltaR(parton1);
	 double DeltaR2_GenWide_hp = Genwj2.DeltaR(parton2);
	 
	 if(DeltaR1_GenWide_hp < DeltaR2_GenWide_hp){
	   DeltaR_GenWide_hp = DeltaR1_GenWide_hp;
	 }     else if(DeltaR1_GenWide_hp > DeltaR2_GenWide_hp){
	   DeltaR_GenWide_hp = DeltaR2_GenWide_hp;
	 }
	 
	 CreateAndFillUserTH1D("H_DeltaR_GenWide_hp", 100, 0. , 10. , DeltaR_GenWide_hp );
	 
	 if(DeltaR_GenWide_hp <= 0.3){
	   //Calcolo ratio con un taglio sensato in DeltaR
	   
	   double     R_PtGenWidejet_Ptjet1hp;
	   /*
	   if(DeltaR1_GenWide_hp <= DeltaR2_GenWide_hp){
	     R_PtGenWidejet_Ptjet1hp = Genwj1.Pt() / Genjet1_hp.Pt() ;
	   }else if(DeltaR1_GenWide_hp > DeltaR2_GenWide_hp){
	     R_PtGenWidejet_Ptjet1hp = Genwj2.Pt() / Genjet1_hp.Pt() ;
	   }
	   CreateAndFillUserTH1D("H_R_PtGenWidejet_Ptjet1hp", 100, 0. , 10. , R_PtGenWidejet_Ptjet1hp );
	   CreateAndFillUserTProfile("Profile_PtGenWide", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenWidejet_Ptjet1hp );  //, Double_t weight)	 
	   
	   if(	fabs( Genjet1_hp.Eta() ) >= 0. && fabs( Genjet1_hp.Eta() ) < 0.5){ 
	     CreateAndFillUserTProfile("Profile_PtGenWide_05", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 0.5 && fabs( Genjet1_hp.Eta() ) < 1.){ 
	     CreateAndFillUserTProfile("Profile_PtGenWide_1", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1. && fabs( Genjet1_hp.Eta() ) < 1.5){ 
	     CreateAndFillUserTProfile("Profile_PtGenWide_15", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1.5 && fabs( Genjet1_hp.Eta() ) < 2.){ 
	     CreateAndFillUserTProfile("Profile_PtGenWide_2", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 2. && fabs( Genjet1_hp.Eta() ) < 2.5){ 
	     CreateAndFillUserTProfile("Profile_PtGenWide_25", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtGenWidejet_Ptjet1hp );  //, Double_t weight)	 
	     }
	     }
     }
   }*/

     //++++++++++++++++++++++++++++++++++++++++++++

     //Livello reco
     //ricostruisce i widejet partendo di recojet_ak4
 
     size_t no_jets_ak4=jetPtAK4->size();


     TLorentzVector jet1, jet2, dijet;
     vector<TLorentzVector> jets;

     TLorentzVector wj1_tmp, wj2_tmp;
     TLorentzVector wj1, wj2, wdijet;
     vector<TLorentzVector> widejets;

     double wideJetDeltaR_ = getPreCutValue1("DeltaR");

     double MJJ = 0; 
     double DeltaEtaJJ = 0;
     double DeltaPhiJJ = 0;
     

     
     ////////////////////////////////////////////////////////////


     if(no_jets_ak4>=2)
       {
	 if(fabs(jetEtaAK4->at(0)) < getPreCutValue1("jetFidRegion") 
	    && jetPtAK4->at(0) > getPreCutValue1("pt0Cut"))
	   {
	     if(fabs(jetEtaAK4->at(1)) < getPreCutValue1("jetFidRegion") 
		&& jetPtAK4->at(1) > getPreCutValue1("pt1Cut"))
	       {
		 //  TLorentzVector jet1, jet2;
		 jet1.SetPtEtaPhiM(jetPtAK4->at(0),jetEtaAK4->at(0),jetPhiAK4->at(0),jetMassAK4->at(0));
		 jet2.SetPtEtaPhiM(jetPtAK4->at(1),jetEtaAK4->at(1),jetPhiAK4->at(1),jetMassAK4->at(1));

		 //costruisco la massa invariante di questi due jet (prima di fare le somme per il widejet)
		 
		 if( jet1.Pt()>0 && jet2.Pt()>0 )
		   {
		     // Create dijet(reco) system
		     dijet = jet1 + jet2;
		     MJJ = dijet.M();
		     DeltaEtaJJ = fabs(jet1.Eta()-jet2.Eta());
		     DeltaPhiJJ = fabs(jet1.DeltaPhi(jet2));

		     CreateAndFillUserTH1D("H_MJJ", 1000, 0. , 10000. , MJJ );
		     CreateAndFillUserTH1D("H_DeltaEtaJJ", 100, 0. , 8. , DeltaEtaJJ );

		     // Put widejets in the container
		     jets.push_back( jet1 );
		     jets.push_back( jet2 );
		   }
		 
		 /*
     double deltaR_Jet_hp_min = 1000;
     double deltaR_currentJet_hp = 1000;
     int IdJet_min = -1000;

     if(no_jets_ak4>=2){
       for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++)
	 { 
	   
	   if( fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
	       && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
	       && jetPtAK4->at(ijet) > getPreCutValue1("ptCut"))
	     {
	       
	       TLorentzVector currentJet;
	       currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   

	       deltaR_currentJet_hp = currentJet.DeltaR(parton1);

	       if( deltaR_currentJet_hp < deltaR_Jet_hp_min){
		 
		 deltaR_Jet_hp_min = deltaR_currentJet_hp;
		 IdJet_min = ijet;

	       }
	     }
	 }	   
     }
     
     
     if(IdJet_min > -1){
       
       CreateAndFillUserTH1D("H_DeltaR_Jet_hp", 100, 0. , 10. , deltaR_Jet_hp_min );
       CreateAndFillUserTH1D("H_Id_min_Jet", 16, -0.5 , 15.5 , IdJet_min );
       
       if(deltaR_Jet_hp_min <= 0.3){
	 //Calcolo ratio con un taglio sensato in DeltaR
	 
	 double     R_Ptjet_Ptjet1hp;
	 
	 R_Ptjet_Ptjet1hp = jetPtAK4->at(IdJet_min) / Genjet1_hp.Pt() ;  	
	 
	 CreateAndFillUserTH1D("H_R_Ptjet_Ptjet1hp", 100, 0. , 10. , R_Ptjet_Ptjet1hp );
	 
	 CreateAndFillUserTProfile("Profile_PtReco", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_Ptjet_Ptjet1hp );  //, Double_t weight)	 

	   if(	fabs( Genjet1_hp.Eta() ) >= 0. && fabs( Genjet1_hp.Eta() ) < 0.5){ 
	     CreateAndFillUserTProfile("Profile_PtReco_05", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_Ptjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 0.5 && fabs( Genjet1_hp.Eta() ) < 1.){ 
	     CreateAndFillUserTProfile("Profile_PtReco_1", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_Ptjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1. && fabs( Genjet1_hp.Eta() ) < 1.5){ 
	     CreateAndFillUserTProfile("Profile_PtReco_15", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_Ptjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1.5 && fabs( Genjet1_hp.Eta() ) < 2.){ 
	     CreateAndFillUserTProfile("Profile_PtReco_2", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_Ptjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 2. && fabs( Genjet1_hp.Eta() ) < 2.5){ 
	     CreateAndFillUserTProfile("Profile_PtReco_25", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_Ptjet_Ptjet1hp );  //, Double_t weight)	 
	   }
	 
       }		   
     }*/
	
     //++++++++++++++++++++++++++++++++++++++++

     //costruisco reco widejet 

		 for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4

		     TLorentzVector currentJet;
		     
		     if(fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
			&& idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
			&& jetPtAK4->at(ijet) > getPreCutValue1("ptCut")){

			 TLorentzVector currentJet;
			 currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   
			 
			 double DeltaR1 = currentJet.DeltaR(jet1);
			 double DeltaR2 = currentJet.DeltaR(jet2);
			 
			 if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
			     wj1_tmp += currentJet;
			 }
			 else if(DeltaR2 < wideJetDeltaR_){
			     wj2_tmp += currentJet;
			 }			 
		     } // if AK4 jet passes fid and jetid.
		 } //end of ak4 jet loop		     
		 
		 // if(wj1_tmp.Pt()==0 && wj2_tmp.Pt() ==0) 
		 // std::cout << " wj1_tmp.Pt() IN  " <<wj1_tmp.Pt()  << " wj2_tmp.Pt() " <<  wj2_tmp.Pt()  << std::endl;		     
   		 
	       } //fid, jet id, pt cut
	   } //fid, jet id, pt cut
       } // end of two jets.
     
     // Re-order the wide jets in pt
     if( wj1_tmp.Pt() > wj2_tmp.Pt())
       {
	 wj1 = wj1_tmp;
	 wj2 = wj2_tmp;
       }
     else
       {
	 wj1 = wj2_tmp;
	 wj2 = wj1_tmp;
       }

     double MJJWide = 0; 
     double DeltaEtaJJWide = 0;
     double DeltaPhiJJWide = 0;

     if( wj1.Pt()>0 && wj2.Pt()>0 )
     {

       // Create dijet system
       wdijet = wj1 + wj2;
       MJJWide = wdijet.M();
       DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
       DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));


       CreateAndFillUserTH1D("H_MJJWide", 1000, 0. , 10000. , MJJWide );
       CreateAndFillUserTH1D("H_DeltaEtaJJWide", 100, 0. , 8. , DeltaEtaJJWide );

       // Put widejets in the container
       widejets.push_back( wj1 );
       widejets.push_back( wj2 );
     }




/*       
       double DeltaR1_Wide_hp = wj1.DeltaR(parton1);
       double DeltaR2_Wide_hp = wj2.DeltaR(parton1);
       double DeltaR_Wide_hp = 0 ;
       
       if(DeltaR1_Wide_hp < DeltaR2_Wide_hp){
	 DeltaR_Wide_hp = DeltaR1_Wide_hp;
       }     else if(DeltaR1_Wide_hp > DeltaR2_Wide_hp){
	 DeltaR_Wide_hp = DeltaR2_Wide_hp;
       }
       
       CreateAndFillUserTH1D("H_DeltaR_Wide_hp", 100, 0. , 10. , DeltaR_Wide_hp );
       
       if(DeltaR_Wide_hp <= 0.3){
	 //Calcolo ratio con un taglio sensato in DeltaR
	 
	 double     R_PtWidejet_Ptjet1hp;
	 /*
	 if(DeltaR1_Wide_hp <= DeltaR2_Wide_hp){
	   R_PtWidejet_Ptjet1hp = wj1.Pt() / Genjet1_hp.Pt() ;
	 }else if(DeltaR1_Wide_hp > DeltaR2_Wide_hp){
	   R_PtWidejet_Ptjet1hp = wj2.Pt() / Genjet1_hp.Pt() ;
	 }
	 CreateAndFillUserTH1D("H_R_PtWidejet_Ptjet1hp", 100, 0. , 10. , R_PtWidejet_Ptjet1hp );
	 
	 CreateAndFillUserTProfile("Profile_PtRecoWide", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtWidejet_Ptjet1hp );  //, Double_t weight)	 

	   if(	fabs( Genjet1_hp.Eta() ) >= 0. && fabs( Genjet1_hp.Eta() ) < 0.5){ 
	     CreateAndFillUserTProfile("Profile_PtRecoWide_05", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 0.5 && fabs( Genjet1_hp.Eta() ) < 1.){ 
	     CreateAndFillUserTProfile("Profile_PtRecoWide_1", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1. && fabs( Genjet1_hp.Eta() ) < 1.5){ 
	     CreateAndFillUserTProfile("Profile_PtRecoWide_15", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 1.5 && fabs( Genjet1_hp.Eta() ) < 2.){ 
	     CreateAndFillUserTProfile("Profile_PtRecoWide_2", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	   if(	fabs( Genjet1_hp.Eta() ) >= 2. && fabs( Genjet1_hp.Eta() ) < 2.5){ 
	     CreateAndFillUserTProfile("Profile_PtRecoWide_25", 1000, 0. , 10000. , 0., 3., Genjet1_hp.Pt(), R_PtWidejet_Ptjet1hp );  //, Double_t weight)	 
	   }
	    
       }
*/
     
     //+++++++++++++++++++++++++++++++++++++++++++++++++++

     //== Fill Variables ==

     fillVariableWithValue("run",runNo);     
     fillVariableWithValue("event",evtNo);     
     fillVariableWithValue("lumi",lumi);     
     fillVariableWithValue("nVtx",nvtx);     
     fillVariableWithValue("nJet",widejets.size());

     // Trigger
     int NtriggerBits = triggerResult->size();
     if( NtriggerBits > 0)
       fillVariableWithValue("passHLT",triggerResult->at(0));// HLT_PFHT900_v*    

     if( no_jets_ak4 >=1 )
       fillVariableWithValue("IdTight_j1",idTAK4->at(0));

     if( no_jets_ak4 >=2 )
       fillVariableWithValue("IdTight_j2",idTAK4->at(1));

     if( widejets.size() >= 1 )
       {
         fillVariableWithValue( "pT_j1", widejets[0].Pt() );
         fillVariableWithValue( "eta_j1", widejets[0].Eta());

	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "phi_j1", widejets[0].Phi());
         fillVariableWithValue( "neutrHadEnFrac_j1", jetNhfAK4->at(0));
         fillVariableWithValue( "chargedHadEnFrac_j1", jetChfAK4->at(0));
         fillVariableWithValue( "photonEnFrac_j1", jetPhfAK4->at(0));
         fillVariableWithValue( "eleEnFract_j1", jetElfAK4->at(0));
         fillVariableWithValue( "muEnFract_j1", jetMufAK4->at(0));
       }

     if( widejets.size() >= 2 )
       {
         fillVariableWithValue( "pT_j2", widejets[1].Pt() );
         fillVariableWithValue( "eta_j2", widejets[1].Eta());
	 fillVariableWithValue( "GendeltaETAjj_hp", diparton_DeltaEtaJJ ) ;
         fillVariableWithValue( "Genmjj_hp", diparton_MJJ ) ;
	 fillVariableWithValue( "GendeltaETAjj", GenDeltaEtaJJ ) ;
         fillVariableWithValue( "Genmjj", GenMJJ ) ; 
	 fillVariableWithValue( "GendeltaETAjj_wide", GenDeltaEtaJJWide ) ;
         fillVariableWithValue( "Genmjj_wide", GenMJJWide ) ;
	 fillVariableWithValue( "deltaETAjj", DeltaEtaJJ ) ;
         fillVariableWithValue( "mjj", MJJ ) ; 
	 fillVariableWithValue( "deltaETAjj_wide", DeltaEtaJJWide ) ;
         fillVariableWithValue( "mjj_wide", MJJWide ) ;


	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "phi_j2", widejets[1].Phi());	
         fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(1));
         fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(1));
         fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(1));
         fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(1));
         fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(1));
	 fillVariableWithValue( "deltaPHIjj", DeltaPhiJJWide ) ;

       }

     //no cuts on these variables, just to store in output
     fillVariableWithValue("trueVtx",PileupInteractions->at(12));
     fillVariableWithValue("MET",met);
     double METoverHTAK4=double(met/htAK4);
     fillVariableWithValue("METoverHTAK4",METoverHTAK4);
     fillVariableWithValue("HTAK4",htAK4);
     fillVariableWithValue("ptHat",ptHat);

     // Evaluate cuts (but do not apply them)
     evaluateCuts();
     
     // optional call to fill a skim with the full content of the input roottuple
     //if( passedCut("nJetFinal") ) fillSkimTree();
     
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
