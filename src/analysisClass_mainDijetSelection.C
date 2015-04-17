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

int n_bin;
// GeV for pt binning
int step_pt=200;
char HistoName[200];
int xmin = 0. ;
int xmax = 2. ;
int histo_bin = 100;


double invariant_mass(double m1, double m2, double E1, double E2, double Angle12) {

  double p1 = sqrt( E1*E1 - m1*m1 );
  double p2 = sqrt( E2*E2 - m2*m2 );

  double invariant_mass =sqrt( m1*m1 + m2*m2 + 2*(E1*E2 - fabs(p1)*fabs(p2)*cos(Angle12) ) );
  return invariant_mass;
}

double invariant_mass_approssimation(double E1, double E2, double Angle12) {

  double invariant_mass =sqrt( 2*E1*E2 *(1- cos(Angle12) ) );
  return invariant_mass;
}

//////////////////////////////////////////////////////////////////////////////////////////
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

   for(int ii=0; ii<5; ii++){
     n_bin = 4500/ step_pt;
     for(int jj=0; jj<n_bin ; jj++){	 
       int pt_bin_max = (jj*step_pt)+step_pt;
       sprintf(HistoName,"Histo_R_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);	 
       CreateUserTH1D(HistoName, histo_bin, xmin, xmax);
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

     //divido in due il dataset -> per calcolare le funzioni di smearing uso i dispari (butto i pari)
     if( (int)jentry%2 == 0 ) continue;

     cout<<" "<<endl;
     cout<<"Evento numero "<<jentry<<endl;

     CreateAndFillUserTH1D("H_step_pt", 1000, 0. , 1000. , step_pt );

     ////////////////////////////////////////////////////////////////////////////

     // Partoni generati

     TLorentzVector p1, p2;
     double p1_pdgId;
     double p2_pdgId;
     double parton1_pdgId;
     double parton2_pdgId;
     TLorentzVector p0; //
     TLorentzVector parton1, parton2;
     TLorentzVector diparton;

     ///////////////////////////////////////////////

     bool Histogram_InitialObjects = false;

     if(Histogram_InitialObjects == true){
     //studio su partoni e risonanza a livello generatore
     double no_objects = gen_pt->size();
     cout<<"Numero di oggetti a livello generatore: "<<no_objects<<endl;     

     for(int nn=0; nn < no_objects; nn++){

       cout<<"Object "<<nn<<endl;
       cout<<"Pt: "<<gen_pt->at(nn)<< " Eta: "<<gen_eta->at(nn)<<" Phi: "<<gen_phi->at(nn)<<" M: "<<" pdgId: "<<gen_pdgId->at(nn)<<endl;
       
       CreateAndFillUserTH1D("H_Objects_size", 6, -0.5, 5.5, no_objects );
       CreateAndFillUserTH1D("H_Objects_pdgId", 50, -4, 4, gen_pdgId->at(nn) );
       CreateAndFillUserTH1D("H_Objects_Eta", 100, -10, 10, gen_eta->at(nn) );
       CreateAndFillUserTH1D("H_Objects_Pt", 100, 0, 2000, gen_pt->at(nn) );
       CreateAndFillUserTH1D("H_Objects_Phi", 50, -4, 4, gen_phi->at(nn) );
     }
     }
     ///////////////////////////////////////////////////////////////////////
     // Parton Level: definisco i partoni prodotti dalla risonanza
     if(gen_pt->size() == 3){
       
       CreateAndFillUserTH1D("H_partonsize_3", 5, 0.5, 5.5, gen_pt->size() );
       
       p1.SetPxPyPzE(gen_px->at(1) , gen_py->at(1) , gen_pz->at(1), gen_energy->at(1) );
       p1_pdgId = gen_pdgId -> at(1);
       
       p2.SetPxPyPzE(gen_px->at(2) , gen_py->at(2) , gen_pz->at(2), gen_energy->at(2) );
       p2_pdgId = gen_pdgId -> at(2);

       //dovrebbe essere la risonanza
       p0.SetPxPyPzE(gen_px->at(0) , gen_py->at(0) , gen_pz->at(0), gen_energy->at(0) );

       CreateAndFillUserTH1D("H_partonsize_3_Eta_resonance", 50, -10, 10, gen_eta->at(0) );
       CreateAndFillUserTH1D("H_partonsize_3_Eta", 50, -10, 10, p1.Eta() );
       CreateAndFillUserTH1D("H_partonsize_3_Eta", 50, -10, 10, p2.Eta() );
       CreateAndFillUserTH1D("H_partonsize_3_Eta_SumPartons",50,-10,10, (p1+p2).Eta() );
       CreateAndFillUserTH1D("H_partonsize_3_Pt_SumPartons",100,0,2000, (p1+p2).Pt() );
     }
     
     if(gen_pt->size() == 4){
       
       CreateAndFillUserTH1D("H_partonsize_4", 5, 0.5, 5.5, gen_pt->size() );
       
       p1.SetPxPyPzE(gen_px->at(2) , gen_py->at(2) , gen_pz->at(2), gen_energy->at(2) );
       p1_pdgId = gen_pdgId -> at(2);
       
       p2.SetPxPyPzE(gen_px->at(3) , gen_py->at(3) , gen_pz->at(3), gen_energy->at(3) );
       p2_pdgId = gen_pdgId -> at(3);     
     }
     
     if(gen_pt->size() == 5){
       
       CreateAndFillUserTH1D("H_partonsize_5", 5, 0.5, 5.5, gen_pt->size() );
       
       p1.SetPxPyPzE(gen_px->at(3) , gen_py->at(3) , gen_pz->at(3), gen_energy->at(3) );
       p1_pdgId = gen_pdgId -> at(3);
       
       p2.SetPxPyPzE(gen_px->at(4) , gen_py->at(4) , gen_pz->at(4), gen_energy->at(4) );
       p2_pdgId = gen_pdgId -> at(4);   

       CreateAndFillUserTH1D("H_partonsize_5_Eta_final", 50, -10, 10, p1.Eta() );
       CreateAndFillUserTH1D("H_partonsize_5_Eta_final", 50, -10, 10, p2.Eta() );
       CreateAndFillUserTH1D("H_partonsize_5_Eta_initial", 50, -10, 10, gen_eta->at(0) );
       CreateAndFillUserTH1D("H_partonsize_5_Eta_initial", 50, -10, 10, gen_eta->at(1) );
       CreateAndFillUserTH1D("H_partonsize_5_Eta_resonance", 50, -10, 10, gen_eta->at(2) );
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
     
     cout<<"PARTONE 1"<<endl;
     cout<<"Pt: "<<parton1.Pt()<< " Eta: "<<parton1.Eta()<<" Phi: "<<parton1.Phi()<<" M: "<<parton1.M()<<" pdgId: "<<parton1_pdgId<<endl;
     cout<<"PARTONE 2"<<endl;
     cout<<"Pt: "<<parton2.Pt()<< " Eta: "<<parton2.Eta()<<" Phi: "<<parton2.Phi()<<" M: "<<parton2.M()<<" pdgId: "<<parton2_pdgId<<endl;

     //calcolo le funzioni di smearing con i genwidejets ora -> rimosso taglio in |eta| dei partoni     
     //     if( fabs(parton1.Eta() ) > 2.5 || fabs(parton2.Eta() ) > 2.5) cout<<"Uno dei due partoni e' fuori range in eta"<<endl;
     //     if( fabs(parton1.Eta() ) > 2.5 || fabs(parton2.Eta() ) > 2.5) continue; //see
     
     double diparton_DeltaEtaJJ;
     double diparton_cosThetaStar;
     
     // Create dijet(reco) system
     diparton = parton1 + parton2;
     diparton_DeltaEtaJJ = fabs( parton1.Eta() - parton2.Eta() );
     diparton_cosThetaStar = tanh ( diparton_DeltaEtaJJ/2 ) ;
     
     //////////////////////////////////////////////////////////////
     // Test: seleziono solo quark down(1) oppure up(2)
     //    if(fabs(parton1_pdgId ==2) || fabs(parton2_pdgId==2)){
     //cout<<parton1_pdgId<<endl;
     //cout<<parton2_pdgId<<endl;
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++ 
     //Livello gen
     //ricostruisce i Gen_widejet partendo dai genjet_ak4
     
     size_t no_Genjets_ak4 = jetPtGenAK4->size();
     TLorentzVector Genjet1, Genjet2;
     TLorentzVector Genwj1_tmp, Genwj2_tmp;
     TLorentzVector Genwj1, Genwj2;
     TLorentzVector Genwdijet;
     double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");
 
     TLorentzVector Genwdijet_11;
     TLorentzVector Genwj1_16, Genwj2_16;
     TLorentzVector Genwdijet_16;
     TLorentzVector Genwj1_21, Genwj2_21;
     TLorentzVector Genwdijet_21;
     TLorentzVector Genwj1_26, Genwj2_26;
     TLorentzVector Genwdijet_26;
     TLorentzVector Genwj1_31, Genwj2_31;
     TLorentzVector Genwdijet_31;
    
     if(no_Genjets_ak4 < 2) cout<<"# Genjet AK4 <2. Butto evento"<<endl;
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
       
       // cout<<"Genjet numero "<<ijet<<endl;
       // cout<<"Pt: "<<currentGenJet.Pt()<<"          Eta: "<<currentGenJet.Eta()<<"         Phi: "<<currentGenJet.Phi()<<"         M: "<<currentGenJet.M()<<endl;
       
       double DeltaR1 = currentGenJet.DeltaR(Genjet1);
       double DeltaR2 = currentGenJet.DeltaR(Genjet2);
       // cout<<"DeltaR tra il genjet numero "<<ijet<<" e il genjet 0: "<<DeltaR1<<endl;
       // cout<<"DeltaR tra il genjet numero "<<ijet<<" e il genjet 1: "<<DeltaR2<<endl;
       
       if(DeltaR1 < DeltaR2 && DeltaR1 < GenwideJetDeltaR_){
	 // cout<<"Sommo il jet numero "<<ijet<<" al genjet numero 0"<<endl;
	 Genwj1_tmp += currentGenJet;
       }
       else if(DeltaR2 < GenwideJetDeltaR_){
	 // cout<<"Sommo il jet numero "<<ijet<<" al genjet numero 1"<<endl;
	 Genwj2_tmp += currentGenJet;
       }
       //++++++++++++++++++++++++++++++++++++++++++++++++
       // Test sulla costruzine dei genwideJet   
       //R = 1.6
       if(DeltaR1 < DeltaR2 && DeltaR1 < 1.6){
	 Genwj1_16 += currentGenJet;
       }
       else if(DeltaR2 < 1.6){
	 Genwj2_16 += currentGenJet;
       }
       //R = 2.1
       if(DeltaR1 < DeltaR2 && DeltaR1 < 2.1){
	 Genwj1_21 += currentGenJet;
       }
       else if(DeltaR2 < 2.1){
	 Genwj2_21 += currentGenJet;
       }
       //R = 2.6
       if(DeltaR1 < DeltaR2 && DeltaR1 < 2.6){
	 Genwj1_26 += currentGenJet;
       }
       else if(DeltaR2 < 2.6){
	 Genwj2_26 += currentGenJet;
       }
       //R = 3.1
       if(DeltaR1 < DeltaR2 && DeltaR1 < 3.1){
	 Genwj1_31 += currentGenJet;
       }
       else if(DeltaR2 < 3.1){
	 Genwj2_31 += currentGenJet;
       }
       //++++++++++++++++++++++++++++
     }//end of ak4 loop
     
     if( Genwj1_tmp.Pt() <0 || Genwj2_tmp.Pt() <0) cout<<"Pt GenWidejet <0"<<endl;
     if( Genwj1_tmp.Pt() <0 || Genwj2_tmp.Pt() <0) continue;
     
     Genwdijet_11 = Genwj1_tmp + Genwj2_tmp;
     Genwdijet_16 = Genwj1_16 + Genwj2_16;
     Genwdijet_21 = Genwj1_21 + Genwj2_21;
     Genwdijet_26 = Genwj1_26 + Genwj2_26;
     Genwdijet_31 = Genwj1_31 + Genwj2_31;

     CreateAndFillUserTH1D("H_GendijetWide_M_R11", 1000, 0. , 10000. , Genwdijet_11.M() );
     CreateAndFillUserTH1D("H_GendijetWide_M_R16", 1000, 0. , 10000. , Genwdijet_16.M() );
     CreateAndFillUserTH1D("H_GendijetWide_M_R21", 1000, 0. , 10000. , Genwdijet_21.M() );
     CreateAndFillUserTH1D("H_GendijetWide_M_R26", 1000, 0. , 10000. , Genwdijet_26.M() );
     CreateAndFillUserTH1D("H_GendijetWide_M_R31", 1000, 0. , 10000. , Genwdijet_31.M() );
     
     //////////////////////////////////////////////////////
     
     //associo i Gen widejets_tmp ai partoni
     double DeltaR_GenWidejet1_tmp_parton1 = Genwj1_tmp.DeltaR(parton1);     
     double DeltaR_GenWidejet2_tmp_parton1 = Genwj2_tmp.DeltaR(parton1);     
     
     if(DeltaR_GenWidejet1_tmp_parton1 < DeltaR_GenWidejet2_tmp_parton1){
       // cout<<"GenWidejet_tmp 1 associato al partone 1 Allora Genwj1 = Genwj1_tmp"<<endl;
       // cout<<"GenWidejet_tmp 2 associato al partone 2 Allora Genwj2 = Genwj2_tmp"<<endl;
       Genwj1 = Genwj1_tmp;
       Genwj2 = Genwj2_tmp;
     }
     else{
       // cout<<"GenWidejet_tmp 2 associato al partone 1 Allora Genwj1 = Genwj2_tmp"<<endl;
       // cout<<"GenWidejet_tmp 1 associato al partone 2 Allora Genwj2 = Genwj1_tmp"<<endl;
       Genwj1 = Genwj2_tmp; 
       Genwj2 = Genwj1_tmp; 
     }

     double DeltaR_GenWideJet1_parton1;   
     double DeltaR_GenWideJet2_parton2;   
     DeltaR_GenWideJet1_parton1 = Genwj1.DeltaR(parton1);
     DeltaR_GenWideJet2_parton2 = Genwj2.DeltaR(parton2);
     cout<< "DeltaR tra Genwidejet 1 e partone 1 = "<<DeltaR_GenWideJet1_parton1<<endl;
     cout<< "DeltaR tra Genwidejet 2 e partone 2 = "<<DeltaR_GenWideJet2_parton2<<endl;
     
    CreateAndFillUserTH1D("H_DeltaR_GenWideJet1_Parton1", 100, 0. , 3. , DeltaR_GenWideJet1_parton1);
    CreateAndFillUserTH1D("H_DeltaR_GenWideJet2_Parton2", 100, 0. , 3. , DeltaR_GenWideJet2_parton2);     

     ////////////////////////////////////////////////////////////////////

     cout<<"GenWideJet 1"<<endl;
     cout<< "Pt: "<<Genwj1.Pt()<<"         Eta: "<< Genwj1.Eta()<<"         Phi: "<<Genwj1.Phi()<<"         M: "<<Genwj1.M()<<"          E: "<<Genwj1.E()<<endl;
     cout<<"GenWideJet 2"<<endl;
     cout<< "Pt: "<<Genwj2.Pt()<<"         Eta: "<< Genwj2.Eta()<<"         Phi: "<<Genwj2.Phi()<<"         M: "<<Genwj2.M()<<"          E: "<<Genwj2.E()<<endl;

     // Funzioni di smearing con i gen widejets nel range |eta|<2.5     
     if( fabs(Genwj1.Eta() ) > 2.5 || fabs(Genwj2.Eta() ) > 2.5) cout<<"Uno dei due GenWideJet e' fuori range in eta"<<endl;
     if( fabs(Genwj1.Eta() ) > 2.5 || fabs(Genwj2.Eta() ) > 2.5) continue; //see

     ///////////////////////////////////////////////////////

     // Create dijet system
     Genwdijet = Genwj1 + Genwj2;
     double GenDeltaEtaJJWide;
     double GenWide_cosThetaStar;    
     GenDeltaEtaJJWide = fabs(Genwj1.Eta()-Genwj2.Eta());
     GenWide_cosThetaStar = tanh( GenDeltaEtaJJWide/2 );
     

     if(GenDeltaEtaJJWide<=1.3){//cut in delta eta solo per questo plot
       CreateAndFillUserTH1D("H_GenDeltaEtaJJWide_Cut13", 100, 0. , 2. , GenDeltaEtaJJWide );
     }
     
     /////////////////////////////////////////////////////////
     //Livello reco
     //ricostruisce i widejet partendo dai recojet_ak4
     
     size_t no_jets_ak4=jetPtAK4->size();    
     TLorentzVector jet1, jet2;
     TLorentzVector wj1_tmp, wj2_tmp;
     TLorentzVector wj1, wj2; 
     TLorentzVector wdijet;
     double wideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     if(no_jets_ak4<2) cout<<"# recojet AK4 <2. Butto evento"<<endl;
     if(no_jets_ak4<2) continue;
     
     //cout<<"Costruisco Widejet"<<endl;    
     //cout<<"Ciclo sui recojet AK4"<<endl;
     
     //    TLorentzVector jet1, jet2;
     jet1.SetPtEtaPhiM(jetPtAK4->at(0), jetEtaAK4->at(0), jetPhiAK4->at(0), jetMassAK4->at(0));
     jet2.SetPtEtaPhiM(jetPtAK4->at(1), jetEtaAK4->at(1), jetPhiAK4->at(1), jetMassAK4->at(1));
     
     //costruisco reco widejet attorno ai jet a pT piu alto
     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
       
       TLorentzVector currentJet;
       currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet), jetEtaAK4->at(ijet), jetPhiAK4->at(ijet), jetMassAK4->at(ijet));   
       
       // cout<<"Jet numero "<<ijet<<endl;
       // cout<<"Pt: "<<currentJet.Pt()<<"         Eta: "<<currentJet.Eta()<<"          Phi: "<<currentJet.Phi()<<"         M: "<<currentJet.M()<<endl;
       
       double DeltaR1 = currentJet.DeltaR(jet1);
       double DeltaR2 = currentJet.DeltaR(jet2);
       // cout<<"DeltaR tra il jet numero "<<ijet<<" e il jet 0 : "<<DeltaR1<<endl;
       // cout<<"DeltaR tra il jet numero "<<ijet<<" e il jet 1 : "<<DeltaR2<<endl;       
       
       if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
	 // cout<<"Sommo il jet numero "<<ijet<<" al jet numero 0"<<endl;
	 wj1_tmp += currentJet;
       }
       else if(DeltaR2 < wideJetDeltaR_){
	 // cout<<"Sommo il jet numero "<<ijet<<" al jet numero 1"<<endl;
	 wj2_tmp += currentJet;
       }			 
     } //end of ak4 jet loop		     
     
     if( wj1_tmp.Pt()<0 || wj2_tmp.Pt()<0 ) cout<<"Pt Widejet <0"<<endl;;
     if( wj1_tmp.Pt()<0 || wj2_tmp.Pt()<0 ) continue;
     
     /////////////////////////////////////////////////////////////////////////
    
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
 
     //check the spatial resolution - Widejet-parton
     double Eta_difference_Widejet1_parton1 = wj1.Eta() - parton1.Eta();
     double Phi_difference_Widejet1_parton1 = wj1.Phi() - parton1.Phi() ;
     double DeltaR_WideJet1_parton1            = wj1.DeltaR(parton1);
     double Eta_difference_Widejet2_parton2 = wj2.Eta() - parton2.Eta() ;
     double Phi_difference_Widejet2_parton2 = wj2.Phi() - parton2.Phi() ;
     double DeltaR_WideJet2_parton2            = wj2.DeltaR(parton2);

     cout<< "DeltaR tra widejet 1 e partone 1 = "<<DeltaR_WideJet1_parton1<<endl;
     cout<< "DeltaR tra widejet 2 e partone 2 = "<<DeltaR_WideJet2_parton2<<endl;
     
     //calcolo smearing function with parton -> adesso lo faccio con i GenWideJets
     //    double R_PtWidejet_Ptparton[2];
     //R_PtWidejet_Ptparton[0] = wj1.Pt() / parton1.Pt() ;
     //R_PtWidejet_Ptparton[1] = wj2.Pt() / parton2.Pt() ;
     
    CreateAndFillUserTH1D("H_Eta_difference_Widejet1_parton1", 100, -2. , 2. , Eta_difference_Widejet1_parton1);
    CreateAndFillUserTH1D("H_Phi_difference_Widejet1_parton1", 100, -2. , 2. , Phi_difference_Widejet1_parton1);
    CreateAndFillUserTH1D("H_Eta_difference_Widejet2_parton2", 100, -2. , 2. , Eta_difference_Widejet2_parton2);
    CreateAndFillUserTH1D("H_Phi_difference_Widejet2_parton2", 100, -2. , 2. , Phi_difference_Widejet2_parton2);       
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_parton1", 100, 0. , 3. , DeltaR_WideJet1_parton1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_parton2", 100, 0. , 3. , DeltaR_WideJet1_parton1);
    
    /// F.
    double Angle_difference_Widejet1_parton1 = wj1.Angle(parton1.Vect());
    double Angle_difference_Widejet2_parton2= wj2.Angle(parton2.Vect());
    double Angle_difference_Partons= parton1.Angle(parton2.Vect());
    CreateAndFillUserTH1D("H_Angle_difference_WideJet1_parton1", 100, -3.15 , 3.15, Angle_difference_Widejet1_parton1 );
    CreateAndFillUserTH1D("H_Angle_difference_WideJet2_parton2", 100, -3.15 , 3.15, Angle_difference_Widejet2_parton2 );
    CreateAndFillUserTH1D("H_Angle_difference_Partons", 100, -3.15 , 3.15, Angle_difference_Partons );
    ///    
    //////////////////////////////////////////////////////////    
    //    cout<<" "<<endl;
    cout<<"WideJet 1"<<endl;
    cout<< "Pt: "<<wj1.Pt()<<"         Eta: "<< wj1.Eta()<<"         Phi: "<<wj1.Phi()<<"         M: "<<wj1.M()<<"          E: "<<wj1.E()<<endl;
    cout<<"WideJet 2"<<endl;
    cout<< "Pt: "<<wj2.Pt()<<"         Eta: "<< wj2.Eta()<<"         Phi: "<<wj2.Phi()<<"         M: "<<wj2.M()<<"          E: "<<wj2.E()<<endl;
    
    //////////////////////////////////////////////////////////////////////////////
    
    // Create dijet system
    wdijet = wj1 + wj2;
    double DeltaEtaJJWide = 1000;
    double RecoWide_cosThetaStar = 1000;    
    DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
    RecoWide_cosThetaStar = tanh( DeltaEtaJJWide/2 );

    if(DeltaEtaJJWide<=1.3){//cut in delta eta solo per questo plot
      CreateAndFillUserTH1D("H_DeltaEtaJJWide_Cut13", 100, 0. , 2. , DeltaEtaJJWide );
    }
    
    ////////////////////////////////////////////////////////////////////////////////////
    
    //Test sulle masse invarianti con la formula invece che con i TLorentzVector
    // A questo punto ho i due partoni parton1 parton2 
    //    i due widejets gen Genwj1 Genwj2
    // e i due widejets reco wj1 wj2
    
    double Angle_partons          = parton2.Angle(parton1.Vect());
    double Angle_Genwidejets  = Genwj2.Angle(Genwj1.Vect());
    double Angle_widejets        = wj2.Angle(wj1.Vect());
    
    CreateAndFillUserTH1D("H_Angle_Partons", 50, 0. , 3.15, Angle_partons );
    CreateAndFillUserTH1D("H_Angle_GenWideJets", 50, 0., 3.15, Angle_Genwidejets );
    CreateAndFillUserTH1D("H_Angle_WideJets", 50, 0., 3.15 , Angle_widejets );

    //caso 1 -> Tutto partoni
    double invariant_mass_AllGenerator = invariant_mass(parton1.M(), parton2.M(), parton1.E(), parton2.E(), Angle_partons);
    CreateAndFillUserTH1D("H_invariant_mass_AllGenerator", 1000, 0. , 10000. , invariant_mass_AllGenerator );
    
    //caso 2 ->  Direction Reco, Energy Partons
    double invariant_mass_PartonEnergy_RecoDirection = invariant_mass(parton1.M(), parton2.M(), parton1.E(), parton2.E(), Angle_widejets);
    CreateAndFillUserTH1D("H_invariant_mass_PartonEnergy_RecoDirection", 1000, 0. , 10000. , invariant_mass_PartonEnergy_RecoDirection );
    
    //caso 3 -> Energy Reco, Direction Partons
    double invariant_mass_PartonDirection_RecoEnergy = invariant_mass(wj1.M(), wj2.M(), wj1.E(), wj2.E(), Angle_partons);
    CreateAndFillUserTH1D("H_invariant_mass_PartonDirection_RecoEnergy", 1000, 0. , 10000. , invariant_mass_PartonDirection_RecoEnergy );
    
    //caso 4 -> Tutto GenWideJets
    double invariant_mass_GenWideJets = invariant_mass(Genwj1.M(), Genwj2.M(), Genwj1.E(), Genwj2.E(), Angle_Genwidejets);
    CreateAndFillUserTH1D("H_invariant_mass_GenWideJets", 1000, 0. , 10000. , invariant_mass_GenWideJets );

    //caso 5 -> Energy GenWideJets Direction RecoWideJets
    double invariant_mass_GenEnergy_RecoDirection = invariant_mass(Genwj1.M(), Genwj2.M(), Genwj1.E(), Genwj2.E(), Angle_widejets);
    CreateAndFillUserTH1D("H_invariant_mass_GenEnergy_RecoDirection", 1000, 0. , 10000. , invariant_mass_GenEnergy_RecoDirection );
    
    //caso 6 -> Direction GenWideJets Energy RecoWideJets
    double invariant_mass_GenDirection_RecoEnergy = invariant_mass(wj1.M(), wj2.M(), wj1.E(), wj2.E(), Angle_Genwidejets);
    CreateAndFillUserTH1D("H_invariant_mass_GenDirection_RecoEnergy", 1000, 0. , 10000. , invariant_mass_GenDirection_RecoEnergy );
    
    //caso 7 -> Tutto reco
    double invariant_mass_AllReco = invariant_mass(wj1.M(), wj2.M(), wj1.E(), wj2.E(), Angle_widejets);
    CreateAndFillUserTH1D("H_invariant_mass_AllReco", 1000, 0. , 10000. , invariant_mass_AllReco );
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Calcolo le masse invarianti di partoni, GenWidejets e RecoWideJets con le due formule (completa e approssimata)
    
    // Formula completa -> Partoni
    double invariant_mass_partons_Complete = invariant_mass(parton1.M(), parton2.M(), parton1.E(), parton2.E(), Angle_partons);
    CreateAndFillUserTH1D("H_invariant_mass_partons_Complete", 1000, 0. , 10000. , invariant_mass_partons_Complete );
    
    // Formula approssimata -> Partoni
    double invariant_mass_partons_Approx = invariant_mass_approssimation( parton1.E(), parton2.E(), Angle_partons);
    CreateAndFillUserTH1D("H_invariant_mass_partons_Approx", 1000, 0. , 10000. , invariant_mass_partons_Approx );
    
    // Formula completa -> GenWideJets
    double invariant_mass_GenWideJets_Complete = invariant_mass(Genwj1.M(), Genwj2.M(), Genwj1.E(), Genwj2.E(), Angle_Genwidejets);
    CreateAndFillUserTH1D("H_invariant_mass_GenWideJets_Complete", 1000, 0. , 10000. , invariant_mass_GenWideJets_Complete );
    
    // Formula approssimata -> GenWideJets
    double invariant_mass_GenWideJets_Approx = invariant_mass_approssimation( Genwj1.E(), Genwj2.E(), Angle_Genwidejets);
    CreateAndFillUserTH1D("H_invariant_mass_GenWideJets_Approx", 1000, 0. , 10000. , invariant_mass_GenWideJets_Approx );
    
    // Formula completa -> Reco WideJets
    double invariant_mass_WideJets_Complete = invariant_mass(wj1.M(), wj2.M(), wj1.E(), wj2.E(), Angle_widejets);
    CreateAndFillUserTH1D("H_invariant_mass_WideJets_Complete", 1000, 0. , 10000. , invariant_mass_WideJets_Complete );
    
    // Formula approssimata -> WideJets
    double invariant_mass_WideJets_Approx = invariant_mass_approssimation( wj1.E(), wj2.E(), Angle_widejets);
    CreateAndFillUserTH1D("H_invariant_mass_WideJets_Approx", 1000, 0. , 10000. , invariant_mass_WideJets_Approx );

    /////////////////////////////////////////////////////////////////////////////////////

    //Questo mi serve per calcolare le funzioni di smearing 
    double DeltaR_WideJet1_GenWideJet1 = wj1.DeltaR(Genwj1);     
    double DeltaR_WideJet2_GenWideJet2 = wj2.DeltaR(Genwj2);         
    cout<< "DeltaR tra widejet 1 e Genwidejet 1 = "<<DeltaR_WideJet1_GenWideJet1<<endl;
    cout<< "DeltaR tra widejet 2 e Genwidejet 2 = "<<DeltaR_WideJet2_GenWideJet2<<endl;
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_GenWideJet1", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_GenWideJet2", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
    
    //calcolo smearing function with GenWideJet
    double R_PtWidejet_PtGenWideJet[2];
    R_PtWidejet_PtGenWideJet[0] = wj1.Pt() / Genwj1.Pt() ;
    R_PtWidejet_PtGenWideJet[1] = wj2.Pt() / Genwj2.Pt() ;
    CreateAndFillUserTH1D("H_R_PtWidejet1_PtGenWideJet1", 100, 0. , 2. , R_PtWidejet_PtGenWideJet[0] );
    CreateAndFillUserTH1D("H_R_PtWidejet2_PtGenWideJet2", 100, 0. , 2. , R_PtWidejet_PtGenWideJet[1] );
   
    //voglio fare un ciclo sui GenWideJets -> costruisco array 
    double GenWideJet_Pt[2];
    double GenWideJet_Eta[2];    
    GenWideJet_Pt[0]   = Genwj1.Pt();
    GenWideJet_Eta[0] = Genwj1.Eta();
    
    GenWideJet_Pt[1]   = Genwj2.Pt();
    GenWideJet_Eta[1] = Genwj2.Eta();    
    
    if(DeltaR_WideJet1_GenWideJet1 > 0.3 || DeltaR_WideJet2_GenWideJet2 > 0.3) cout<<"Uno dei due widejet non matchea il genwidejet . Butto evento"<<endl;
    if(DeltaR_WideJet1_GenWideJet1 > 0.3 || DeltaR_WideJet2_GenWideJet2 > 0.3) continue;
    
    CreateAndFillUserTH1D("H_DeltaR_WideJet1_GenWideJet1_Match", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
    CreateAndFillUserTH1D("H_DeltaR_WideJet2_GenWideJet2_Match", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);

    //se commentato -> sto usando i jet in ogni bin     
    // see // considero i GenWideJets in un solo bin di |Eta| e pT
    // if( fabs(GenWideJet_Eta[0] ) >=0.5 || fabs( GenWideJet_Eta[1] ) >= 0.5) continue;
    // if( GenWideJet_Pt[0] <1400 || GenWideJet_Pt[0] >=1600) continue;
    // if( GenWideJet_Pt[1] <1400 || GenWideJet_Pt[1] >=1600) continue;

    cout<<"CALCOLO LE FUNZIONI DI SMEARING"<<endl;
    
    for(int kk=0; kk<2; kk++){//ciclo sui GenWideJets
      
      if(kk == 0) CreateAndFillUserTH1D("H_DeltaR_Totale", 100, 0. , 3. , DeltaR_WideJet1_GenWideJet1);
      if(kk == 1) CreateAndFillUserTH1D("H_DeltaR_Totale", 100, 0. , 3. , DeltaR_WideJet2_GenWideJet2);
      CreateAndFillUserTH1D("H_R_PtWidejet_PtGenWideJet", 100, 0. , 2. , R_PtWidejet_PtGenWideJet[kk] );
      
      cout<<"GenWideJet "<< kk+1<<endl;
      cout<<"Pt: "<<GenWideJet_Pt[kk]<<"         |Eta|: "<< fabs(GenWideJet_Eta[kk])<<endl;
      
      for(int ii=0; ii<5; ii++){	  
	double eta_bin_min = ii/2.;       
	double eta_bin_max = ii/2. +0.5;
	
	if( fabs(GenWideJet_Eta[kk]) >=eta_bin_min && fabs(GenWideJet_Eta[kk]) < eta_bin_max){	    
	  
	  cout<<"Passed bin in eta"<<endl;	
	  cout << "eta_bin: "<< eta_bin_min <<" - "<<eta_bin_max << endl;
	
	  for(int jj=0; jj<n_bin ; jj++){	 	      
	    int pt_bin_min = ( (jj-1)*step_pt)+step_pt;
	    int pt_bin_max = (jj*step_pt)+step_pt;
	
	    if(GenWideJet_Pt[kk] >=pt_bin_min && GenWideJet_Pt[kk] < pt_bin_max){ 
	      
	      cout<<"Passed bin in pt"<<endl;
	      cout << "pt_bin: "<< pt_bin_min <<" - "<<pt_bin_max<< endl;
	      
	      sprintf(HistoName,"Histo_R_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	      
	      CreateAndFillUserTH1D(HistoName, histo_bin, xmin, xmax, R_PtWidejet_PtGenWideJet[kk] );
	      cout<<"Riempio TH1D "<<HistoName<<" con R: "<<R_PtWidejet_PtGenWideJet[kk]<<endl;
	      
	    }//if parton_Pt
	  }//ciclo sui bin in pT
	}//if parton_Eta
      }//ciclo sui bin di eta
    }//ciclo sui partoni

    ///////////////////////////////////////////////////////////////////////////////////////////
    //Riempio tutti gli istogrammi qua

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
     CreateAndFillUserTH1D("H_diparton_DeltaEtaJJ", 100, 0. , 8. , diparton_DeltaEtaJJ );
     CreateAndFillUserTH1D("H_diparton_cosThetaStar", 50, 0. , 1. , diparton_cosThetaStar);
     
     //+++++++++++++++++++++++++++++++++++++++++++++
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
     CreateAndFillUserTH1D("H_GenDeltaEtaJJWide", 100, 0. , 8. , GenDeltaEtaJJWide );
     CreateAndFillUserTH1D("H_GendijetWide_cosThetaStar", 50, 0. , 1. , GenWide_cosThetaStar );
   
     //+++++++++++++++++++++++++++++++++++++++++++++
     // Reco WideJet
    CreateAndFillUserTH1D("H_WideJet1_Pt", 1000, 0. , 10000. , wj1.Pt()); 
    CreateAndFillUserTH1D("H_WideJet1_Eta", 100, -4. , 4. , wj1.Eta()); 
    CreateAndFillUserTH1D("H_WideJet1_Phi", 100, -4. , 4. , wj1.Phi());
    CreateAndFillUserTH1D("H_WideJet1_M", 1000, 0. , 10000. , wj1.M());  
    
    CreateAndFillUserTH1D("H_WideJet2_Pt", 1000, 0. , 10000. , wj2.Pt()); 
    CreateAndFillUserTH1D("H_WideJet2_Eta", 100, -4. , 4. , wj2.Eta()); 
    CreateAndFillUserTH1D("H_WideJet2_Phi", 100, -4. , 4. , wj2.Phi());
    CreateAndFillUserTH1D("H_WideJet2_M", 1000, 0. , 10000. , wj2.M());  
    
    //Variabili del sistema Reco Dijet Wide	   
    CreateAndFillUserTH1D("H_dijetWide_Pt", 1000, 0. , 10000. , wdijet.Pt() );
    CreateAndFillUserTH1D("H_dijetWide_Eta", 100, -4. , 4. , wdijet.Eta() );
    CreateAndFillUserTH1D("H_dijetWide_Phi", 100, -4. , 4. , wdijet.Phi() );
    CreateAndFillUserTH1D("H_dijetWide_M", 1000, 0. , 10000. , wdijet.M() );
    CreateAndFillUserTH1D("H_DeltaEtaJJWide", 100, 0. , 8. , DeltaEtaJJWide );
    CreateAndFillUserTH1D("H_dijetWide_cosThetaStar", 50, 0. , 1. , RecoWide_cosThetaStar );
    
    //////////////////////////////////////////////////////////////////////////////////////

    //== Fill Variables ==
    
    fillVariableWithValue("run",runNo);     
    fillVariableWithValue("event",evtNo);     
    fillVariableWithValue("lumi",lumi);     
    fillVariableWithValue("nVtx",nvtx);     
    //     fillVariableWithValue("nJet",widejets.size());
    
    // Trigger
    int NtriggerBits = triggerResult->size();
    if( NtriggerBits > 0)
      fillVariableWithValue("passHLT",triggerResult->at(0));// HLT_PFHT900_v*    
    
    if( no_jets_ak4 >=1 )
      fillVariableWithValue("IdTight_j1",idTAK4->at(0));
    
    if( no_jets_ak4 >=2 )
       fillVariableWithValue("IdTight_j2",idTAK4->at(1));

     //     if( widejets.size() >= 1 )
     //{
         fillVariableWithValue( "pT_j1", wj1.Pt() );
         fillVariableWithValue( "eta_j1", wj1.Eta());

	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "phi_j1", wj1.Phi());
         fillVariableWithValue( "neutrHadEnFrac_j1", jetNhfAK4->at(0));
         fillVariableWithValue( "chargedHadEnFrac_j1", jetChfAK4->at(0));
         fillVariableWithValue( "photonEnFrac_j1", jetPhfAK4->at(0));
         fillVariableWithValue( "eleEnFract_j1", jetElfAK4->at(0));
         fillVariableWithValue( "muEnFract_j1", jetMufAK4->at(0));
	 //}

	 //     if( widejets.size() >= 2 )
	 //{
         fillVariableWithValue( "pT_j2", wj2.Pt() );
         fillVariableWithValue( "eta_j2", wj2.Eta());
	 fillVariableWithValue( "deltaETAjj", DeltaEtaJJWide ) ;
         fillVariableWithValue( "mjj", wdijet.M() ) ;

	 //no cuts on these variables, just to store in output
         fillVariableWithValue( "phi_j2", wj2.Phi());	
         fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(1));
         fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(1));
         fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(1));
         fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(1));
         fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(1));
	 //	 fillVariableWithValue( "deltaPHIjj", DeltaPhiJJWide ) ;

	 //       }

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

////////////////////////////////////////////////////////////////////////

// rapporto con i genjet AK4/AK8 non servono piu perche abbiamo visto che hanno una risoluzione peggiore rispetto ai widejet
// sia generatore che reco non servono

    /*    
    double deltaR_currentGenjet_parton1 = 1000;
    double deltaRmin_Genjet_parton1 = 1000;
    int Id1_min = -1000;
    
    //       cout<< " " << endl;
    //cout<< "COMINCIA QUI"<<endl;
    
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



 /*    costruisco widejet prendendo i due jet che matcheano meglio i partoni (sui dati non si puo fare)

       double DeltaR1_min = 100;
       int index_jet1 =100;
       double DeltaR2_min = 100;
       int index_jet2 =100; 
       
       for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
      
      TLorentzVector currentJet;
      currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   
      
      //cout<<"Jet "<< ijet<<endl;
      //cout<<"Pt: "<<currentJet.Pt()<< " Eta: "<<currentJet.Eta()<<" Phi: "<<currentJet.Phi()<<endl;
      
      double DeltaR1 = currentJet.DeltaR(parton1);
      //      cout<<"DeltaR currentJet " << ijet<< " e Partone1 "<< DeltaR1<<endl;
      
      if(DeltaR1 < DeltaR1_min){
	DeltaR1_min = DeltaR1;
	jet1 = currentJet ; 
	index_jet1 = ijet;
      }

      double DeltaR2 = currentJet.DeltaR(parton2);
      //      cout<<"DeltaR currentJet " << ijet<< " e Partone2 "<< DeltaR2<<endl;

      if(DeltaR2 < DeltaR2_min){
	DeltaR2_min = DeltaR2;
	jet2 = currentJet ; 
	index_jet2 = ijet;
      }
    }
    //	cout<< "Jet numero "<<index_jet1<<" matchea il partone 1 con dR "<<DeltaR1_min<<endl;
    //	cout<< "Jet numero "<<index_jet2<<" matchea il partone 2 con dR "<<DeltaR2_min<<endl;
    
      */
   


  /*  //associo in modo diverso i widejets -> non c'e' bisogno
  //associo i Reco WideJets ai GenWideJets
    double DeltaR_Widejet1_tmp_GenWideJets1 = wj1_tmp.DeltaR(Genwj1);     
    double DeltaR_Widejet2_tmp_GenWideJets1 = wj2_tmp.DeltaR(Genwj1);       
    if(DeltaR_Widejet1_tmp_GenWideJets1 < DeltaR_Widejet2_tmp_GenWideJets1){
      cout<<"Widejet_tmp 1 associato al GenWideJet 1 Allora wj1 = wj1_tmp"<<endl;
      cout<<"Widejet_tmp 2 associato al GenWideJet 2 Allora wj2 = wj2_tmp"<<endl;
      wj1 = wj1_tmp;
      wj2 = wj2_tmp;
    }
    else{
      cout<<"Widejet_tmp 2 associato al GenWideJet 1 Allora wj1 = wj2_tmp"<<endl;
      cout<<"Widejet_tmp 1 associato al GenWideJet 2 Allora wj2 = wj1_tmp"<<endl;
      wj1 = wj2_tmp; 
      wj2 = wj1_tmp; 
    }
    double DeltaR_Gen_Reco[2];   
    DeltaR_Gen_Reco[0] = wj1.DeltaR(Genwj1);
    DeltaR_Gen_Reco[1] = wj2.DeltaR(Genwj2);
    cout<< "DeltaR tra widejet 1 e Genwidejet 1 = "<<DeltaR_Gen_Reco[0]<<endl;
    cout<< "DeltaR tra widejet 2 e Genwidejet 2 = "<<DeltaR_Gen_Reco[1]<<endl;
*/    


/*   //Funzioni di smearing calcolate con i partoni
   
     double parton_Pt[2];
    double parton_Eta[2];    
    parton_Pt[0]   = parton1.Pt();
    parton_Eta[0] = parton1.Eta();
    
    parton_Pt[1]   = parton2.Pt();
    parton_Eta[1] = parton2.Eta();    
    
    if(DeltaR[0] > 0.3 || DeltaR[1] > 0.3) cout<<"Uno dei due jet non matchea il partone. Butto evento"<<endl;
    if(DeltaR[0] > 0.3 || DeltaR[1] > 0.3) continue;
    
    CreateAndFillUserTH1D("H_DeltaR1_Match", 100, 0. , 3. , DeltaR[0]);
    CreateAndFillUserTH1D("H_DeltaR2_Match", 100, 0. , 3. , DeltaR[1]);
    
    for(int kk=0; kk<2; kk++){
      
      CreateAndFillUserTH1D("H_DeltaR_Totale", 100, 0. , 3. , DeltaR[kk]);
      CreateAndFillUserTH1D("H_R_PtWidejet_Ptparton", 100, 0. , 2. , R_PtWidejet_Ptparton[kk] );
      
      cout<<"Parton "<< kk+1<<endl;
      cout<<"Pt: "<<parton_Pt[kk]<<"         |Eta|: "<< fabs(parton_Eta[kk])<<endl;

      for(int ii=0; ii<5; ii++){	  
	double eta_bin_min = ii/2.;       
	double eta_bin_max = ii/2. +0.5;
	
	if( fabs(parton_Eta[kk]) >=eta_bin_min && fabs(parton_Eta[kk]) < eta_bin_max){	    

	  //mi fisso in un bin solo
	  //	if( fabs(parton_Eta[kk]) >=0 && fabs(parton_Eta[kk]) < 0.5){	    

	  cout<<"Passed bin in eta"<<endl;	
	  cout << "eta_bin: "<< eta_bin_min <<" - "<<eta_bin_max << endl;
	  
	  for(int jj=0; jj<n_bin ; jj++){	 	      
	    char HistoName[200];	 
	    int pt_bin_min = ( (jj-1)*step_pt)+step_pt;
	    int pt_bin_max = (jj*step_pt)+step_pt;
	    
	    if(parton_Pt[kk] >=pt_bin_min && parton_Pt[kk] < pt_bin_max){ 
	      cout<<"Passed bin in pt"<<endl;
	      cout << "pt_bin: "<< pt_bin_min <<" - "<<pt_bin_max<< endl;
	      
	      sprintf(HistoName,"Histo_RRecowide_Parton_%d_%d",ii,pt_bin_max);	 
	      
	      CreateAndFillUserTH1D(HistoName, histo_bin, xmin, xmax, R_PtWidejet_Ptparton[kk] );
	      cout<<"Riempio TH1D "<<HistoName<<" con R: "<<R_PtWidejet_Ptparton[kk]<<endl;
	      
	    }//if parton_Pt
	  }//ciclo sui bin in pT
	}//if parton_Eta
      }//ciclo sui bin di eta
    }//ciclo sui partoni
*/
