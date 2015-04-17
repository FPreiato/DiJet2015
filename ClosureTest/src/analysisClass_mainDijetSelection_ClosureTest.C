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


TFile file1("output/rootFile_QuarkQuark.root");
//Gev for pT binning
int step_pt = 500;

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
   
   char HistoName[200];	 
   TH1D *h[100][100];
   for(int ii=0; ii<5; ii++){
     //     step_pt = 500;
     int n_bin = 4500 / step_pt ;   
     for(int jj=0; jj<n_bin ; jj++){	 	 
       int pt_bin_max = (jj*step_pt)+step_pt;
       sprintf(HistoName,"Histo_RRecowide_Parton_%d_%d",ii,pt_bin_max);
       h[ii][jj] = new TH1D(HistoName, HistoName, 50,0.,2.);   
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

     if( (int)jentry%2 ==! 0 ) continue;

     //     cout << jentry << endl;

     ////////////////////////////////////////////////////////////////////////////
    // Partoni generati

     TLorentzVector p1, p2;
     TLorentzVector parton1_tmp, parton2_tmp;
     TLorentzVector parton1, parton2;
     TLorentzVector diparton;
     vector<TLorentzVector> partons;

     double diparton_MJJ;
     double diparton_DeltaEtaJJ;
     double diparton_cosThetaStar;
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

     parton1_tmp.SetPtEtaPhiM(parton1_pt, parton1_eta, parton1_phi, parton1_mass );
     parton2_tmp.SetPtEtaPhiM(parton2_pt, parton2_eta, parton2_phi, parton2_mass );

     // Re-order the partons in pt
     if( parton1_tmp.Pt() > parton2_tmp.Pt() ){
       parton1 = parton1_tmp;
       parton2 = parton2_tmp;
     }else{
       parton1 = parton2_tmp;
       parton2 = parton1_tmp;
     }
     
     if( parton1.Pt()<0 || parton2.Pt()<0 ) continue ;
     //     if( fabs(parton1.Eta() ) > 0.5 || fabs(parton2.Eta() ) > 0.5) continue; // see -> aggiunto questo
       
     // Create dijet(reco) system
     diparton = parton1 + parton2;
     diparton_MJJ = diparton.M();
     diparton_DeltaEtaJJ = fabs( parton1.Eta() - parton2.Eta() );
     diparton_cosThetaStar = tanh ( diparton_DeltaEtaJJ/2 ) ;
     diparton_DeltaPhiJJ = fabs(parton1.DeltaPhi(parton2));
      
     //   Put parton in the container
     partons.push_back( parton1 );
     partons.push_back( parton2 );
     
     //+++++++++++++++++++++++++++++++++++++++++++++++++++

    
    // Creo, riempio e salvo gli istogrammi
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
    CreateAndFillUserTH1D("H_diparton_M", 1000, 0. , 10000. , diparton_MJJ );
    CreateAndFillUserTH1D("H_diparton_DeltaEtaJJ", 100, 0. , 8. , diparton_DeltaEtaJJ );
    CreateAndFillUserTH1D("H_diparton_cosThetaStar", 50, 0. , 1. , diparton_cosThetaStar);
    //    CreateAndFillUserTH2D("H_diparton_Pt_Eta", 1000, 0., 10000., 120, -10. , 10. , diparton.Pt(), diparton.Eta() );
    
    //////////////////////////////////////////////////////////////////////////////////////////

    // mi calcolo il widejet partendo dai partoni con le funzioni di smearing
    
    double factor;
    //    TH1D *h;

    double parton_Pt[3];
    double parton_Eta[3];
    double recowidejet_Pt[3];

    parton_Pt[0] = parton1.Pt();
    parton_Pt[1] = parton2.Pt();
    parton_Eta[0] = parton1.Eta();
    parton_Eta[1] = parton2.Eta();
    
    for(int kk = 0 ; kk < 2 ; kk++){
      
      for(int ii=0; ii<5; ii++){
	
	double eta_bin_min = ii/2.;       
	double eta_bin_max = ii/2. +0.5;
	
	//	int step_pt = 500;

	int n_bin = 4500/ step_pt;
	
	if( fabs(parton_Eta[kk])>=eta_bin_min && fabs(parton_Eta[kk]) < eta_bin_max){
	  
	  for(int jj=0; jj<n_bin ; jj++){	 
	    
	    int pt_bin_min = ( (jj-1)*step_pt)+step_pt;
	    int pt_bin_max = (jj*step_pt)+step_pt;
	    
	    if(parton_Pt[kk] >=pt_bin_min && parton_Pt[kk] < pt_bin_max){ 
	      
	      sprintf(HistoName,"Histo_RRecowide_Parton_%d_%d",ii,pt_bin_max);	    
	      
	      //	      TH1D *h[100][100];
	      //h[ii][jj] = new TH1D(HistoName, HistoName, 50,0.,2.);   
	      h[ii][jj]= (TH1D*)file1.Get(HistoName);
	      
	      if( h[ii][jj]->Integral() == 0){
		
		factor = 1;
	      }else{
		
		factor = h[ii][jj]->GetRandom();
		
	      }
	      
	      //	      cout<< factor<<endl;
	      
	      recowidejet_Pt[kk] = parton_Pt[kk] * factor ;
	      
	    }
	  }
	}
      }
    }
    
    double RecoWideJet1_Pt = recowidejet_Pt[0];
    double RecoWideJet2_Pt = recowidejet_Pt[1];
    
    //    TLorentzVector RecoWideJet1 , RecoWideJet2 ;
    
    //    RecoWideJet1 = SetPtEtaPhiM()
    
    CreateAndFillUserTH1D("H_RecoWideJet1_Pt_Mio", 1000, 0. , 10000. ,RecoWideJet1_Pt); 
    CreateAndFillUserTH1D("H_RecoWideJet2_Pt_Mio", 1000, 0. , 10000. ,RecoWideJet2_Pt); 

    //////////////////////////////////////////////////////////////////////////

    //confronto quello uscito dalle funzioni di smearing con il MC
     //Livello reco dentro la ntupla -> MC
     //ricostruisce i widejet partendo di recojet_ak4
    
     size_t no_jets_ak4=jetPtAK4->size();
     
     TLorentzVector jet1, jet2, dijet;
     vector<TLorentzVector> jets;
     
     TLorentzVector wj1_tmp, wj2_tmp;
     TLorentzVector wj1, wj2, wdijet;
     vector<TLorentzVector> widejets;
     
     double wideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     double MJJ = -1000; 
     double DeltaEtaJJ = 1000;
     double Reco_cosThetaStar = 1000;
     double DeltaPhiJJ = 1000;
     
     ////////////////////////////////////////////////////////////

     if(no_jets_ak4<2) continue;

     //  TLorentzVector jet1, jet2;
     jet1.SetPtEtaPhiM(jetPtAK4->at(0),jetEtaAK4->at(0),jetPhiAK4->at(0),jetMassAK4->at(0));
     jet2.SetPtEtaPhiM(jetPtAK4->at(1),jetEtaAK4->at(1),jetPhiAK4->at(1),jetMassAK4->at(1));

     if( jet1.Pt()<0 || jet2.Pt()<0 ) continue;
     
     //costruisco reco widejet 
     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4

       //rimossi ragli anche qui       
       // if( fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
       // && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
       // && jetPtAK4->at(ijet) > getPreCutValue1("ptCut") ){
       
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
       //       } // if AK4 jet passes fid and jetid.
     } //end of ak4 jet loop		     
     
     // Re-order the wide jets in pt
     if( wj1_tmp.Pt() > wj2_tmp.Pt()){
       wj1 = wj1_tmp;
       wj2 = wj2_tmp;
     }
     else{
       wj1 = wj2_tmp;
       wj2 = wj1_tmp;
     }
     
     double MJJWide = -1000; 
     double DeltaEtaJJWide = 1000;
     double RecoWide_cosThetaStar = 1000;
     double DeltaPhiJJWide = 1000;
     
     if( wj1.Pt()<0 || wj2.Pt()<0 ) continue;
     
     // Create dijet system
     wdijet = wj1 + wj2;
     MJJWide = wdijet.M();
     DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
     RecoWide_cosThetaStar = tanh( DeltaEtaJJWide/2 );
     DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));
     
     // Put widejets in the container
     widejets.push_back( wj1 );
     widejets.push_back( wj2 );



     double difference_jet1 = (RecoWideJet1_Pt - wj1.Pt() ) / wj1.Pt() ;
     double difference_jet2 = (RecoWideJet2_Pt - wj2.Pt() ) / wj1.Pt() ;

    CreateAndFillUserTH1D("H_difference_jet1", 50, -1. , 5. , difference_jet1); 
    CreateAndFillUserTH1D("H_difference_jet2", 50, -1. , 5. , difference_jet2); 

    /////////////////////////////////////////////////////////////////////////////////////
    
    // Creo, riempio e salvo gli istogrammi     
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
     CreateAndFillUserTH1D("H_dijetWide_M", 1000, 0. , 10000. , MJJWide );
     CreateAndFillUserTH1D("H_DeltaEtaJJWide", 100, 0. , 8. , DeltaEtaJJWide );
     CreateAndFillUserTH1D("H_dijetWide_cosThetaStar", 50, 0. , 1. , RecoWide_cosThetaStar );
     
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
