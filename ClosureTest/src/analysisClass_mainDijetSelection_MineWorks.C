#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TMath.h>

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
   
   //////////book histos here

   // TH1F *h_nJetFinal = new TH1F ("h_nJetFinal","",10,0,10);
   // h_nJetFinal->Sumw2();      
   // TH1F *h_nVtx = new TH1F ("h_nVtx","",30,0,30);
   // h_nVtx->Sumw2(); 
   // TH1F *h_trueVtx = new TH1F ("h_trueVtx","",40,0,40);
   // h_trueVtx->Sumw2();  
   // TH1F *h_pT1stJet = new TH1F ("h_pT1stJet","",100,0,3000);
   // h_pT1stJet->Sumw2();
   // TH1F *h_pT2ndJet = new TH1F ("h_pT2ndJet","",100,0,3000);
   // h_pT2ndJet->Sumw2();
   // TH1F *h_eta1stJet = new TH1F ("h_eta1stJet","",5,-2.5,2.5);
   // h_eta1stJet->Sumw2();
   // TH1F *h_eta2ndJet = new TH1F ("h_eta2ndJet","",5,-2.5,2.5);
   // h_eta2ndJet->Sumw2();
   // TH1F *h_DijetMass = new TH1F ("h_DijetMass","",600,0,6000);
   // h_DijetMass->Sumw2();
   // TH1F *h_DeltaETAjj = new TH1F ("h_DeltaETAjj","",120,0,3.);
   // h_DeltaETAjj->Sumw2();

   TH1F *H_GenMJJ_hp = new TH1F ("H_GenMJJ_hp","",1000,0,10000);
   TH1F *H_GenMJJ = new TH1F ("H_GenMJJ","",1000,0,10000);
   TH1F *H_GenMJJWide = new TH1F ("H_GenMJJWide","",1000,0,10000);
   TH1F *H_MJJ = new TH1F ("H_MJJ","",1000,0,10000);
   TH1F *H_MJJWide = new TH1F ("H_MJJWide","",1000,0,10000);

   TH1F *H_GenDeltaEtaJJ_hp = new TH1F ("H_GenDeltaEtaJJ_hp","",100,0,3.);
   TH1F *H_GenDeltaEtaJJ = new TH1F ("H_GenDeltaEtaJJ","",100,0,3.);
   TH1F *H_GenDeltaEtaJJWide = new TH1F ("H_GenDeltaEtaJJWide","",100,0,3.);
   TH1F *H_DeltaEtaJJ = new TH1F ("H_DeltaEtaJJ","",100,0,3.);
   TH1F *H_DeltaEtaJJWide = new TH1F ("H_DeltaEtaJJWide","",100,0,3.);

   TH1F *H_Genjet1Eta_hp = new TH1F ("H_Genjet1Eta_hp","",100,-3.,3.);
   TH1F *H_Genjet2Eta_hp = new TH1F ("H_Genjet2Eta_hp","",100,-3.,3.);

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

     //     if(gen_pt->size() == 3){
     //cout<< "gen_pt->size():" << gen_pt->size() << endl;
     // cout<< "gen_pt->at(0):" << gen_pt->at(0) << endl;
     //cout<< "gen_pt->at(1):" << gen_pt->at(1) << endl;
     //cout<< "gen_pt->at(2):" << gen_pt->at(2) << endl;
     //     cout<< "gen_pt->at(3):" << gen_pt->at(3) << endl;
     //     cout<< "gen_pt->at(4):" << gen_pt->at(4) << endl;
     //}

     
     // Generated jet from hard process

     TLorentzVector Genjet1_hp, Genjet2_hp;
     TLorentzVector Gendijet_hp;
     double GenMJJ_hp=0;
     double GenDeltaEtaJJ_hp =0;
     double GenDeltaPhiJJ_hp=0;
     float gen_mass[5];

     if(gen_pt->size() == 3){

       gen_mass[1] =sqrt( gen_energy->at(1)*gen_energy->at(1) - gen_p->at(1)*gen_p->at(1) );
       gen_mass[2] =sqrt( gen_energy->at(2)*gen_energy->at(2) - gen_p->at(2)*gen_p->at(2) );

       Genjet1_hp.SetPtEtaPhiM(gen_pt->at(1),gen_eta->at(1),gen_phi->at(1),gen_mass[1] );
       Genjet2_hp.SetPtEtaPhiM(gen_pt->at(2),gen_eta->at(2),gen_phi->at(2),gen_mass[2] );
     }

     if(gen_pt->size() == 5){

       gen_mass[3] = gen_energy->at(3)*gen_energy->at(3) - gen_p->at(3)*gen_p->at(3);
       gen_mass[4] = gen_energy->at(4)*gen_energy->at(4) - gen_p->at(4)*gen_p->at(4);

       Genjet1_hp.SetPtEtaPhiM(gen_pt->at(3),gen_eta->at(3),gen_phi->at(3),gen_mass[3] );
       Genjet2_hp.SetPtEtaPhiM(gen_pt->at(4),gen_eta->at(4),gen_phi->at(4),gen_mass[4] );
     }

     if( Genjet1_hp.Pt()>0 && Genjet2_hp.Pt()>0 )
       { 
		     // Create dijet(reco) system
		     Gendijet_hp = Genjet1_hp + Genjet2_hp;
		     GenMJJ_hp = Gendijet_hp.M();
		     GenDeltaEtaJJ_hp = fabs( Genjet1_hp.Eta() - Genjet2_hp.Eta() );
		     GenDeltaPhiJJ_hp = fabs(Genjet1_hp.DeltaPhi(Genjet2_hp));

		     //		     if(Genjet1_hp.Eta()<0.5 && Genjet1_hp.Eta()>-0.5 || Genjet2_hp.Eta()<0.5 && Genjet2_hp.Eta()>-0.5){
		     //		     if(GenDeltaEtaJJ_hp < 0.5){
		     //cout<< "Genjet1_hp.Eta(): " << Genjet1_hp.Eta() << endl;
		     //cout<< "Genjet2_hp.Eta(): " << Genjet2_hp.Eta() << endl;
		     //cout<< "DeltaEtaJJ: " << GenDeltaEtaJJ_hp << endl;
		     //}

		     H_Genjet1Eta_hp-> Fill(Genjet1_hp.Eta());
		     H_Genjet2Eta_hp-> Fill(Genjet2_hp.Eta());
		    
		     H_GenMJJ_hp -> Fill(GenMJJ_hp);
		     H_GenDeltaEtaJJ_hp -> Fill(GenDeltaEtaJJ_hp);
  
		     // Put widejets in the container -> see
		     //		     Genjets.push_back( Genjet1 );
		     //Genjets.push_back( Genjet2 );
		   }


     // Generated jet from parton shower

     size_t no_Genjets_ak4=jetPtGenAK4->size();
     vector<TLorentzVector> Genwidejets;
     TLorentzVector Genwj1, Genwj2, Genwdijet;
     vector<TLorentzVector> Genjets;
     TLorentzVector Genjet1, Genjet2, Gendijet;

     TLorentzVector Genwj1_tmp, Genwj2_tmp;
     double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");

     double GenMJJ = 0; 
     double GenDeltaEtaJJ = 0;
     double GenDeltaPhiJJ = 0;

     if(no_Genjets_ak4>=2)
       {
	 if(fabs(jetEtaGenAK4->at(0)) < getPreCutValue1("jetFidRegion") 
	    && jetPtGenAK4->at(0) > getPreCutValue1("pt0Cut"))
	   {
	     if(fabs(jetEtaGenAK4->at(1)) < getPreCutValue1("jetFidRegion") 
		&& jetPtGenAK4->at(1) > getPreCutValue1("pt1Cut"))
	       {
		 //  TLorentzVector Genjet1, Genjet2;
		 Genjet1.SetPtEtaPhiM(jetPtGenAK4->at(0),jetEtaGenAK4->at(0),jetPhiGenAK4->at(0),jetMassGenAK4->at(0));
		 Genjet2.SetPtEtaPhiM(jetPtGenAK4->at(1),jetEtaGenAK4->at(1),jetPhiGenAK4->at(1),jetMassGenAK4->at(1));

		 //costruisco la massa invariante di questi due jet (prima di fare le somme per il widejet)
		 
		 if( Genjet1.Pt()>0 && Genjet2.Pt()>0 )
		   {
		     // Create dijet(reco) system
		     Gendijet = Genjet1 + Genjet2;
		     GenMJJ = Gendijet.M();
		     GenDeltaEtaJJ = fabs(Genjet1.Eta()-Genjet2.Eta());
		     GenDeltaPhiJJ = fabs(Genjet1.DeltaPhi(Genjet2));
		 
		     H_GenMJJ -> Fill(GenMJJ);
		     H_GenDeltaEtaJJ -> Fill(GenDeltaEtaJJ);
    
		     // Put widejets in the container
		     Genjets.push_back( Genjet1 );
		     Genjets.push_back( Genjet2 );
		   }
		 
		 // Generated jet clustering -> wide genjet
		 for(Long64_t ijet=0; ijet<no_Genjets_ak4; ijet++)
		   { //Genjet loop for ak4
		     TLorentzVector currentGenJet;
		     
		     if(fabs(jetEtaGenAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
			//&& idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
		       	&& jetPtGenAK4->at(ijet) > getPreCutValue1("ptCut"))
		       {
			 TLorentzVector currentGenJet;
			 currentGenJet.SetPtEtaPhiM(jetPtGenAK4->at(ijet),jetEtaGenAK4->at(ijet),jetPhiGenAK4->at(ijet),jetMassGenAK4->at(ijet));   
			 
			 double DeltaR1 = currentGenJet.DeltaR(Genjet1);
			 double DeltaR2 = currentGenJet.DeltaR(Genjet2);

			 if(DeltaR1 < DeltaR2 && DeltaR1 < GenwideJetDeltaR_)
			   {
			     Genwj1_tmp += currentGenJet;
			   }
			 else if(DeltaR2 < GenwideJetDeltaR_)
			   {
			     Genwj2_tmp += currentGenJet;
			   }			 
		       } // if AK4 jet passes fid and jetid.
		   } //end of ak4 jet loop		     

	       } //fid, jet id, pt cut
	   } //fid, jet id, pt cut
       } // end of two jets.

     // Re-order the wide genjets in pt
     if( Genwj1_tmp.Pt() > Genwj2_tmp.Pt())
       {
	 Genwj1 = Genwj1_tmp;
	 Genwj2 = Genwj2_tmp;
       }
     else
       {
	 Genwj1 = Genwj2_tmp;
	 Genwj2 = Genwj1_tmp;
       }

     double GenMJJWide = 0; 
     double GenDeltaEtaJJWide = 0;
     double GenDeltaPhiJJWide = 0;
     if( Genwj1.Pt()>0 && Genwj2.Pt()>0 )
     {
       // Create dijet system
       Genwdijet = Genwj1 + Genwj2;
       GenMJJWide = Genwdijet.M();
       GenDeltaEtaJJWide = fabs(Genwj1.Eta() - Genwj2.Eta());
       GenDeltaPhiJJWide = fabs(Genwj1.DeltaPhi(Genwj2));

       H_GenMJJWide -> Fill(GenMJJWide);
       H_GenDeltaEtaJJWide -> Fill(GenDeltaEtaJJWide);

       // Put widejets in the container
       Genwidejets.push_back( Genwj1 );
       Genwidejets.push_back( Genwj2 );
     }

     //Livello reco
     //ricostruisce i widejet partendo di recojet_ak4
 
     size_t no_jets_ak4=jetPtAK4->size();
     vector<TLorentzVector> widejets;
     TLorentzVector wj1, wj2, wdijet;
     vector<TLorentzVector> jets;
     TLorentzVector jet1, jet2, dijet;
   
     TLorentzVector wj1_tmp, wj2_tmp;
     double wideJetDeltaR_ = getPreCutValue1("DeltaR");

     double MJJ = 0; 
     double DeltaEtaJJ = 0;
     double DeltaPhiJJ = 0;
     
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

		     H_MJJ -> Fill(MJJ);		     
		     H_DeltaEtaJJ -> Fill(DeltaEtaJJ);		     

		     // Put widejets in the container
		     jets.push_back( jet1 );
		     jets.push_back( jet2 );
		   }
		 
		 for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++)
		   { //jet loop for ak4
		     TLorentzVector currentJet;
		     
		     if(fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
			&& idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
			&& jetPtAK4->at(ijet) > getPreCutValue1("ptCut"))
		       {
			 TLorentzVector currentJet;
			 currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet),jetEtaAK4->at(ijet),jetPhiAK4->at(ijet),jetMassAK4->at(ijet));   
			 
			 double DeltaR1 = currentJet.DeltaR(jet1);
			 double DeltaR2 = currentJet.DeltaR(jet2);
			 
			 if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_)
			   {
			     wj1_tmp += currentJet;
			   }
			 else if(DeltaR2 < wideJetDeltaR_)
			   {
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

       H_MJJWide -> Fill(MJJWide);
       H_DeltaEtaJJWide -> Fill(DeltaEtaJJWide);

       // Put widejets in the container
       widejets.push_back( wj1 );
       widejets.push_back( wj2 );
     }

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
	 fillVariableWithValue( "GendeltaETAjj_hp", GenDeltaEtaJJ_hp ) ;
         fillVariableWithValue( "Genmjj_hp", GenMJJ_hp ) ;
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
   

 
  TCanvas *c1 = new TCanvas("c1","",600,400);
   H_GenMJJ_hp->Draw();
   c1 -> SaveAs("output/Plot/GenMJJ_hp.png");

  TCanvas *c1b = new TCanvas("c1b","",600,400);
   H_GenDeltaEtaJJ_hp->Draw();
   c1b -> SaveAs("output/Plot/GenDeltaEtaJJ_hp.png");

  TCanvas *c1c = new TCanvas("c1c","",600,400);
   H_Genjet1Eta_hp->Draw();
   c1c -> SaveAs("output/Plot/Genjet1Eta_hp.png");

  TCanvas *c1d = new TCanvas("c1d","",600,400);
   H_Genjet2Eta_hp->Draw();
   c1d -> SaveAs("output/Plot/Genjet2Eta_hp.png");

  TCanvas *c2 = new TCanvas("c2","",600,400);
   H_GenMJJ->Draw();
   c2 -> SaveAs("output/Plot/GenMJJ.png");

   TCanvas *c2b = new TCanvas("c2b","",600,400);
   H_GenDeltaEtaJJ->Draw();
   c2b -> SaveAs("output/Plot/GenDeltaEtaJJ.png");

  TCanvas *c3 = new TCanvas("c3","",600,400);
   H_GenMJJWide->Draw();
   c3 -> SaveAs("output/Plot/GenMJJWide.png");

   TCanvas *c3b = new TCanvas("c3b","",600,400);
   H_DeltaEtaJJWide->Draw();
   c3b -> SaveAs("output/Plot/DeltaEtaJJWide.png");

  TCanvas *c4 = new TCanvas("c4","",600,400);
   H_MJJ->Draw();
   c4 -> SaveAs("output/Plot/MJJ.png");

   TCanvas *c4b = new TCanvas("c4b","",600,400);
   H_DeltaEtaJJ->Draw();
   c4b -> SaveAs("output/Plot/DeltaEtaJJ.png");

  TCanvas *c5 = new TCanvas("c5","",600,400);
   H_MJJWide -> Draw();
   c5 -> SaveAs("output/Plot/MJJWide.png");

   TCanvas *c5b = new TCanvas("c5b","",600,400);
   H_DeltaEtaJJWide->Draw();
   c5b -> SaveAs("output/Plot/DeltaEtaJJWide.png");


   TLegend *leg1 = new TLegend( 0.15, 0.9, 0.5, 0.6  );
   leg1 -> AddEntry(H_GenMJJ_hp, " Genjets Hard Process","L");
   leg1 -> AddEntry(H_GenMJJ, " Genjets ak4","L");
   leg1 -> AddEntry(H_GenMJJWide, " Genjets Wide","L");
   leg1 -> AddEntry(H_MJJ, "Recojets ak4","L");
   leg1 -> AddEntry(H_MJJWide, "Recojets Wide","L");

  TCanvas *c6 = new TCanvas("c6","",600,400);
  H_GenMJJ_hp -> SetStats(0);
  H_GenMJJ_hp -> SetYTitle("Normalized to unit");
  H_GenMJJ_hp -> SetXTitle("M(jj) [GeV]");
  H_GenMJJ_hp -> SetTitle(" ");

  H_GenMJJ_hp     -> SetLineColor(kBlack);
  H_GenMJJ           -> SetLineColor(kGreen);
  H_GenMJJWide   -> SetLineColor(kBlue);
  H_MJJ                 -> SetLineColor(kRed);
  H_MJJWide         -> SetLineColor(kYellow+8);

  H_GenMJJ_hp     -> DrawNormalized();
  H_GenMJJ           -> DrawNormalized("same");
  H_GenMJJWide   -> DrawNormalized("same");
  H_MJJ                 -> DrawNormalized("same");
  H_MJJWide         -> DrawNormalized("same");

  leg1 -> Draw();

  c6 -> SaveAs("output/Plot/All_MJJ.png");

   TLegend *leg2 = new TLegend( 0.9, 0.9, 0.7, 0.6  );
   leg2 -> AddEntry(H_GenDeltaEtaJJ_hp, " Genjets Hard Process","L");
   leg2 -> AddEntry(H_GenDeltaEtaJJ, " Genjets ak4","L");
   leg2 -> AddEntry(H_GenDeltaEtaJJWide, " Genjets Wide","L");
   leg2 -> AddEntry(H_DeltaEtaJJ, "Recojets ak4","L");
   leg2 -> AddEntry(H_DeltaEtaJJWide, "Recojets Wide","L");


  TCanvas *c6b = new TCanvas("c6b","",600,400);
  H_GenDeltaEtaJJ_hp -> SetStats(0);
  H_GenDeltaEtaJJ_hp -> SetYTitle("Normalized to unit");
  H_GenDeltaEtaJJ_hp -> SetXTitle("Delta eta");
  H_GenDeltaEtaJJ_hp -> SetTitle(" ");

  H_GenDeltaEtaJJ_hp     -> SetLineColor(kBlack);
  H_GenDeltaEtaJJ           -> SetLineColor(kGreen);
  H_GenDeltaEtaJJWide   -> SetLineColor(kBlue);
  H_DeltaEtaJJ                 -> SetLineColor(kRed);
  H_DeltaEtaJJWide         -> SetLineColor(kYellow+8);

  H_GenDeltaEtaJJ_hp     -> DrawNormalized();
  H_GenDeltaEtaJJ           -> DrawNormalized("same");
  H_GenDeltaEtaJJWide   -> DrawNormalized("same");
  H_DeltaEtaJJ                 -> DrawNormalized("same");
  H_DeltaEtaJJWide         -> DrawNormalized("same");

  leg2 -> Draw();

  c6b -> SaveAs("output/Plot/All_DeltaEta.png");

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
