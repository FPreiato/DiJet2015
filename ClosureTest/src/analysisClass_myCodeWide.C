#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

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

   TH1F *h_nJetFinal = new TH1F ("h_nJetFinal","",10,0,10);
   h_nJetFinal->Sumw2();      
   TH1F *h_nVtx = new TH1F ("h_nVtx","",40,0,40);
   h_nVtx->Sumw2(); 
   TH1F *h_trueVtx = new TH1F ("h_trueVtx","",50,0,50);
   h_trueVtx->Sumw2();  
   TH1F *h_pT1stJet = new TH1F ("h_pT1stJet","",100,0,3000);
   h_pT1stJet->Sumw2();
   TH1F *h_pT2ndJet = new TH1F ("h_pT2ndJet","",100,0,3000);
   h_pT2ndJet->Sumw2();
   TH1F *h_eta1stJet = new TH1F ("h_eta1stJet","",5,-2.5,2.5);
   h_eta1stJet->Sumw2();
   TH1F *h_eta2ndJet = new TH1F ("h_eta2ndJet","",5,-2.5,2.5);
   h_eta2ndJet->Sumw2();
   TH1F *h_DijetMass = new TH1F ("h_DijetMass","",600,0,6000);
   h_DijetMass->Sumw2();
   TH1F *h_DeltaETAjj = new TH1F ("h_DeltaETAjj","",120,0,3.);
   h_DeltaETAjj->Sumw2();

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //for (Long64_t jentry=0; jentry<10;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done for every event

     Int_t no_jets_ak4=jetPtAK4->size();
     vector<TLorentzVector> widejets;
     TLorentzVector currentJet, wj1_tmp, wj1, wj2_tmp, wj2, v_jet_ak4,wdijet,wj1math,wj2math;
     double wideJetDeltaR_ = getPreCutValue1("DeltaR");

     resetCuts();

/*   
     if(no_jets_ak4>=2){
       if(!(fabs(jetEtaAK4->at(0)) < getPreCutValue1("jetFidRegion") && idTAK4->at(0) == getPreCutValue1("tightJetID"))){
         std::cout << " JET 0 FAIL " << jetEtaAK4->at(0) << " JET 0  ID " << idTAK4->at(0) << std::endl;
       }
       if(!(fabs(jetEtaAK4->at(1)) < getPreCutValue1("jetFidRegion") && idTAK4->at(1) == getPreCutValue1("tightJetID"))){
         std::cout << " JET 1 FAIL " << jetEtaAK4->at(1) << " JET 1  ID " << idTAK4->at(1) << std::endl;
       }  
     }
*/
     if(no_jets_ak4>=2){
       if(abs(jetEtaAK4->at(0)-jetEtaAK4->at(1)) < getPreCutValue1("deltaEtaCut")){
         if(fabs(jetEtaAK4->at(0)) < getPreCutValue1("jetFidRegion") && idTAK4->at(0) == getPreCutValue1("tightJetID") && jetPtAK4->at(0) > getPreCutValue1("pt0Cut")){
          if(fabs(jetEtaAK4->at(1)) < getPreCutValue1("jetFidRegion") && idTAK4->at(1) == getPreCutValue1("tightJetID") && jetPtAK4->at(1) > getPreCutValue1("pt1Cut")){
            TLorentzVector jet1, jet2;
            jet1.SetPtEtaPhiM(jetPtAK4->at(0),jetEtaAK4->at(0),jetPhiAK4->at(0),jetMassAK4->at(0));
            jet2.SetPtEtaPhiM(jetPtAK4->at(1),jetEtaAK4->at(1),jetPhiAK4->at(1),jetMassAK4->at(1));

            for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
              TLorentzVector currentJet;

              if(fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") && idTAK4->at(ijet) == getPreCutValue1("tightJetID") && jetPtAK4->at(ijet) > getPreCutValue1("ptCut")){
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

            //if(wj1_tmp.Pt()==0 && wj2_tmp.Pt() ==0) 
            //std::cout << " wj1_tmp.Pt() IN  " <<wj1_tmp.Pt()  << " wj2_tmp.Pt() " <<  wj2_tmp.Pt()  << std::endl;

         } //fid, jet id, pt cut
        } //fid, jet id, pt cut
      }//deltaEta cut
    } // end of two jets.

    //if(wj1_tmp.Pt()==0 && wj2_tmp.Pt() ==0)     std::cout << " wj1_tmp.Pt()  " <<wj1_tmp.Pt()  << " wj2_tmp.Pt() " << wj2_tmp.Pt()  << "  no_jets_ak4 " << no_jets_ak4 << std::endl;

    double MJJWide,DeltaEtaJJWide,DeltaPhiJJWide;

    if( wj1_tmp.Pt()>0 && wj2_tmp.Pt() >0 ){
       // Re-order the wide jets in pt
        if( wj1_tmp.Pt() > wj2_tmp.Pt()){
                wj1 = wj1_tmp;
                wj2 = wj2_tmp;
        }
        else{
                wj1 = wj2_tmp;
                wj2 = wj1_tmp;
        }

        // Create dijet system
        wdijet = wj1 + wj2;
        MJJWide = wdijet.M();
        DeltaEtaJJWide = fabs(wj1.Eta()-wj2.Eta());
        DeltaPhiJJWide = fabs(wj1.DeltaPhi(wj2));

        // Put widejets in the container
        wj1math.SetPtEtaPhiM(wj1.Pt(), wj1.Eta(), wj1.Phi(), wj1.M());
        wj2math.SetPtEtaPhiM(wj2.Pt(), wj2.Eta(), wj2.Phi(), wj2.M());  
        widejets.push_back( wj1math );
        widejets.push_back( wj2math );
   }

     fillVariableWithValue("nJetFinal",widejets.size());
     fillVariableWithValue("nVtx",nvtx);
     //cout<<PileupInteractions->at(12)<<endl;
     //fillVariableWithValue("trueVtx",PileupInteractions->at(12));

     if( widejets.size() >= 1 )
       {
         fillVariableWithValue( "pT1stJet", widejets[0].Pt() );
         fillVariableWithValue( "eta1stJet", widejets[0].Eta());

         fillVariableWithValue( "phi_j1", widejets[0].Phi());
	 fillVariableWithValue( "neutrHadEnFrac_j1", jetNhfAK4->at(0));
	 fillVariableWithValue( "chargedHadEnFrac_j1", jetChfAK4->at(0));
	 fillVariableWithValue( "photonEnFrac_j1", jetPhfAK4->at(0));
	 fillVariableWithValue( "eleEnFract_j1", jetElfAK4->at(0));
	 fillVariableWithValue( "muEnFract_j1", jetMufAK4->at(0));
       }

     if( widejets.size() >= 2 )
       {
         fillVariableWithValue( "pT2ndJet", widejets[1].Pt() );
         fillVariableWithValue( "eta2ndJet", widejets[1].Eta());

         // Calculate deltaETAjj
	 fillVariableWithValue( "deltaETAjj", DeltaEtaJJWide ) ;

         // Calculate Mjj
         fillVariableWithValue( "Dijet_MassW", wdijet.M() ) ;
         fillVariableWithValue( "Dijet_MassA", mjjAK8 ) ;  
         fillVariableWithValue( "Dijet_MassC", mjjCA8 ) ;
         if(wdijet.M()<1){
           std::cout << " INV MASS IS " << wdijet.M() << std::endl;
           std::cout << " Delta Eta IS " << DeltaEtaJJWide << " n is  " << widejets.size() << std::endl;
           std::cout << " INV MASS FROM NTUPLE AK8 " << mjjAK8 << std::endl;
           std::cout << " INV MASS FROM NTUPLE CA8 " << mjjCA8 << std::endl;
          }

        fillVariableWithValue( "phi_j2", widejets[1].Phi());
        fillVariableWithValue( "neutrHadEnFrac_j2", jetNhfAK4->at(1));
        fillVariableWithValue( "chargedHadEnFrac_j2", jetChfAK4->at(1));
        fillVariableWithValue( "photonEnFrac_j2", jetPhfAK4->at(1));
        fillVariableWithValue( "eleEnFract_j2", jetElfAK4->at(1));
        fillVariableWithValue( "muEnFract_j2", jetMufAK4->at(1));
        fillVariableWithValue( "deltaPHIjj", DeltaPhiJJWide ) ;

     }

     // Evaluate cuts (but do not apply them)
     evaluateCuts();

     // optional call to fill a skim with the full content of the input roottuple
     if( passedCut("nJetFinal") ) fillSkimTree();

     // optional call to fill a skim with a subset of the variables defined in the cutFile (use flag SAVE)
     if( passedCut("nJetFinal") ) fillReducedSkimTree();

     // reject events that did not pass level 0 cuts
     if( !passedCut("0") ) continue;
     // ......

     // reject events that did not pass level 1 cuts
     if( !passedCut("1") ) continue;
     // ......

     // reject events that did not pass the full cut list
     if( !passedCut("all") ) continue;
     // ......

     if( widejets.size() >= 2) {
      h_nJetFinal->Fill(widejets.size());
      h_DijetMass->Fill(wdijet.M());
      h_pT1stJet->Fill(widejets[0].Pt());
      h_pT2ndJet->Fill(widejets[1].Pt());
      h_eta1stJet->Fill(widejets[0].Eta());
      h_eta2ndJet->Fill(widejets[1].Eta());
     }
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 

   h_nVtx->Write();
   h_trueVtx->Write();
   h_nJetFinal->Write();
   h_pT1stJet->Write();
   h_pT2ndJet->Write();
   h_DijetMass->Write();
   h_eta1stJet->Write();
   h_eta2ndJet->Write();

   //pT of both jets, to be built using the histograms produced automatically by baseClass
   TH1F * h_pTJets = new TH1F ("h_pTJets","", getHistoNBins("pT1stJet"), getHistoMin("pT1stJet"), getHistoMax("pT1stJet"));
   h_pTJets->Add( & getHisto_noCuts_or_skim("pT1stJet") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   h_pTJets->Add( & getHisto_noCuts_or_skim("pT2ndJet") );
   //one could also do:  *h_pTJets = getHisto_noCuts_or_skim("pT1stJet") + getHisto_noCuts_or_skim("pT2ndJet");
   h_pTJets->Write();
   //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
