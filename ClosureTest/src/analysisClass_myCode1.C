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

  //////////book histos here

   TH1F *h_nJetFinal = new TH1F ("h_nJetFinal","",10,0,10);
   h_nJetFinal->Sumw2();      
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

   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done every event
 
     Int_t no_jets_ak4=jetPtAK4->size();
     Int_t no_jets_ak8=jetPtAK8->size();
  
     vector<int> v_idx_jet_final;
     resetCuts();

     for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
      // h_ptak4->Fill(jetPtAK4->at(ijet));
     } //end of ak4 jet loop
   
     for(Long64_t ijet=0; ijet<no_jets_ak8; ijet++){ //jet loop for ak8
        if( fabs(jetEtaAK8->at(ijet)) < getPreCutValue1("jetFidRegion") && idTAK8->at(ijet) == getPreCutValue1("tightJetID")) v_idx_jet_final.push_back(ijet);
     } //end of ak8 jet loop

     fillVariableWithValue("nJetFinal", v_idx_jet_final.size()) ;

     TLorentzVector v_jet, jet1, jet2;

     if( v_idx_jet_final.size() >= 1 )

       {
         fillVariableWithValue( "pT1stJet", jetPtAK8->at(v_idx_jet_final[0]) );
         fillVariableWithValue( "eta1stJet", jetEtaAK8->at(v_idx_jet_final[0]) );
       }

     if( v_idx_jet_final.size() >= 2 )
       {
         fillVariableWithValue( "pT2ndJet", jetPtAK8->at(v_idx_jet_final[1]) );
         fillVariableWithValue( "eta2ndJet", jetEtaAK8->at(v_idx_jet_final[1]) );


         // Calculate Mjj
         jet1.SetPtEtaPhiM(jetPtAK8->at(v_idx_jet_final[0]),jetEtaAK8->at(v_idx_jet_final[0]),jetPhiAK8->at(v_idx_jet_final[0]),0);
         jet2.SetPtEtaPhiM(jetPtAK8->at(v_idx_jet_final[1]),jetEtaAK8->at(v_idx_jet_final[1]),jetPhiAK8->at(v_idx_jet_final[1]),0);
         v_jet = jet1 + jet2;
         fillVariableWithValue( "Dijet_Mass", v_jet.M() ) ;
         std::cout << " INV MASS IS " << v_jet.M() << std::endl;
         std::cout << " INV MASS FROM NTUPLE " << mjjAK8 << std::endl;

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

     if( v_idx_jet_final.size() >= 2) {
      h_nJetFinal->Fill(v_idx_jet_final.size());
      h_DijetMass->Fill(v_jet.M());
      h_pT1stJet->Fill(jetPtAK8->at(v_idx_jet_final[0]));
      h_pT2ndJet->Fill(jetPtAK8->at(v_idx_jet_final[1]));
      h_eta1stJet->Fill(jetEtaAK8->at(v_idx_jet_final[0]));
      h_eta2ndJet->Fill(jetEtaAK8->at(v_idx_jet_final[1]));


     }
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 


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
