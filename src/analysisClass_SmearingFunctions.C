#define analysisClass_cxx
#include "analysisClass.h"
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TParameter.h>
#include <TH2D.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TMath.h>
#include "ptBinning.h"
#include "etaBinning.h"

#include "TFileService.h"

bool verbose = false;
int n_categories = 2;
int Ev_Initial;
int Ev_NPartons;
int Ev_PartonsOK;
int Ev_NGenJet;
int Ev_GenJetOK;
int Ev_NRecoJet;
int Ev_RecoJetOK;
int Ev_Matching;

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

  outputFileNameSmearing_ = outputFileName;

 // For JECs
  if( int(getPreCutValue1("useJECs"))==1 )
  {
    std::cout << "Reapplying JECs on the fly" << std::endl;
    //std::string L1Path = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_6_patch6/src/CMSDIJET/DijetRootTreeMaker/data/Summer15_V5/Summer15_V5_MC_L1FastJet_AK4PFchs.txt";
    //std::string L2Path = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_6_patch6/src/CMSDIJET/DijetRootTreeMaker/data/Summer15_V5/Summer15_V5_MC_L2Relative_AK4PFchs.txt"; 
    //std::string L3Path = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_6_patch6/src/CMSDIJET/DijetRootTreeMaker/data/Summer15_V5/Summer15_V5_MC_L3Absolute_AK4PFchs.txt";
    //
    std::string L1Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt";
    std::string L2Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt"; 
    std::string L3Path = "data/Summer15_25nsV6_MC/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt";
    std::string L1DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L1FastJet_AK4PFchs.txt";
    std::string L2DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2Relative_AK4PFchs.txt"; 
    std::string L3DATAPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L3Absolute_AK4PFchs.txt";
    std::string L2L3ResidualPath = "data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_L2L3Residual_AK4PFchs.txt" ;
    
    L1Par = new JetCorrectorParameters(L1Path);
    L2Par = new JetCorrectorParameters(L2Path);
    L3Par = new JetCorrectorParameters(L3Path);
    L1DATAPar = new JetCorrectorParameters(L1DATAPath);
    L2DATAPar = new JetCorrectorParameters(L2DATAPath);
    L3DATAPar = new JetCorrectorParameters(L3DATAPath);
    L2L3Residual = new JetCorrectorParameters(L2L3ResidualPath);

    std::vector<JetCorrectorParameters> vPar;
    std::vector<JetCorrectorParameters> vPar_data;
    vPar.push_back(*L1Par);
    vPar.push_back(*L2Par);
    vPar.push_back(*L3Par);
   
    //residuals are applied only to data
    vPar_data.push_back(*L1DATAPar);
    vPar_data.push_back(*L2DATAPar);
    vPar_data.push_back(*L3DATAPar);
    vPar_data.push_back(*L2L3Residual);

    JetCorrector = new FactorizedJetCorrector(vPar);
    JetCorrector_data = new FactorizedJetCorrector(vPar_data);

    //uncertainty
    //    unc = new JetCorrectionUncertainty("data/Summer15_25nsV6_DATA/Summer15_25nsV6_DATA_Uncertainty_AK4PFchs.txt");
  }

  std::cout << "analysisClass::analysisClass(): ends " << std::endl;

}

analysisClass::~analysisClass()
{
  TH1D *EventCounter = new TH1D("EventCounter","",10, 0 , 10);

  EventCounter -> GetXaxis()->SetBinLabel(1,"Initial ev.");
  EventCounter -> GetXaxis()->SetBinLabel(2,"N(partons)>2");
  EventCounter -> GetXaxis()->SetBinLabel(3,"Partons ok");
  EventCounter -> GetXaxis()->SetBinLabel(4,"N(GenJet)>2");
  EventCounter -> GetXaxis()->SetBinLabel(5,"GenJet OK");
  EventCounter -> GetXaxis()->SetBinLabel(6,"N(RecoJet)>2");
  EventCounter -> GetXaxis()->SetBinLabel(7,"RecoJet OK");
  EventCounter -> GetXaxis()->SetBinLabel(8,"Matching");
  
  EventCounter -> SetBinContent(1,Ev_Initial);
  EventCounter -> SetBinContent(2,Ev_NPartons);
  EventCounter -> SetBinContent(3,Ev_PartonsOK);
  EventCounter -> SetBinContent(4,Ev_NGenJet);
  EventCounter -> SetBinContent(5,Ev_GenJetOK);
  EventCounter -> SetBinContent(6,Ev_NRecoJet);
  EventCounter -> SetBinContent(7,Ev_RecoJetOK);
  EventCounter -> SetBinContent(8,Ev_Matching);

  EventCounter -> Write();

  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;
  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
   std::cout << "Running smearing function analysis" <<std::endl;   
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
    
   if (fChain == 0) return;


   if( n_categories != 2 && n_categories !=3 && n_categories != 4 && n_categories != 6){
     cout<<" "<<endl;
     cout<<"Wrong number of categories"<<endl;
     cout<<"Choose between 2, 3, 4 or 6"<<endl;
     exit(1);
   } else {
     cout<<"Number of categories choosen: "<<n_categories << endl;
   }  

   TFile *outputSmearing_root = new TFile((*outputFileNameSmearing_ + "_smearing.root").c_str(),"RECREATE");
   fwlite::TFileService fs(outputSmearing_root);
   Step_pt.SetVal(mPtBinning.get_PtStep()) ;
   cout<<Step_pt.GetVal() <<endl;
   NCategory.SetVal(n_categories) ;
   cout<<NCategory.GetVal() <<endl;
   NFile.SetVal(1) ;
   cout<<NFile.GetVal()<<endl;

   Step_pt.Write("PtStep");
   NCategory.Write("n_categories");
   NFile.Write("n_files");
   TFileDirectory smearingDir = fs.mkdir("smearingFunction");
   std::vector<std::vector<TH1F*> > smearingFunction_Q ;
   std::vector<std::vector<TH1F*> > smearingFunction_QUp ;
   std::vector<std::vector<TH1F*> > smearingFunction_QDown ;
   std::vector<std::vector<TH1F*> > smearingFunction_QStrange ;
   std::vector<std::vector<TH1F*> > smearingFunction_QBottom ;
   std::vector<std::vector<TH1F*> > smearingFunction_QCharm ;
   std::vector<std::vector<TH1F*> > smearingFunction_QLight ;
   std::vector<std::vector<TH1F*> > smearingFunction_QHeavy ;

   if(n_categories ==2){
     smearingFunction_Q= buildEtaPtVector<TH1F>(smearingDir, "smearing_Quark", 150, 0., 2.);   
   } else if(n_categories ==3){
     smearingFunction_QLight   = buildEtaPtVector<TH1F>(smearingDir, "smearing_LightQuark", 150, 0., 2.);   
     smearingFunction_QHeavy = buildEtaPtVector<TH1F>(smearingDir, "smearing_HeavyQuark", 150, 0., 2.);   
   } else if(n_categories ==4){
     smearingFunction_QLight    = buildEtaPtVector<TH1F>(smearingDir, "smearing_LightQuark", 150, 0., 2.);   
     smearingFunction_QCharm  = buildEtaPtVector<TH1F>(smearingDir, "smearing_CharmQuark", 150, 0., 2.);   
     smearingFunction_QBottom = buildEtaPtVector<TH1F>(smearingDir, "smearing_BottomQuark", 150, 0., 2.);   
   } else if(n_categories ==6){
     smearingFunction_QUp       = buildEtaPtVector<TH1F>(smearingDir, "smearing_UpQuark", 150, 0., 2.);   
     smearingFunction_QDown   = buildEtaPtVector<TH1F>(smearingDir, "smearing_DownQuark", 150, 0., 2.);   
     smearingFunction_QStrange = buildEtaPtVector<TH1F>(smearingDir, "smearing_StrangeQuark", 150, 0., 2.);   
     smearingFunction_QCharm  = buildEtaPtVector<TH1F>(smearingDir, "smearing_CharmQuark", 150, 0., 2.);   
     smearingFunction_QBottom = buildEtaPtVector<TH1F>(smearingDir, "smearing_BottomQuark", 150, 0., 2.);   
   }

   std::vector<std::vector<TH1F*> > smearingFunction_G = buildEtaPtVector<TH1F>(smearingDir, "smearing_Gluon", 150, 0., 2.);   
   output_root_ -> cd();   
   
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

     //find intime BX
     int idx_InTimeBX=-1;
     //     cout<<"PileupOriginBX.size() = "<< PileupOriginBX->size() << endl;
     for(size_t j=0; j<PileupOriginBX->size(); ++j)
       {
	 //	 cout <<"PileUpOriginBX->at(j)" << PileupOriginBX->at(j) << endl;	 
	 if(PileupOriginBX->at(j)==0)
	   {
	     idx_InTimeBX = j;
	     //cout << "idx_InTimeBX: " << idx_InTimeBX << endl; 
	   }
       }
     
     std::vector<double> jecFactors;
     std::vector<double> jecUncertainty;
     // new JECs could change the jet pT ordering. the vector below
     // holds sorted jet indices after the new JECs had been applied
     std::vector<unsigned> sortedJetIdx;
     bool isData = 0;
     //  cout<<"idx_InTimeBX = "<< idx_InTimeBX << endl;
     if(idx_InTimeBX > -1 ) isData = 0;
     else isData = 1;
     
     //     cout<<"isData ?"<< isData << endl;

     Ev_Initial++;

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

     Ev_NPartons++;

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
     
     Ev_PartonsOK++;

     if(verbose){
       cout<<endl;
       cout<<"PARTON 1"<<endl;
       cout<<"Pt: "<<parton1.Pt()<< " Eta: "<<parton1.Eta()<<" Phi: "<<parton1.Phi()<<" M: "<<parton1.M()<<" pdgId: "<<parton1_pdgId<<endl;
       cout<<"PARTON 2"<<endl;
       cout<<"Pt: "<<parton2.Pt()<< " Eta: "<<parton2.Eta()<<" Phi: "<<parton2.Phi()<<" M: "<<parton2.M()<<" pdgId: "<<parton2_pdgId<<endl;
     }
     // Create dijet(reco) system
     diparton = parton1 + parton2;
     
     //++++++++++++++++++++++++++++
     // Gen-level : Construct the GenWidejet from genjet_ak4
     size_t no_Genjets_ak4 = jetPtGenAK4->size();
     TLorentzVector Genjet1, Genjet2;
     TLorentzVector Genwj1_tmp, Genwj2_tmp;
     TLorentzVector Genwj1, Genwj2;
     TLorentzVector Genwdijet;
     double GenwideJetDeltaR_ = getPreCutValue1("DeltaR");
     
     if(no_Genjets_ak4 < 2) continue;

     Ev_NGenJet++;

     if( fabs(jetEtaGenAK4->at(0)) > getPreCutValue1("jetFidRegion") || jetPtGenAK4->at(0) < getPreCutValue1("pt0Cut") ) continue;
     if( fabs(jetEtaGenAK4->at(1)) > getPreCutValue1("jetFidRegion") || jetPtGenAK4->at(1) < getPreCutValue1("pt1Cut") ) continue;

     Ev_GenJetOK++;

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
       
       //+++++++++++++++++++++++++
       // Reco level : Construct  widejet from recojet_ak4
       size_t no_jets_ak4=jetPtAK4->size();    
       TLorentzVector jet1, jet2;
       TLorentzVector wj1_tmp, wj2_tmp;
       TLorentzVector wj1, wj2; 
       TLorentzVector wdijet;
       double wideJetDeltaR_ = getPreCutValue1("DeltaR");
      
       
       //////////////////// JEC
       if( int(getPreCutValue1("useJECs"))==1 )
	 {
	   // sort jets by increasing pT
	   std::multimap<double, unsigned> sortedJets;
	   for(size_t j=0; j<no_jets_ak4; ++j)
	     {
	       JetCorrector->setJetEta(jetEtaAK4->at(j));
	       JetCorrector->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)); //pTraw
	       JetCorrector->setJetA(jetAreaAK4->at(j));
	       JetCorrector->setRho(rho);
	       
	       JetCorrector_data->setJetEta(jetEtaAK4->at(j));
	       JetCorrector_data->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)); //pTraw
	       JetCorrector_data->setJetA(jetAreaAK4->at(j));
	       JetCorrector_data->setRho(rho);
	       
	       //nominal value of JECs
	       double correction;//, old_correction, nominal_correction;
	       //if( int(getPreCutValue1("shiftJECs"))==0 ){
	       if (isData == 1) correction = JetCorrector_data->getCorrection();
	       else correction = JetCorrector->getCorrection();
	       //nominal_correction=correction;
	       //old_correction = jetJecAK4->at(j);
	       //}
	       //JEC uncertainties
	       //	       unc->setJetEta(jetEtaAK4->at(j));
	       //  unc->setJetPt(jetPtAK4->at(j)/jetJecAK4->at(j)*correction);
	       // double uncertainty = unc->getUncertainty(true);
	       // jecUncertainty.push_back(uncertainty); 
	       
	       // std::cout << "run:" << runNo << "    lumi:" << lumi << "   event:" << evtNo << "   jet pt:" << jetPtAK4->at(j)/jetJecAK4->at(j)*correction << "   correction:" << correction <<   "   uncertainty:" <<  uncertainty  << "  nominal correction:" << nominal_correction  << " old correction: " << old_correction << std::endl;
	       //use "shifted" JECs for study of systematic uncertainties 
	       if( int(getPreCutValue1("shiftJECs"))==1 ){
	       //flat shift
	       //if (isData == 1) correction = JetCorrector_data->getCorrection() * getPreCutValue2("shiftJECs");
	       //else correction = JetCorrector->getCorrection() * getPreCutValue2("shiftJECs");
	       //shift of the corresponding unc
	       // ///////////correction = correction + getPreCutValue2("shiftJECs")*uncertainty*correction;
	       //  std::cout << "run:" << runNo << "    lumi:" << lumi << "   event:" << evtNo << "   jet pt:" << jetPtAK3->at(j)/jetJecAK4->at(j)*correction << "   correction:" << correction << "   uncertainty:" <<  uncertainty  << std::endl << std::endl;
	       
	       }
	       
	     jecFactors.push_back(correction);
	     sortedJets.insert(std::make_pair((jetPtAK4->at(j)/jetJecAK4->at(j))*correction, j));
	     
	     }
	 // get jet indices in decreasing pT order
	   for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	     sortedJetIdx.push_back(it->second);
	   
	 }
       else if( int(getPreCutValue1("noJECs"))==1  )
	 {
	   // sort jets by increasing pT
	   std::multimap<double, unsigned> sortedJets;
	   for(size_t j=0; j<no_jets_ak4; ++j) //same ordering of original root trees
	     {
	       //	       jecUncertainty.push_back(0.); 
	       jecFactors.push_back(1.);
	       sortedJets.insert(std::make_pair((jetPtAK4->at(j)/jetJecAK4->at(j)), j)); //raw
	     }       
	   // get jet indices in decreasing pT order
	   for(std::multimap<double, unsigned>::const_reverse_iterator it = sortedJets.rbegin(); it != sortedJets.rend(); ++it)
	     sortedJetIdx.push_back(it->second);
	 }
       else
	 {
	   for(size_t j=0; j<no_jets_ak4; ++j) //same ordering of original root trees
	     {
	       jecFactors.push_back(jetJecAK4->at(j));
	       //	       jecUncertainty.push_back(0.); 
	       sortedJetIdx.push_back(j);
	     }
	 }
       ///////////////////////////////////////////////////

       // JEC applied
       
       if(no_jets_ak4<2) continue;

       Ev_NRecoJet++;
  
       if(fabs(jetEtaAK4->at(sortedJetIdx[0])) > getPreCutValue1("jetFidRegion") || (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0]))*jetPtAK4->at(sortedJetIdx[sortedJetIdx[0]]) < getPreCutValue1("pt0Cut") ) continue;
       if(fabs(jetEtaAK4->at(sortedJetIdx[1])) > getPreCutValue1("jetFidRegion")  || (jecFactors[sortedJetIdx[1]]/jetJecAK4->at(sortedJetIdx[1]))*jetPtAK4->at(sortedJetIdx[1]) < getPreCutValue1("pt1Cut") ) continue;

       //       if(fabs(jetEtaAK4->at(0)) > getPreCutValue1("jetFidRegion") || jetPtAK4->at(0) < getPreCutValue1("pt0Cut")) continue;
       //       if(fabs(jetEtaAK4->at(1)) > getPreCutValue1("jetFidRegion") || jetPtAK4->at(1) < getPreCutValue1("pt1Cut")) continue;
       
       Ev_RecoJetOK++;
       
       // TLorentzVector jet1, jet2 leading reco jet
       jet1.SetPtEtaPhiM( (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0])) *jetPtAK4->at(sortedJetIdx[0])
			  ,jetEtaAK4->at(sortedJetIdx[0]),jetPhiAK4->at(sortedJetIdx[0])
			  , (jecFactors[sortedJetIdx[0]]/jetJecAK4->at(sortedJetIdx[0])) * jetMassAK4->at(sortedJetIdx[0]));
       jet2.SetPtEtaPhiM( (jecFactors[sortedJetIdx[1]]/jetJecAK4->at(sortedJetIdx[1])) *jetPtAK4->at(sortedJetIdx[1])
			  ,jetEtaAK4->at(sortedJetIdx[1]),jetPhiAK4->at(sortedJetIdx[1])
			  , (jecFactors[sortedJetIdx[1]]/jetJecAK4->at(sortedJetIdx[1])) * jetMassAK4->at(sortedJetIdx[1]));
       
       //       jet1.SetPtEtaPhiM(jetPtAK4->at(0), jetEtaAK4->at(0), jetPhiAK4->at(0), jetMassAK4->at(0));
       //       jet2.SetPtEtaPhiM(jetPtAK4->at(1), jetEtaAK4->at(1), jetPhiAK4->at(1), jetMassAK4->at(1));
       
       for(Long64_t ijet=0; ijet<no_jets_ak4; ijet++){ //jet loop for ak4
	 TLorentzVector currentJet;
	 
	 //	 if(fabs(jetEtaAK4->at(ijet)) < getPreCutValue1("jetFidRegion") 
	 //	    && idTAK4->at(ijet) == getPreCutValue1("tightJetID") 
	 //	    && jetPtAK4->at(ijet) > getPreCutValue1("ptCut")) {
	 if(fabs(jetEtaAK4->at(sortedJetIdx[ijet])) < getPreCutValue1("jetFidRegion") 
	    && idTAK4->at(sortedJetIdx[ijet]) == getPreCutValue1("tightJetID") 
	    && (jecFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet]))*jetPtAK4->at(sortedJetIdx[ijet]) > getPreCutValue1("ptCut"))
	   {
	     
	     TLorentzVector currentJet;
	     currentJet.SetPtEtaPhiM( (jecFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet])) *jetPtAK4->at(sortedJetIdx[ijet])
				      ,jetEtaAK4->at(sortedJetIdx[ijet]),jetPhiAK4->at(sortedJetIdx[ijet])
				      , (jecFactors[sortedJetIdx[ijet]]/jetJecAK4->at(sortedJetIdx[ijet])) *jetMassAK4->at(sortedJetIdx[ijet]));   
	     //	     currentJet.SetPtEtaPhiM(jetPtAK4->at(ijet), jetEtaAK4->at(ijet), jetPhiAK4->at(ijet), jetMassAK4->at(ijet));   
	     
	     double DeltaR1 = currentJet.DeltaR(jet1);
	     double DeltaR2 = currentJet.DeltaR(jet2);
	     
	     if(DeltaR1 < DeltaR2 && DeltaR1 < wideJetDeltaR_){
	       wj1_tmp += currentJet;
	     }
	     else if(DeltaR2 < wideJetDeltaR_){
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
       
     //+++++++++++++++++++
     // Utilities for SMEARING FUNCTIONS
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

       Ev_Matching++;

       for(int kk=0; kk<2; kk++){// loop on GenWideJets
	 
	 int etaBin = mEtaBinning.getBin( fabs(GenWideJet_Eta[kk]) );
	 int ptBin = mPtBinning.getPtBin( GenWideJet_Pt[kk] );
	 
	 if(verbose){
	   const std::pair<float, float> etaBins = mEtaBinning.getBinValue(etaBin);
	   const std::pair<float, float> ptBins = mPtBinning.getBinValue(ptBin);
	   cout<<"etaGen  "<< fabs(GenWideJet_Eta[kk]) << "    pTGen   "<<  GenWideJet_Pt[kk] <<endl;
	   cout<<"etaBin  "<< etaBin << "    pTBin   "<<ptBin<<endl;
	   cout<<"etaBin.first  "<< etaBins.first << "    etaBin.second   "<<etaBins.second<<endl;
	   cout<<"ptBin.first  "<< ptBins.first << "    ptBin.second   "<<ptBins.second<<endl;
	 }

	 //	 if(Parton_pdgId[kk] == 21) smearingFunction_G[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk]) ;	 

	 if(n_categories ==2){
	   if(Parton_pdgId[kk] == 21){
	     smearingFunction_G[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk]) ;	 
	   }else{
	     smearingFunction_Q[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   } 
	 } else if(n_categories == 3 ){
	   if(Parton_pdgId[kk] == 21){
	     smearingFunction_G[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk]) ;	 
	   }else if( fabs(Parton_pdgId[kk]) == 1 || fabs(Parton_pdgId[kk]) == 2  || fabs(Parton_pdgId[kk]) == 3 ){ 
	     smearingFunction_QLight[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 4 || fabs(Parton_pdgId[kk]) == 5 ){ 
	     smearingFunction_QHeavy[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }
	 }else if(n_categories == 4 ){
	   if(Parton_pdgId[kk] == 21){
	     smearingFunction_G[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk]) ;	 
	   }else if( fabs(Parton_pdgId[kk]) == 1 || fabs(Parton_pdgId[kk]) == 2  || fabs(Parton_pdgId[kk]) == 3 ){ 
	     smearingFunction_QLight[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 4 ){ 
	     smearingFunction_QCharm[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 5 ){ 
	     smearingFunction_QBottom[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }
	 }else if(n_categories == 6 ){
       	   if(Parton_pdgId[kk] == 21){
	     smearingFunction_G[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk]) ;	 
	   }else if( fabs(Parton_pdgId[kk]) == 1){ 
	     smearingFunction_QUp[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 2){ 
	     smearingFunction_QDown[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 3){ 
	     smearingFunction_QStrange[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 4 ){ 
	     smearingFunction_QCharm[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }else if( fabs(Parton_pdgId[kk]) == 5 ){ 
	     smearingFunction_QBottom[etaBin][ptBin]->Fill( R_PtWidejet_PtGenWideJet[kk] );
	   }
	 }  
       }

     //++++++++++++++++++++++++++++++++
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

//////////////////////////////////////////////////////////////

template<typename T>
std::vector<T*> baseClass::buildPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {
  //std::vector<T*> baseClass::buildPtVector( const std::string& branchName, int nBins, double xMin, double xMax) {

  bool appendText = (xMin >= 0 && xMax >= 0);
  std::vector<T*> vector;
  size_t ptBinningSize = mPtBinning.size();
  for (size_t j = 0; j < ptBinningSize; j++) {

    const std::pair<float, float> bin = mPtBinning.getBinValue(j);
    std::stringstream ss;
    if (appendText)
      ss << branchName << "_pt_" << (int) bin.first << "_" << (int) bin.second;
    else
      ss << branchName << "_" << (int) bin.first << "_" << (int) bin.second;

    if (!appendText) {
      xMin = bin.first;
    }

    if (!appendText) {
      xMax = bin.second;
    }

    T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
    vector.push_back(object);
  }

  return vector;
}
////////////////////////

template<typename T>
std::vector<T*> baseClass::buildPtVector(TFileDirectory dir,  const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {
  //std::vector<T*> baseClass::buildPtVector(  const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {
  return buildPtVector<T>(dir, branchName + "_" + etaName, nBins, xMin, xMax);
}
///////////////////////////////
template<typename T>
std::vector<std::vector<T*> > baseClass::buildEtaPtVector(TFileDirectory dir,  const std::string& branchName, int nBins, double xMin, double xMax) {
  //std::vector<std::vector<T*> > baseClass::buildEtaPtVector(  const std::string& branchName, int nBins, double xMin, double xMax) {

  size_t etaBinningSize = mEtaBinning.size();
  std::vector<std::vector<T*> > etaBinning;

  
  for (size_t i = 0; i < etaBinningSize; i++) {
    const std::string etaName = mEtaBinning.getBinName(i);
    etaBinning.push_back(buildPtVector<T>(dir, branchName, etaName, nBins, xMin, xMax));
  }

  return etaBinning;
}
