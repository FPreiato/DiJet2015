//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 29 17:18:52 2015 by ROOT version 5.34/18
// from TChain dijets/events/
//////////////////////////////////////////////////////////

#ifndef rootNtupleClass_h
#define rootNtupleClass_h

//// Lines added by make_rootNtupleClass.sh - BEGIN 
#include <vector> 
using namespace std; 
//// Lines added by make_rootNtupleClass.sh - END 

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxrun = 1;
const Int_t kMaxevt = 1;
const Int_t kMaxlumi = 1;
const Int_t kMaxnVtx = 1;
const Int_t kMaxrho = 1;
const Int_t kMaxmet = 1;
const Int_t kMaxmetSig = 1;
const Int_t kMaxnJetsAK4 = 1;
const Int_t kMaxhtAK4 = 1;
const Int_t kMaxmjjAK4 = 1;
const Int_t kMaxdEtajjAK4 = 1;
const Int_t kMaxdPhijjAK4 = 1;
const Int_t kMaxnJetsAK8 = 1;
const Int_t kMaxhtAK8 = 1;
const Int_t kMaxmjjAK8 = 1;
const Int_t kMaxdEtajjAK8 = 1;
const Int_t kMaxdPhijjAK8 = 1;
const Int_t kMaxptHat = 1;
const Int_t kMaxprocessID = 1;
const Int_t kMaxweight = 1;
const Int_t kMaxnGenJetsAK4 = 1;
const Int_t kMaxnGenJetsAK8 = 1;

class rootNtupleClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runNo;
   Int_t           evtNo;
   Int_t           lumi;
   Int_t           nvtx;
   Float_t         rho;
   Float_t         met;
   Float_t         metSig;
   vector<float>   *gen_eta;
   vector<float>   *gen_phi;
   vector<float>   *gen_p;
   vector<float>   *gen_px;
   vector<float>   *gen_py;
   vector<float>   *gen_pz;
   vector<float>   *gen_pt;
   vector<float>   *gen_energy;
   vector<int>     *gen_pdgId;
   vector<float>   *gen_vx;
   vector<float>   *gen_vy;
   vector<float>   *gen_vz;
   vector<int>     *gen_numDaught;
   vector<int>     *gen_status;
   vector<int>     *gen_index;
   vector<int>     *gen_motherIndex;
   Int_t           nJetsAK4;
   Float_t         htAK4;
   Float_t         mjjAK4;
   Float_t         dEtajjAK4;
   Float_t         dPhijjAK4;
   Int_t           nJetsAK8;
   Float_t         htAK8;
   Float_t         mjjAK8;
   Float_t         dEtajjAK8;
   Float_t         dPhijjAK8;
   vector<float>   *jetPtAK4;
   vector<float>   *jetJecAK4;
   vector<float>   *jetEtaAK4;
   vector<float>   *jetPhiAK4;
   vector<float>   *jetMassAK4;
   vector<float>   *jetEnergyAK4;
   vector<float>   *jetAreaAK4;
   vector<float>   *jetChfAK4;
   vector<float>   *jetNhfAK4;
   vector<float>   *jetPhfAK4;
   vector<float>   *jetMufAK4;
   vector<float>   *jetElfAK4;
   vector<int>     *idLAK4;
   vector<int>     *idTAK4;
   vector<float>   *jetPtAK8;
   vector<float>   *jetJecAK8;
   vector<float>   *jetEtaAK8;
   vector<float>   *jetPhiAK8;
   vector<float>   *jetMassAK8;
   vector<float>   *jetEnergyAK8;
   vector<float>   *jetAreaAK8;
   vector<float>   *jetChfAK8;
   vector<float>   *jetNhfAK8;
   vector<float>   *jetPhfAK8;
   vector<float>   *jetMufAK8;
   vector<float>   *jetElfAK8;
   vector<int>     *idLAK8;
   vector<int>     *idTAK8;
   vector<float>   *jetMassPrunedAK8;
   vector<float>   *jetTau1AK8;
   vector<float>   *jetTau2AK8;
   vector<float>   *jetTau3AK8;
   vector<bool>    *triggerResult;
   vector<float>   *npu;
   vector<int>     *PileupInteractions;
   vector<int>     *PileupOriginBX;
   Float_t         ptHat;
   Int_t           processID;
   Float_t         weight;
   Int_t           nGenJetsAK4;
   Int_t           nGenJetsAK8;
   vector<float>   *jetPtGenAK4;
   vector<float>   *jetEtaGenAK4;
   vector<float>   *jetPhiGenAK4;
   vector<float>   *jetMassGenAK4;
   vector<float>   *jetEnergyGenAK4;
   vector<float>   *jetPtGenAK8;
   vector<float>   *jetEtaGenAK8;
   vector<float>   *jetPhiGenAK8;
   vector<float>   *jetMassGenAK8;
   vector<float>   *jetEnergyGenAK8;

   // List of branches
   TBranch        *b_run_;   //!
   TBranch        *b_evt_;   //!
   TBranch        *b_lumi_;   //!
   TBranch        *b_nVtx_;   //!
   TBranch        *b_rho_;   //!
   TBranch        *b_met_;   //!
   TBranch        *b_metSig_;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_p;   //!
   TBranch        *b_gen_px;   //!
   TBranch        *b_gen_py;   //!
   TBranch        *b_gen_pz;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_energy;   //!
   TBranch        *b_gen_pdgId;   //!
   TBranch        *b_gen_vx;   //!
   TBranch        *b_gen_vy;   //!
   TBranch        *b_gen_vz;   //!
   TBranch        *b_gen_numDaught;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_index;   //!
   TBranch        *b_gen_motherIndex;   //!
   TBranch        *b_nJetsAK4_;   //!
   TBranch        *b_htAK4_;   //!
   TBranch        *b_mjjAK4_;   //!
   TBranch        *b_dEtajjAK4_;   //!
   TBranch        *b_dPhijjAK4_;   //!
   TBranch        *b_nJetsAK8_;   //!
   TBranch        *b_htAK8_;   //!
   TBranch        *b_mjjAK8_;   //!
   TBranch        *b_dEtajjAK8_;   //!
   TBranch        *b_dPhijjAK8_;   //!
   TBranch        *b_jetPtAK4;   //!
   TBranch        *b_jetJecAK4;   //!
   TBranch        *b_jetEtaAK4;   //!
   TBranch        *b_jetPhiAK4;   //!
   TBranch        *b_jetMassAK4;   //!
   TBranch        *b_jetEnergyAK4;   //!
   TBranch        *b_jetAreaAK4;   //!
   TBranch        *b_jetChfAK4;   //!
   TBranch        *b_jetNhfAK4;   //!
   TBranch        *b_jetPhfAK4;   //!
   TBranch        *b_jetMufAK4;   //!
   TBranch        *b_jetElfAK4;   //!
   TBranch        *b_idLAK4;   //!
   TBranch        *b_idTAK4;   //!
   TBranch        *b_jetPtAK8;   //!
   TBranch        *b_jetJecAK8;   //!
   TBranch        *b_jetEtaAK8;   //!
   TBranch        *b_jetPhiAK8;   //!
   TBranch        *b_jetMassAK8;   //!
   TBranch        *b_jetEnergyAK8;   //!
   TBranch        *b_jetAreaAK8;   //!
   TBranch        *b_jetChfAK8;   //!
   TBranch        *b_jetNhfAK8;   //!
   TBranch        *b_jetPhfAK8;   //!
   TBranch        *b_jetMufAK8;   //!
   TBranch        *b_jetElfAK8;   //!
   TBranch        *b_idLAK8;   //!
   TBranch        *b_idTAK8;   //!
   TBranch        *b_jetMassPrunedAK8;   //!
   TBranch        *b_jetTau1AK8;   //!
   TBranch        *b_jetTau2AK8;   //!
   TBranch        *b_jetTau3AK8;   //!
   TBranch        *b_triggerResult;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_PileupInteractions;   //!
   TBranch        *b_PileupOriginBX;   //!
   TBranch        *b_ptHat_;   //!
   TBranch        *b_processID_;   //!
   TBranch        *b_weight_;   //!
   TBranch        *b_nGenJetsAK4_;   //!
   TBranch        *b_nGenJetsAK8_;   //!
   TBranch        *b_jetPtGenAK4;   //!
   TBranch        *b_jetEtaGenAK4;   //!
   TBranch        *b_jetPhiGenAK4;   //!
   TBranch        *b_jetMassGenAK4;   //!
   TBranch        *b_jetEnergyGenAK4;   //!
   TBranch        *b_jetPtGenAK8;   //!
   TBranch        *b_jetEtaGenAK8;   //!
   TBranch        *b_jetPhiGenAK8;   //!
   TBranch        *b_jetMassGenAK8;   //!
   TBranch        *b_jetEnergyGenAK8;   //!

   rootNtupleClass(TTree *tree=0);
   virtual ~rootNtupleClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef rootNtupleClass_cxx
rootNtupleClass::rootNtupleClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("dijets/events",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("dijets/events","");
      chain->Add("../ntuple/RSGravitonToQuarkQuark_kMpl01_M_3000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM_10_1_WRL.root/dijets/events");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

rootNtupleClass::~rootNtupleClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t rootNtupleClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t rootNtupleClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void rootNtupleClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   gen_eta = 0;
   gen_phi = 0;
   gen_p = 0;
   gen_px = 0;
   gen_py = 0;
   gen_pz = 0;
   gen_pt = 0;
   gen_energy = 0;
   gen_pdgId = 0;
   gen_vx = 0;
   gen_vy = 0;
   gen_vz = 0;
   gen_numDaught = 0;
   gen_status = 0;
   gen_index = 0;
   gen_motherIndex = 0;
   jetPtAK4 = 0;
   jetJecAK4 = 0;
   jetEtaAK4 = 0;
   jetPhiAK4 = 0;
   jetMassAK4 = 0;
   jetEnergyAK4 = 0;
   jetAreaAK4 = 0;
   jetChfAK4 = 0;
   jetNhfAK4 = 0;
   jetPhfAK4 = 0;
   jetMufAK4 = 0;
   jetElfAK4 = 0;
   idLAK4 = 0;
   idTAK4 = 0;
   jetPtAK8 = 0;
   jetJecAK8 = 0;
   jetEtaAK8 = 0;
   jetPhiAK8 = 0;
   jetMassAK8 = 0;
   jetEnergyAK8 = 0;
   jetAreaAK8 = 0;
   jetChfAK8 = 0;
   jetNhfAK8 = 0;
   jetPhfAK8 = 0;
   jetMufAK8 = 0;
   jetElfAK8 = 0;
   idLAK8 = 0;
   idTAK8 = 0;
   jetMassPrunedAK8 = 0;
   jetTau1AK8 = 0;
   jetTau2AK8 = 0;
   jetTau3AK8 = 0;
   triggerResult = 0;
   npu = 0;
   PileupInteractions = 0;
   PileupOriginBX = 0;
   jetPtGenAK4 = 0;
   jetEtaGenAK4 = 0;
   jetPhiGenAK4 = 0;
   jetMassGenAK4 = 0;
   jetEnergyGenAK4 = 0;
   jetPtGenAK8 = 0;
   jetEtaGenAK8 = 0;
   jetPhiGenAK8 = 0;
   jetMassGenAK8 = 0;
   jetEnergyGenAK8 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNo", &runNo, &b_run_);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evt_);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi_);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nVtx_);
   fChain->SetBranchAddress("rho", &rho, &b_rho_);
   fChain->SetBranchAddress("met", &met, &b_met_);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig_);
   fChain->SetBranchAddress("gen_eta", &gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_p", &gen_p, &b_gen_p);
   fChain->SetBranchAddress("gen_px", &gen_px, &b_gen_px);
   fChain->SetBranchAddress("gen_py", &gen_py, &b_gen_py);
   fChain->SetBranchAddress("gen_pz", &gen_pz, &b_gen_pz);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_energy", &gen_energy, &b_gen_energy);
   fChain->SetBranchAddress("gen_pdgId", &gen_pdgId, &b_gen_pdgId);
   fChain->SetBranchAddress("gen_vx", &gen_vx, &b_gen_vx);
   fChain->SetBranchAddress("gen_vy", &gen_vy, &b_gen_vy);
   fChain->SetBranchAddress("gen_vz", &gen_vz, &b_gen_vz);
   fChain->SetBranchAddress("gen_numDaught", &gen_numDaught, &b_gen_numDaught);
   fChain->SetBranchAddress("gen_status", &gen_status, &b_gen_status);
   fChain->SetBranchAddress("gen_index", &gen_index, &b_gen_index);
   fChain->SetBranchAddress("gen_motherIndex", &gen_motherIndex, &b_gen_motherIndex);
   fChain->SetBranchAddress("nJetsAK4", &nJetsAK4, &b_nJetsAK4_);
   fChain->SetBranchAddress("htAK4", &htAK4, &b_htAK4_);
   fChain->SetBranchAddress("mjjAK4", &mjjAK4, &b_mjjAK4_);
   fChain->SetBranchAddress("dEtajjAK4", &dEtajjAK4, &b_dEtajjAK4_);
   fChain->SetBranchAddress("dPhijjAK4", &dPhijjAK4, &b_dPhijjAK4_);
   fChain->SetBranchAddress("nJetsAK8", &nJetsAK8, &b_nJetsAK8_);
   fChain->SetBranchAddress("htAK8", &htAK8, &b_htAK8_);
   fChain->SetBranchAddress("mjjAK8", &mjjAK8, &b_mjjAK8_);
   fChain->SetBranchAddress("dEtajjAK8", &dEtajjAK8, &b_dEtajjAK8_);
   fChain->SetBranchAddress("dPhijjAK8", &dPhijjAK8, &b_dPhijjAK8_);
   fChain->SetBranchAddress("jetPtAK4", &jetPtAK4, &b_jetPtAK4);
   fChain->SetBranchAddress("jetJecAK4", &jetJecAK4, &b_jetJecAK4);
   fChain->SetBranchAddress("jetEtaAK4", &jetEtaAK4, &b_jetEtaAK4);
   fChain->SetBranchAddress("jetPhiAK4", &jetPhiAK4, &b_jetPhiAK4);
   fChain->SetBranchAddress("jetMassAK4", &jetMassAK4, &b_jetMassAK4);
   fChain->SetBranchAddress("jetEnergyAK4", &jetEnergyAK4, &b_jetEnergyAK4);
   fChain->SetBranchAddress("jetAreaAK4", &jetAreaAK4, &b_jetAreaAK4);
   fChain->SetBranchAddress("jetChfAK4", &jetChfAK4, &b_jetChfAK4);
   fChain->SetBranchAddress("jetNhfAK4", &jetNhfAK4, &b_jetNhfAK4);
   fChain->SetBranchAddress("jetPhfAK4", &jetPhfAK4, &b_jetPhfAK4);
   fChain->SetBranchAddress("jetMufAK4", &jetMufAK4, &b_jetMufAK4);
   fChain->SetBranchAddress("jetElfAK4", &jetElfAK4, &b_jetElfAK4);
   fChain->SetBranchAddress("idLAK4", &idLAK4, &b_idLAK4);
   fChain->SetBranchAddress("idTAK4", &idTAK4, &b_idTAK4);
   fChain->SetBranchAddress("jetPtAK8", &jetPtAK8, &b_jetPtAK8);
   fChain->SetBranchAddress("jetJecAK8", &jetJecAK8, &b_jetJecAK8);
   fChain->SetBranchAddress("jetEtaAK8", &jetEtaAK8, &b_jetEtaAK8);
   fChain->SetBranchAddress("jetPhiAK8", &jetPhiAK8, &b_jetPhiAK8);
   fChain->SetBranchAddress("jetMassAK8", &jetMassAK8, &b_jetMassAK8);
   fChain->SetBranchAddress("jetEnergyAK8", &jetEnergyAK8, &b_jetEnergyAK8);
   fChain->SetBranchAddress("jetAreaAK8", &jetAreaAK8, &b_jetAreaAK8);
   fChain->SetBranchAddress("jetChfAK8", &jetChfAK8, &b_jetChfAK8);
   fChain->SetBranchAddress("jetNhfAK8", &jetNhfAK8, &b_jetNhfAK8);
   fChain->SetBranchAddress("jetPhfAK8", &jetPhfAK8, &b_jetPhfAK8);
   fChain->SetBranchAddress("jetMufAK8", &jetMufAK8, &b_jetMufAK8);
   fChain->SetBranchAddress("jetElfAK8", &jetElfAK8, &b_jetElfAK8);
   fChain->SetBranchAddress("idLAK8", &idLAK8, &b_idLAK8);
   fChain->SetBranchAddress("idTAK8", &idTAK8, &b_idTAK8);
   fChain->SetBranchAddress("jetMassPrunedAK8", &jetMassPrunedAK8, &b_jetMassPrunedAK8);
   fChain->SetBranchAddress("jetTau1AK8", &jetTau1AK8, &b_jetTau1AK8);
   fChain->SetBranchAddress("jetTau2AK8", &jetTau2AK8, &b_jetTau2AK8);
   fChain->SetBranchAddress("jetTau3AK8", &jetTau3AK8, &b_jetTau3AK8);
   fChain->SetBranchAddress("triggerResult", &triggerResult, &b_triggerResult);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("PileupInteractions", &PileupInteractions, &b_PileupInteractions);
   fChain->SetBranchAddress("PileupOriginBX", &PileupOriginBX, &b_PileupOriginBX);
   fChain->SetBranchAddress("ptHat", &ptHat, &b_ptHat_);
   fChain->SetBranchAddress("processID", &processID, &b_processID_);
   fChain->SetBranchAddress("weight", &weight, &b_weight_);
   fChain->SetBranchAddress("nGenJetsAK4", &nGenJetsAK4, &b_nGenJetsAK4_);
   fChain->SetBranchAddress("nGenJetsAK8", &nGenJetsAK8, &b_nGenJetsAK8_);
   fChain->SetBranchAddress("jetPtGenAK4", &jetPtGenAK4, &b_jetPtGenAK4);
   fChain->SetBranchAddress("jetEtaGenAK4", &jetEtaGenAK4, &b_jetEtaGenAK4);
   fChain->SetBranchAddress("jetPhiGenAK4", &jetPhiGenAK4, &b_jetPhiGenAK4);
   fChain->SetBranchAddress("jetMassGenAK4", &jetMassGenAK4, &b_jetMassGenAK4);
   fChain->SetBranchAddress("jetEnergyGenAK4", &jetEnergyGenAK4, &b_jetEnergyGenAK4);
   fChain->SetBranchAddress("jetPtGenAK8", &jetPtGenAK8, &b_jetPtGenAK8);
   fChain->SetBranchAddress("jetEtaGenAK8", &jetEtaGenAK8, &b_jetEtaGenAK8);
   fChain->SetBranchAddress("jetPhiGenAK8", &jetPhiGenAK8, &b_jetPhiGenAK8);
   fChain->SetBranchAddress("jetMassGenAK8", &jetMassGenAK8, &b_jetMassGenAK8);
   fChain->SetBranchAddress("jetEnergyGenAK8", &jetEnergyGenAK8, &b_jetEnergyGenAK8);
   Notify();
}

Bool_t rootNtupleClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void rootNtupleClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t rootNtupleClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef rootNtupleClass_cxx
