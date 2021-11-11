//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Aug  1 11:02:43 2021 by ROOT version 6.18/04
// from TTree evt/a Tree with SDHCAL frame storage
// found on file: /data/root_trees/Run_1032.root
//////////////////////////////////////////////////////////

#ifndef ana_FEBv2_h
#define ana_FEBv2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ana_FEBv2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   Int_t	   hv_;
   Int_t	   sn_;
   Int_t	   mt_;
//   Int_t	   muW1_;
//   Int_t	   muW2_;
   const char *    loc_;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          event;
   UInt_t          bc0;
   UInt_t          nframe;
   ULong64_t       frame[65535];   //[nframe]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bc0;   //!
   TBranch        *b_nframe;   //!
   TBranch        *b_frame;   //!

   ana_FEBv2(Int_t hv, Int_t sn , Int_t mt, /*Int_t muW1, Int_t muW2 ,*/ const char * loc , TTree *tree=0);
   virtual ~ana_FEBv2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ana_FEBv2_cxx
ana_FEBv2::ana_FEBv2(Int_t hv, Int_t sn , Int_t mt, /*Int_t muW1, Int_t muW2,*/ const char * loc , TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   hv_ = hv;
   sn_ = sn;
   mt_ = mt;
 //  muW1_ = muW1;
 //  muW2_ = muW2;
   loc_ = loc;
   TString 	   s("");
   s.Form("%s_HV_%d_SN_%d_MaxTrig_%d_Run1176.root",loc,hv,sn,mt);
   std::cout<<"Formatting the string"<<std::endl;
   if (tree == 0) {
	   TChain *chain = new TChain("evt"); 
	   chain->Add(s);
           //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(s);
           //if (!f || !f->IsOpen()) {
           //   f = new TFile(s);
           //}
           //f->GetObject("evt",tree);
   	   tree = chain;
   }
   Init(tree);
}

ana_FEBv2::~ana_FEBv2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ana_FEBv2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ana_FEBv2::LoadTree(Long64_t entry)
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

void ana_FEBv2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bc0", &bc0, &b_bc0);
   fChain->SetBranchAddress("nframe", &nframe, &b_nframe);
   fChain->SetBranchAddress("frame", frame, &b_frame);
   Notify();
}

Bool_t ana_FEBv2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ana_FEBv2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ana_FEBv2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ana_FEBv2_cxx
