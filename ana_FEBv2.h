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
		//   Int_t	   mt_;
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

		ana_FEBv2(Int_t hv, Int_t sn , /*Int_t mt, Int_t muW1, Int_t muW2 ,*/ const char * loc, TTree *tree=0);
		virtual ~ana_FEBv2();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop();
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		double active_area = 10577.5/2; //re41
		bool isFEBv2r2 = true;
		bool plot_clust_par = false;
		bool print_gamma_histos = false;
		bool no_trig_case = false;
		bool align_XY = false;
		bool retrig_study = false;
	        double V_conector = 299.792*0.568;//0.6;// mm/ns  - check it

                TH1F* hLR_ns = new TH1F("hLR_ns", "number of LR signals", 3500, 0.0, 35000);
                TH1F* hHR_ns = new TH1F("hHR_ns", "number of HR signals", 3500, 0.0, 35000);
		TH1F* hLR = new TH1F("hLR", "LR", 50, 0.0, 50);
		TH1F* hHR = new TH1F("hHR", "HR", 50, 0.0, 50);
		TH1F* hLRn = new TH1F("hLRn", "hLRn", 50, 0.0, 50);
		TH1F* hHRn = new TH1F("hHRn", "hHRn", 50, 0.0, 50);
                TH1F* hLHRn = new TH1F("hLHRn", "hLHRn", 50, 0.0, 50);
                TH1F* hLHRn_dif = new TH1F("hLHRn_dif", "hLHRn_dif", 50, 0.0, 50);
                TH1F* hLHR_dif = new TH1F("hLHR_dif", "LHR_dif", 50, 0.0, 50);
		TH2F* nFiredHR_per_strip = new TH2F("nFiredHR_per_strip", "# fired strip ; strip; #", 50,0,50,10, 1, 11);
		TH2F* nFiredLR_per_strip = new TH2F("nFiredLR_per_strip", "# fired strip ; strip; #", 50,0,50,10, 1, 11);
		TH2F* hdeltaT = new TH2F("deltaT", "T (HR - LR)", 50,0,50,300, -30, 30);
		TH2F* hsumT = new TH2F("sumT", "T (HR + LR)", 50,0,50,1000, -5000, 0);
		TH1F* hdeltaT_1D = new TH1F("deltaT", "T (HR - LR)", 300, -30, 30);
		TH1F* hsumT_1D = new TH1F("sumT", "T (HR + LR)", 1000, -5000, 0);
		TH1F* hLRT_ = new TH1F("hLRT_","T (LR - Trig)" ,10000,-9000,1000);
		TH1F* hHRT_ = new TH1F("hHRT_","T (HR - Trig)" ,10000,-9000,1000);
                TH1F* hLRT_strip23 = new TH1F("hLRT_strip23","T (LR - Trig)" ,10000,-9000,1000);
                TH1F* hHRT_strip23 = new TH1F("hHRT_strip23","T (HR - Trig)" ,10000,-9000,1000);
		TH1F* hLRT_temp = new TH1F("hLRT_muon_window","T (LR - Trig)" ,10000,-9000,1000);
		TH1F* hHRT_temp = new TH1F("hHRT_muon_window","T (HR - Trig)" ,10000,-9000,1000);
		TH1F* hLRT_2 = new TH1F("hLRT_2","n strip LR" ,20,0,20);
		TH1F* hHRT_2 = new TH1F("hHRT_2","n strip HR" ,20,0,20);
		TH1F* hnPairs = new TH1F("hnPairs","n paired strips" ,20,0,20);
		TH1F* hdeltaTCluster_1D = new TH1F("hdeltaTCluster_1D", "Cluster T (HR - LR)", 300, -30, 30);
		TH1F* hsumTCluster_1D = new TH1F("hsumTCluster_1D", "Cluster T (HR + LR)", 1000, -5000, 0);
		TH1F* hClusterSize = new TH1F("hClusterSize","Cluster Size" ,10,0,10);
		TH1F* hNClusters = new TH1F("hNClusters","Number of Clusters" ,10,0,10);
		TH2F* hdeltaClusterT = new TH2F("hdeltaClusterT", "ClusterT (HR - LR)", 50,0,50,300, -30, 30);
		TH1F* strip_diff_h = new TH1F("strip_diff_h", " HR-LR strip cluster difference ", 19, -5, 5);
		TGraph *g_cluster_size = new TGraph();
		TGraph *g_cluster_number = new TGraph();

//		TH2F* h2_XY = new TH2F("h2_XY"," X - Y ",400,-60,720,400,-20,1300);
                TH2F* h2_XY = new TH2F("h2_XY"," X - Y ",2000,-100,760,6000,-200,1520);
                TH2F* h2_XY_cls = new TH2F("h2_XY_cls"," X - Y ",2000,-100,760,6000,-200,1520);


		Long64_t nbytes = 0, nb = 0;
		uint64_t t_bc0[3];
		float last_bc0 = 0;
		float run_duration_in_sec = 0;
		int any_strip_fired=0;
		int trig_exists=0;
		int ntrig_all=0;
		float deltaT = -999;
		float sumT = 0;
		int ntrig_allevent=0;
		int ntrig_allevent_muon_window=0;
                int ntrig_allevent_gamma_window=0;
		int ntriggerHR_signal=0;
		int ntriggerLR_signal=0;
		int bad_trig=0;
		int ntrig_event=0;
		int nANDstrip_fired=0;
		int has_paired_strips=0;
		int nANDstrip_noCut_fired=0;
		int nANDstripHR_noCut_fired=0;
		int nANDstripLR_noCut_fired=0;
		int N_cluster_good=0;
		int n_paired_srip = 0;
		int n_paired_srip_med = 0;

};

#endif

#ifdef ana_FEBv2_cxx
ana_FEBv2::ana_FEBv2(Int_t hv, Int_t sn , /*Int_t mt, Int_t muW1, Int_t muW2,*/ const char * loc, TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	hv_ = hv;
	sn_ = sn;

	//  mt_ = mt;
	//  muW1_ = muW1;
	//  muW2_ = muW2;
	loc_ = loc;
	TString 	   s("");
	s.Form("%s_HV_%d_SN_%d_MaxTrig*.root",loc,hv,sn/*,mt*/);
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
