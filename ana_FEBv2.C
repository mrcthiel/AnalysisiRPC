#define ana_FEBv2_cxx
#include "ana_FEBv2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
using namespace std;

struct stripTimes
{
	void addLRframe(int i , float time, int f,  int p) {LRframe.push_back(i); LRtime.push_back(time) ; LRframeOK.push_back(true); LRframeBkg.push_back(false); LRfpga_num.push_back(f); LRpetiroc_num.push_back(p);}
	void addHRframe(int i, float time, int f, int p) {HRframe.push_back(i); HRtime.push_back(time) ; HRframeOK.push_back(true); HRframeBkg.push_back(false); HRfpga_num.push_back(f);  HRpetiroc_num.push_back(p);}
	std::vector<int> LRframe;
	std::vector<float> LRtime;
	std::vector<bool> LRframeOK;
	std::vector<int> HRframe;
	std::vector<int> LRfpga_num;
	std::vector<int> HRfpga_num;
	std::vector<int> LRpetiroc_num;
	std::vector<int> HRpetiroc_num;
	std::vector<float> HRtime;
	std::vector<bool> HRframeOK;
	std::vector<bool> HRframeBkg;
	std::vector<bool> LRframeBkg;

	std::vector<std::pair<int,int> > deltas;
	void buildDeltas()
	{
		deltas.clear();
		for (int iLR=0; iLR<LRframe.size();++iLR)
		{
			if ( not (LRframeOK[iLR])) continue;
			for (int iHR=0; iHR<HRframe.size(); ++iHR)
			{
				if (HRframeOK[iHR]) deltas.push_back(std::make_pair(LRframe[iLR],HRframe[iHR]));
			}
		}
	}
};


void ana_FEBv2::Loop()
{
	if (fChain == 0) return;
	TString         s("");
	s.Form("will save plots for RUN: , _HV_%d_SN_%d_MaxTrig_%d ",hv_,sn_,mt_);

	std::cout<<s<<std::endl;
	//std::ifstream infile("corr_1167.txt");
	std::ifstream infile("corr_1216.txt");
	string Load;
	float corr_LR[50] , corr_HR[50];
	float corr_lr=0, corr_hr=0;
	int i=0;
	while (getline(infile,Load))
	{
		istringstream is(Load);
		is >> corr_lr >> corr_hr;
		corr_LR[i] = corr_lr;
		corr_HR[i] = corr_hr;
		i++;
	}
	std::cout << corr_LR[0] << " " << corr_HR[0] << std::endl;
	//std::cout << "Muon Window:" << muW1_ << " " << muW2_ << std::endl;
	//std::cout << "Bkg Window:" << muW1_-4000-(9*(abs(muW1_)-abs(muW2_))) << " " << muW2_-4000 << std::endl;
	//double time_corr[3] = {12.39 ,12.24,12.5 };
	double time_corr[3] = {0 ,0,0 };
	double time_min = 4000;
	double time_max = 87000;


	//double time_corr_fine[3][32] = {
	//{0.00, -0.37, 1.06, 5.59, 0.59, 3.83, 0.15, 2.84, 2.19, 0.95, 2.52, 3.89, 9.14, 2.62, 4.05, 7.11, 3.48, 8.20, 7.89, 2.97, 7.50, 4.08, 4.79, 3.10, 4.83, 2.84, 8.87, 7.16, -0.01, 1.85, 3.49, 0.57},
	//{0.00, -0.47, 1.11, 5.64, 0.57, 3.80, 0.13, 2.86, 2.23, 0.82, 2.61, 4.06, 9.23, 2.59, 4.10, 7.06, 3.32, 8.13, 8.62, 2.60, 7.67, 3.87, 4.64, 3.05, 4.64, 2.58, 8.72, 7.11, -0.15, 1.62, 3.45, 0.49},
	//{0.00, -0.53, 1.22, 5.74, 0.67, 4.06, 0.19, 2.88, 2.30, 1.08, 2.63, 4.15, 9.32, 2.78, 4.32, 7.33, 3.81, 8.44, 8.03, 3.14, 7.87, 4.32, 5.20, 3.42, 4.99, 3.03, 9.22, 7.59, 0.19, 2.05, 3.76, 0.91}
	//};

	double time_corr_fine[3][32] = {
		{15.3147,13.3131,13.2824,15.2839,14.2824,15.3131,12.3131,15.3131,11.3131,14.3131,15.3131,16.3131,19.074,15.3131,17.3131,19.3131,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19,-1148.19},
		{15.3131,13.3131,13.3131,16.2528,14.3131,15.3131,12.3131,15.3131,11.3131,14.3131,15.3131,16.3131,19.3131,15.3131,18.3127,19.3131,17.3129,18.3131,18.3131,15.3131,18.3131,17.4229,17.303,16.3131,17.3131,17.3131,21.3131,19.3131,17.3131,16.3131,15.3131,15.3131},
		{15.3131,13.3131,13.3131,15.3131,14.3131,15.3131,12.3131,15.3131,11.3838,14.3131,15.3131,16.3131,18.3131,15.3131,16.324,19.0838,16.3131,18.3131,18.3131,15.3131,17.3132,17.3131,16.3131,16.3131,17.3105,17.3131,20.3089,18.3093,17.3089,15.3092,14.9573,15.3093}
	};

	uint64_t nentries;
	nentries = fChain->GetEntriesFast();
	std::cout<<"N Events: " << nentries << std::endl;
	//Long64_t nentries = 100;

	Long64_t nbytes = 0, nb = 0;
	uint64_t t_bc0[3];
	bool debug=false;
	float y = 0;
	float last_bc0 = 0;
	float run_duration_in_sec = 0;
	int any_strip_fired=0;
	int trig_exists=0;
	int ntrig_all=0;
	float deltaT = -999;
	int ntrig_allevent=0;
	int ntriggerHR_signal=0;
	int ntriggerLR_signal=0;
	int ntriggerHR3ns_signal=0;
	int ntriggerLR3ns_signal=0;
	Bool_t histo_exists=0;
	int bad_trig=0; 
	int nbad_trig=0; 
	int ntrig_event=0; 
	int Rstrip_fired=0; 
	int Dstrip_fired=0; 
	int Rstrip_bkg_fired=0; 
	int Dstrip_bkg_fired=0; 
	int nRstrip_fired=0; 
	int nDstrip_fired=0; 
	int nRstrip_bkg_fired=0; 
	int nORstrip_fired=0; 
	int nANDstrip_fired=0; 
	int nDstrip_bkg_fired=0; 
	int has_paired_strips=0;
	int nANDstrip_tight_fired=0; 
	int nANDstrip_medium_fired=0; 
	int tight_num=0; 
	int medium_num=0;
	int nANDstrip_noCut_fired=0; 
	Long64_t trig_index = 0;
	float return_SR_int = 0;
	float return_BR_int = 0;
	float direct_SR_int = 0;
	float direct_BR_int = 0;
	float eff_medium = 0;
	int n_paired_srip = 0;

	// The timing of the first hit in the muon window. Considered as original hit
	int first_hit[6];
	float first_hit_time[6];

	TH1F* hLR = new TH1F("hLR", "LR", 50, 0.0, 50);
	TH1F* hHR = new TH1F("hHR", "HR", 50, 0.0, 50);
	//Adding new histograms
	TH1F* hLRn = new TH1F("hLRn", "hLRn", 50, 0.0, 50);
	TH1F* hHRn = new TH1F("hHRn", "hHRn", 50, 0.0, 50);
	//
	TH2F* hdeltaT = new TH2F("deltaT", "T (HR - LR)", 50,0,50,600, -30, 30);
	TH1F* hdeltaT_1D = new TH1F("deltaT", "T (HR - LR)", 600, -30, 30);
	TH2F* hsumT = new TH2F("sumT", "T (HR + LR)", 50,0,50,200, 0, 200);
	TH1F* hHRT_Xtalk = new TH1F("hHRT_Xtalk", "Time diff at HR for Xtalk", 35, -25.0, 10);
//	TH2F* hHRT = new TH2F("hHRT", "T (HR - Trig)", 50,0,50,abs(muW1_-muW2_), muW1_, muW2_);
//	TH2F* hLRT = new TH2F("hLRT", "T (LR - Trig)", 50,0,50,abs(muW1_-muW2_), muW1_, muW2_);
	TH1F* hLRT_ = new TH1F("hLRT_","T (LR - Trig)" ,10000,-9000,1000);
	TH1F* hHRT_ = new TH1F("hHRT_","T (HR - Trig)" ,10000,-9000,1000);

        TH1F* hLRT_temp = new TH1F("hLRT_muon_window","T (LR - Trig)" ,10000,-9000,1000);
        TH1F* hHRT_temp = new TH1F("hHRT_muon_window","T (HR - Trig)" ,10000,-9000,1000);

	TH1F* hLRT_2 = new TH1F("hLRT_2","n strip LR" ,20,0,20);
	TH1F* hHRT_2 = new TH1F("hHRT_2","n strip HR" ,20,0,20);
	TH1F* hLRT_3 = new TH1F("hLRT_3","n strip LR (LR - HR > 3 ns)" ,20,0,20);
	TH1F* hHRT_3 = new TH1F("hHRT_3","n strip HR (LR - HR > 3 ns)" ,20,0,20);
	TH1F* hnPairs = new TH1F("hnPairs","n paired strips" ,20,0,20);
	//s.Form("plots_904/events/_HV_%d_SN_%d_MaxTrig_%d_.root",hv_,sn_,mt_);

	//TFile *outputFile = new TFile(s,"RECREATE");
	///Get Final bc0 to calculate overall noise
	fChain->Draw("bc0>>histbc0");
	TH1F *hbc0 = (TH1F*)gDirectory->Get("histbc0");
	std::cout << "============>>>>>>>>>>>>>>> hbc0->Integral(): " << hbc0->Integral() << std::endl;

	std::cout << "============>>>>>>>>>>>>>>> hbc0->GetXaxis()->GetXmax(): " << hbc0->GetXaxis()->GetXmax() << std::endl;

	last_bc0 = hbc0->GetXaxis()->GetXmax();
	std::cout<<"integral/3 : "<<hbc0->Integral()/3<<" " << last_bc0 <<" "<<hbc0->FindLastBinAbove(1)<<std::endl;
//	run_duration_in_sec = 10*(abs(muW1_)-abs(muW2_))*pow(10,-9);
	//run_duration_in_sec = ((hbc0->Integral()/3)*90*pow(10,-6));
	TH1F *htrig= new TH1F("htrig", "trigger time",91000,0,91000);
	fChain->Draw("m_time(m_traw(frame))>>htrig","(m_channel(frame)==33)");
	std::cout<<"All triggers /3 : "<< (htrig->Integral())/3 << std::endl;
	//run_duration_in_sec = ((htrig->Integral())/3 )*90*pow(10,-6); 

	std::cout<<"The data taking time is : "<< run_duration_in_sec <<std::endl;
	//fChain->Draw("(1)>>histTrigger","(m_channel(frame)==33)");
	//TH1F *histTrig = (TH1F*)gDirectory->Get("histTrigger");
	//TH2F *hHREvt[100]; //For histograms for each event
	//TH2F *hLREvt[100];
/*	for (uint32_t i=0;i<100;i++) {
		hHREvt[i] = new TH2F(s, "T (return - Trig)", 50,0,50,abs(muW1_-muW2_),muW1_,muW2_);
		hLREvt[i] = new TH2F(s, "T (direct - Trig)", 50,0,50,abs(muW1_-muW2_),muW1_,muW2_);
	}
*/	ofstream myfile;
	s.Form("outputs/_HV_%d_SN_%d_MaxTrig_%d_.txt",hv_,sn_,mt_);
	myfile.open(s);

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		stripTimes allStrips_temp[48];

		trig_exists=0;
		bad_trig = 0;
		ntrig_event = 0;
                int trig[3];
                std::fill_n(trig, 3, 0);

		for (uint32_t i=0;i<nframe;i++)
		{
			if (m_channel(frame[i])==32) {
				t_bc0[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
			}
			else if (m_channel(frame[i])==33){
				if (trig[m_fpga(frame[i])] != 0) bad_trig = 1; // prevents firing of one of the fpga's more than once
				trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
				if ((trig[m_fpga(frame[i])]) < time_max and trig[m_fpga(frame[i])] > time_min) {
					trig_exists=1; // trigger exits in this event
					ntrig_event++; // This value must reach 3 for all events.
				}
			}
			else {
				any_strip_fired = 1; // In this event there is at least one stripped fired. Without any requirement.
				if (c_side(m_channel(frame[i])) > 0.5){
					allStrips_temp[m_strip(frame[i])].addLRframe(i,m_time(m_traw(frame[i]))-(time_corr[m_fpga(frame[i])]+time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]),m_fpga(frame[i]),int((c_petiroc(m_channel(frame[i]))/16)) );
				}
				if (c_side(m_channel(frame[i])) < 0.5){
					allStrips_temp[m_strip(frame[i])].addHRframe(i, m_time(m_traw(frame[i]))-(time_corr[m_fpga(frame[i])]+time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]),m_fpga(frame[i]),int((c_petiroc(m_channel(frame[i]))/16)));
				}
			}
		} // end of first frame loop
		if (!(trig_exists && ntrig_event ==3) ) continue;
		if (bad_trig) continue;

		for (uint32_t i=0;i<48;i++) {

			for (uint32_t fhr=0; fhr<allStrips_temp[i].HRframe.size(); fhr++){
				hHRT_temp->Fill(allStrips_temp[i].HRtime[fhr]-trig[m_fpga(frame[allStrips_temp[i].HRframe[fhr]])]);
			}
			for (uint32_t flr=0; flr<allStrips_temp[i].LRframe.size(); flr++) { // frame LR loop
				hLRT_temp->Fill(allStrips_temp[i].LRtime[flr]-trig[m_fpga(frame[allStrips_temp[i].LRframe[flr]])]);
			}
		}
	}

	hHRT_temp->GetXaxis()->SetRangeUser(hHRT_temp->GetXaxis()->GetBinCenter(hHRT_temp->GetMaximumBin())-50,hHRT_temp->GetXaxis()->GetBinCenter(hHRT_temp->GetMaximumBin())+50);
        hLRT_temp->GetXaxis()->SetRangeUser(hLRT_temp->GetXaxis()->GetBinCenter(hLRT_temp->GetMaximumBin())-50,hLRT_temp->GetXaxis()->GetBinCenter(hLRT_temp->GetMaximumBin())+50);

	hHRT_temp->Fit("gaus");
	hLRT_temp->Fit("gaus");

	std::cout << "hHRT_temp->GetFunction(\"gaus\")->GetParameter(1): " << hHRT_temp->GetFunction("gaus")->GetParameter(1) << std::endl;
	std::cout << "hLRT_temp->GetFunction(\"gaus\")->GetParameter(1): " << hLRT_temp->GetFunction("gaus")->GetParameter(1) << std::endl;

	double muW1_HR = hHRT_temp->GetFunction("gaus")->GetParameter(1)-20.;
        double muW2_HR = hHRT_temp->GetFunction("gaus")->GetParameter(1)+20.;
        double muW1_LR = hLRT_temp->GetFunction("gaus")->GetParameter(1)-20.;
        double muW2_LR = hLRT_temp->GetFunction("gaus")->GetParameter(1)+20.;

        double muW1_bkg = hHRT_temp->GetFunction("gaus")->GetParameter(1)-20.-1000.;
        double muW2_bkg = hHRT_temp->GetFunction("gaus")->GetParameter(1)+20.-1000.;


        TH2F* hHRT = new TH2F("hHRT", "T (HR - Trig)", 50,0,50,abs(muW1_HR-muW2_HR), muW1_HR, muW2_HR);
        TH2F* hLRT = new TH2F("hLRT", "T (LR - Trig)", 50,0,50,abs(muW1_LR-muW2_LR), muW1_LR, muW2_LR);
        run_duration_in_sec = 10*(abs(muW1_HR)-abs(muW2_HR))*pow(10,-9);
        //for (uint32_t i=0;i<100;i++) {
        //        hHREvt[i] = new TH2F(s, "T (return - Trig)", 50,0,50,abs(muW1_HR-muW2_HR),muW1_HR,muW2_HR);
        //        hLREvt[i] = new TH2F(s, "T (direct - Trig)", 50,0,50,abs(muW1_LR-muW2_LR),muW1_LR,muW2_LR);
        //}


	for (int jj=0;jj<2;jj++){
	if(jj==1){
                muW1_HR = muW1_bkg;
                muW2_HR = muW2_bkg;
                muW1_LR = muW1_bkg;
                muW2_LR = muW2_bkg;
	}

	ntrig_allevent=0;
	nANDstrip_fired=0;
	nANDstrip_noCut_fired=0;
	nANDstrip_medium_fired=0;
	nANDstrip_tight_fired=0;

      
	
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		//std::cout<<"Hello "<<jentry<<" bc0 "<<bc0<<std::endl;
		// if (Cut(ientry) < 0) continue;
		// float delta[3][100];
		//std::vector<float> delta[48];
		deltaT = -999;
		float sum[3][100];
		std::vector<float> time_hr[3][100];
		std::vector<float> time_lr[3][100];
		int trig[3];
		std::fill_n(trig, 3, 0);
		for (int i =0; i<3;i++) {
			std::fill_n(sum[i], 100, -999999);
		}
		std::vector<int> Rstrip_fired_f[48];

		stripTimes allStrips[48];


		any_strip_fired=0;
		trig_exists=0;
		bad_trig = 0;
		ntrig_event = 0;
		Rstrip_fired=0; 
		Dstrip_fired=0; 
		Rstrip_bkg_fired=0; 
		Dstrip_bkg_fired=0; 
		ntriggerHR_signal=0;
		ntriggerLR_signal=0;
		ntriggerHR3ns_signal=0;
		ntriggerLR3ns_signal=0;
		has_paired_strips=0;
		n_paired_srip=0;
		medium_num=0;
		tight_num=0;

		std::fill_n(first_hit, 6, -1);
		std::fill_n(first_hit_time, 6, 9999999);


		for (uint32_t i=0;i<nframe;i++)
		{
			if (m_channel(frame[i])==32) {
				t_bc0[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
			}
			else if (m_channel(frame[i])==33){
				if (trig[m_fpga(frame[i])] != 0) bad_trig = 1; // prevents firing of one of the fpga's more than once
				trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
				if ((trig[m_fpga(frame[i])]) < time_max and trig[m_fpga(frame[i])] > time_min) {
					//std::cout<<"fpga: "<< m_fpga(frame[i]) << "Trig Exists !!!!! " << std::endl;
					trig_exists=1; // trigger exits in this event
					ntrig_event++; // This value must reach 3 for all events.
				}
				if (debug) myfile<<"\t "<<"trig time "<< m_time(m_traw(frame[i])) <<"\n";
			}
			else {
				any_strip_fired = 1; // In this event there is at least one stripped fired. Without any requirement.
				if (c_side(m_channel(frame[i])) > 0.5){
					//LR stips
					allStrips[m_strip(frame[i])].addLRframe(i,m_time(m_traw(frame[i]))-(time_corr[m_fpga(frame[i])]+time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]),m_fpga(frame[i]), int((c_petiroc(m_channel(frame[i]))/16)) );
					if(jj==0) hLR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min)*1e-9*120)); //For Noise
				}
				if (c_side(m_channel(frame[i])) < 0.5){
					//HR strips
					allStrips[m_strip(frame[i])].addHRframe(i, m_time(m_traw(frame[i]))-(time_corr[m_fpga(frame[i])]+time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]),m_fpga(frame[i]),int((c_petiroc(m_channel(frame[i]))/16)) );
					if(jj==0) hHR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min)*1e-9*120)); //For noise
				}
			}
		} // end of first frame loop

		//if (ntrig_event!=3 and ntrig_event!=0) std::cout<< "nTRIG in this event : "<< ntrig_event  <<std::endl;
		if (!(trig_exists && ntrig_event ==3) ) continue;
		//if (!(trig_exists) ) continue;
		if (bad_trig) continue;





		for (uint32_t i=0;i<48;i++) {

			//std::cout<< "LR: "<< allStrips[i].LRframe.size()<<" HR: "<<allStrips[i].HRframe.size() << std::endl;

			for (uint32_t fhr=0; fhr<allStrips[i].HRframe.size(); fhr++)
			{ // frame HR loop
				if(jj==0) hHRT->Fill(i,allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]);
				if(jj==0) hHRT_->Fill(allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]);

				if ( !((allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) > muW1_HR and (allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) < muW2_HR ))
				{
					// HR not in muon window
					allStrips[i].HRframeOK[fhr] = false;
					if ( ((allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) > muW1_HR-1000-(9*(abs(muW1_HR)-abs(muW2_HR))) and (allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) < muW2_HR-1000 ))
						//trig = 2000 and 87000, signal is 1000,         
					{
						allStrips[i].HRframeBkg[fhr] = true;
						if(jj==0) hLRn->Fill(m_strip(frame[i]),1/(run_duration_in_sec*120)); //For Noise; newly added histogram
					}
				}
				if (allStrips[i].HRframeOK[fhr]) {
				  ntriggerHR_signal++;
				  int asic = allStrips[i].HRfpga_num[fhr]*2+allStrips[i].HRpetiroc_num[fhr];
				  double time_corr = allStrips[i].HRtime[fhr];

				  if (first_hit_time[asic] > time_corr ){
				    first_hit_time[asic] = time_corr;
				    first_hit[asic] = m_channel(frame[allStrips[i].HRframe[fhr]]);
				    cout << first_hit_time[asic] << " asic " << asic << " hit " << first_hit[asic] << endl;
				  }
				}

			}
			for (uint32_t flr=0; flr<allStrips[i].LRframe.size(); flr++) { // frame LR loop

				if(jj==0) hLRT->Fill(i,allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]);
				if(jj==0) hLRT_->Fill(allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]);
				if ( !((allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) > muW1_LR and (allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) < muW2_LR ))
				{
					// LR not in muon window
					allStrips[i].LRframeOK[flr] = false;
					if ( ((allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) > muW1_LR-1000-(9*(abs(muW1_LR)-abs(muW2_LR))) and (allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) < muW2_LR-1000 ))
					{ 
						allStrips[i].LRframeBkg[flr] = true;
						if(jj==0) hHRn->Fill(m_strip(frame[i]),1/(run_duration_in_sec*120)); //For Noise; newly added histogram
					}   //THIS AND THE PREVIOUS FOUR LINES ADDED BY ME
				}
				if (allStrips[i].LRframeOK[flr]) ntriggerLR_signal++;
			} // LR Frame loop ends
			allStrips[i].buildDeltas();
			//std::cout<< "how many delta computed : "<<allStrips[i].deltas.size() << std::endl;
			if (allStrips[i].deltas.size() > 0 )
			{
				has_paired_strips = 1; //In the event there is at least 1 paired strip
				n_paired_srip++;
				
				
			}
			for (uint32_t d=0; d<allStrips[i].deltas.size(); d++)
			{

			  int asic = m_fpga(frame[allStrips[i].deltas[d].first]) + m_channel(frame[allStrips[i].deltas[d].first])/16;
			  cout << "fpga is " << m_fpga(frame[allStrips[i].deltas[d].first]) << "asic is" << asic << " channel is "<<  m_channel(frame[allStrips[i].deltas[d].first]) << endl;
			  cout << "asic lowest time " << first_hit_time[asic] << endl;


			  double hr_hit_time_corr = m_time(m_traw(frame[allStrips[i].deltas[d].first]))
			    +time_corr_fine[m_fpga(frame[allStrips[i].deltas[d].first])][m_channel(frame[allStrips[i].deltas[d].first])];
			
			  
			  if (first_hit_time[asic]< 999999) hHRT_Xtalk->Fill(first_hit_time[asic]-hr_hit_time_corr);


			  
				deltaT=(
						(
						 m_time(m_traw(frame[allStrips[i].deltas[d].first]))
						 -
						 (time_corr[m_fpga(frame[allStrips[i].deltas[d].first])]+time_corr_fine[m_fpga(frame[allStrips[i].deltas[d].first])][m_channel(frame[allStrips[i].deltas[d].first])])
						)
						-
						(
						 m_time(m_traw(frame[allStrips[i].deltas[d].second]))
						 -
						 (time_corr[m_fpga(frame[allStrips[i].deltas[d].second])]+time_corr_fine[m_fpga(frame[allStrips[i].deltas[d].second])][m_channel(frame[allStrips[i].deltas[d].second])])
						)
				       );
				//std::cout << "LR : " << allStrips[i].deltas[d].first << " HR : " << allStrips[i].deltas[d].second << " delta T " << deltaT << std::endl;
				if(jj==0) hdeltaT->Fill(i,deltaT);
				if(jj==0) hdeltaT_1D->Fill(deltaT);
				if (deltaT > 3 ) {medium_num = 1; }
				if (deltaT > 8 ) {tight_num = 1; }
			}
		} // strip loop ends

		if((ntriggerHR_signal>0) || (ntriggerLR_signal>0)) nANDstrip_noCut_fired++;

		ntrig_allevent++; // count denom of eff
		if (has_paired_strips) {
			nANDstrip_fired++;
			if(jj==0) hHRT_2->Fill(ntriggerHR_signal);
			if(jj==0) hLRT_2->Fill(ntriggerLR_signal);
		}
		if (medium_num) {nANDstrip_medium_fired++;}
		if (tight_num) {nANDstrip_tight_fired++;}
		if(jj==0) hnPairs->Fill(n_paired_srip);
		//std::cout << "ntrig_allevent: " << ntrig_allevent << std::endl;
	} //Event loop ends

	//outputFile->Close();
	myfile<<"Efficiency AND no cut: " << (nANDstrip_noCut_fired)/float(ntrig_allevent) <<"\n";
	myfile<<"Efficiency AND Loose: " << (nANDstrip_fired)/float(ntrig_allevent) <<"\n";
	myfile<<"Efficiency AND Medium: " << (nANDstrip_medium_fired)/float(ntrig_allevent) <<"\n";
	//eff_medium = (nANDstrip_medium_fired)/float(ntrig_allevent);
	myfile<<"Efficiency AND Tight: " << (nANDstrip_tight_fired)/float(ntrig_allevent) <<"\n";
	//myfile<<"Efficiency AND Medium err : " <<  sqrt((1-eff_medium)*eff_medium/float(ntrig_allevent)) <<"\n";
	//myfile<<"hHRn->Integral()/ntrig_allevent/48: " << hHRn->Integral()/ntrig_allevent/48 <<"\n";
	//myfile<<"hLRn->Integral()/ntrig_allevent/48: " << hLRn->Integral()/ntrig_allevent/48 <<"\n";
	//myfile.close();
        myfile<<"===================================="<<"\n";
	}

        myfile.close();

	gStyle->SetOptStat(0);
	TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",200,10,780,780);
	TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.05,0.98,0.55,21);
	TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.55,0.98,0.98,21);
	pad1->UseCurrentStyle(); pad2->UseCurrentStyle();
	pad1->Draw(); pad2->Draw();
	pad1->cd();
	pad1->SetGrid();
	pad1->SetLogy();
	//hLR->Scale(1/ntrig_allevent);
	hLR->Draw("Hist");
	hLR->GetYaxis()->SetTitle("Noise Hz/cm");
	hLR->GetXaxis()->SetTitle("Strips");
	pad2->cd();
	pad2->SetGrid();
	pad2->SetLogy();
	hHR->Draw("Hist");
	hHR->GetYaxis()->SetTitle("Noise Hz/cm");
	hHR->GetXaxis()->SetTitle("Strips");
	c1->cd();
	c1->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__avarage_noise_side.root",hv_,sn_,mt_);
	c1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__avarage_noise_side.pdf",hv_,sn_,mt_);
	c1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__avarage_noise_side.png",hv_,sn_,mt_);
	c1->SaveAs(s);

	//The newly added noise histograms go here
	gStyle->SetOptStat(0);
	TCanvas *can1 = new TCanvas("can1","The Ntuple canvas",200,10,780,780);
	TPad *newpad1 = new TPad("newpad1","This is newpad1",0.05,0.05,0.98,0.55,21);
	TPad *newpad2 = new TPad("newpad2","This is newpad2",0.05,0.55,0.98,0.98,21);
	newpad1->UseCurrentStyle(); newpad2->UseCurrentStyle();
	newpad1->Draw(); newpad2->Draw();
	newpad1->cd();
	newpad1->SetGrid();
	hHRT_Xtalk->Draw();
	newpad2->cd();
	newpad2->SetGrid();
	//newpad2->SetLogy();
	hHRn->Scale(1./ntrig_allevent);
	hHRn->Draw("Hist");
	hHRn->GetYaxis()->SetTitle("Noise Hz/cm");
	hHRn->GetXaxis()->SetTitle("Strips");
	can1->cd();
	can1->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__new_noise_histo.root",hv_,sn_,mt_);
	can1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__new_noise_histo.pdf",hv_,sn_,mt_);
	can1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__new_noise_histo.png",hv_,sn_,mt_);
	can1->SaveAs(s);

	TCanvas *c2 = new TCanvas("c4","The Ntuple canvas",200,10,900,780);
	TPad *pad4 = new TPad("pad4","This is pad4",0.02,0.02,0.98,0.98,21);
	pad4->UseCurrentStyle();
	pad4->Draw();
	pad4->cd();
	pad4->SetGrid();
	pad4->SetRightMargin(0.12);
	hdeltaT->GetYaxis()->SetTitle("delta time [ns]");
	hdeltaT->GetXaxis()->SetTitle("Strips");
	hdeltaT->GetYaxis()->SetLabelOffset(0.007);
	hdeltaT->GetYaxis()->SetTitleOffset(1.5);
	hdeltaT->SetLabelSize(0.02,"x");
	hdeltaT->SetTitle("");
	hdeltaT->Draw("P");
	c2->cd();
	c2->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__deltaT_srip.root",hv_,sn_,mt_);
	c2->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__deltaT_srip.pdf",hv_,sn_,mt_);
	c2->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__deltaT_srip.png",hv_,sn_,mt_);
	c2->SaveAs(s);
	TCanvas *cPairs = new TCanvas("cPairs","The Ntuple canvas",200,10,900,780);
	TPad *padP = new TPad("padP","This is padP",0.02,0.02,0.98,0.98,21);
	padP->UseCurrentStyle();
	padP->Draw();
	padP->cd();
	padP->SetGrid();
	padP->SetRightMargin(0.12);
	hnPairs->GetYaxis()->SetTitle("Events");
	hnPairs->GetXaxis()->SetTitle("n paired strips");
	hnPairs->GetYaxis()->SetLabelOffset(0.007);
	hnPairs->GetYaxis()->SetTitleOffset(1.5);
	hnPairs->SetLabelSize(0.02,"x");
	hnPairs->SetTitle("");
	hnPairs->Draw("Histo");
	cPairs->cd();
	cPairs->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__n_paired_srip.root",hv_,sn_,mt_);
	cPairs->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__n_paired_srip.pdf",hv_,sn_,mt_);
	cPairs->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__n_paired_srip.png",hv_,sn_,mt_);
	cPairs->SaveAs(s);    
	TCanvas *chdeltaT_1D = new TCanvas("chdeltaT_1D","The Ntuple canvas",200,10,900,780);
	TPad *padhdeltaT_1D = new TPad("padhdeltaT_1D","This is hdeltaT_1D",0.02,0.02,0.98,0.98,21);
	padhdeltaT_1D->UseCurrentStyle();
	padhdeltaT_1D->Draw();
	padhdeltaT_1D->cd();
	padhdeltaT_1D->SetGrid();
	padhdeltaT_1D->SetRightMargin(0.12);
	hdeltaT_1D->GetXaxis()->SetTitle("delta time [ns]");
	hdeltaT_1D->GetYaxis()->SetTitle("Events");
	hdeltaT_1D->GetYaxis()->SetLabelOffset(0.007);
	hdeltaT_1D->GetYaxis()->SetTitleOffset(1.5);
	hdeltaT_1D->SetLabelSize(0.02,"x");
	hdeltaT_1D->SetTitle("");
	hdeltaT_1D->Scale(1/(hdeltaT_1D->Integral(hdeltaT_1D->FindBin(9),hdeltaT_1D->FindBin(20))));
	//hdeltaT_1D->Scale(1/(hdeltaT_1D->Integral()));
	hdeltaT_1D->Draw("Histo");
	chdeltaT_1D->cd();
	chdeltaT_1D->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__deltaT_1D_srip.root",hv_,sn_,mt_);
	chdeltaT_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__deltaT_1D_srip.pdf",hv_,sn_,mt_);
	chdeltaT_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__deltaT_1D_srip.png",hv_,sn_,mt_);
	chdeltaT_1D->SaveAs(s);

	TCanvas *c3 = new TCanvas("c3","The Ntuple canvas",200,10,780,780);
	TPad *pad3 = new TPad("pad3","This is pad3",0.05,0.05,0.98,0.55,21);
	TPad *pad5 = new TPad("pad5","This is pad5",0.05,0.55,0.98,0.98,21);
	pad3->UseCurrentStyle(); pad5->UseCurrentStyle();
	pad3->Draw(); pad5->Draw();
	pad3->cd();
	pad3->SetGrid();
	hLRT->Draw("COLZ");
	pad5->cd();
	pad5->SetGrid();
	hHRT->Draw("COLZ");
	c3->cd();
	c3->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side.root",hv_,sn_,mt_);
	c3->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side.pdf",hv_,sn_,mt_);
	c3->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side.png",hv_,sn_,mt_);
	c3->SaveAs(s);


        TCanvas *ct3 = new TCanvas("ct3","The Ntuple canvas",200,10,780,780);
        TPad *padt3 = new TPad("padt3","This is padt3",0.05,0.05,0.98,0.55,21);
        TPad *padt5 = new TPad("padt5","This is padt5",0.05,0.55,0.98,0.98,21);
        padt3->UseCurrentStyle(); padt5->UseCurrentStyle();
        padt3->Draw(); padt5->Draw();
        padt3->cd();
        hLRT_temp->Draw();
        padt5->cd();
        hHRT_temp->Draw();
        ct3->cd();
        ct3->Update();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side_forTimeWindow.png",hv_,sn_,mt_);
        ct3->SaveAs(s);



	TCanvas *c10 = new TCanvas("c10","The Ntuple canvas",200,10,780,780);
	TPad *pad10 = new TPad("pad10","This is pad10",0.05,0.05,0.98,0.55,21);
	TPad *pad11 = new TPad("pad11","This is pad11",0.05,0.55,0.98,0.98,21);
	pad10->UseCurrentStyle(); pad11->UseCurrentStyle();
	pad10->Draw(); pad11->Draw();
	pad10->cd();
	pad10->SetGrid();
	hLRT_->Draw("Hist");
	pad11->cd();
	pad11->SetGrid();
	hHRT_->Draw("Hist");
	c10->cd();
	c10->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side_1d.root",hv_,sn_,mt_);
	c10->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side_1d.pdf",hv_,sn_,mt_);
	c10->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side_1d.png",hv_,sn_,mt_);
	c10->SaveAs(s);

	TCanvas *c11 = new TCanvas("c11","The Ntuple canvas",200,10,780,780);
	TPad *pad20 = new TPad("pad20","This is pad10",0.05,0.05,0.98,0.55,21);
	TPad *pad21 = new TPad("pad21","This is pad11",0.05,0.55,0.98,0.98,21);
	pad20->UseCurrentStyle(); pad21->UseCurrentStyle();
	pad20->Draw(); pad21->Draw();
	pad20->cd();
	pad20->SetGrid();
	hLRT_2->Draw("Hist");
	pad21->cd();
	pad21->SetGrid();
	hHRT_2->Draw("Hist");
	c11->cd();
	c11->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__nStrip_side_1d.root",hv_,sn_,mt_);
	c11->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__nStrip_side_1d.pdf",hv_,sn_,mt_);
	c11->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_%d__nStrip_side_1d.png",hv_,sn_,mt_);
	c11->SaveAs(s);


	gApplication->Terminate();

}
