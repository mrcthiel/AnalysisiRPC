#define ana_FEBv2_cxx
#include "ana_FEBv2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "MuonWindow.C"
#include "stripTimes.h"
#include "CalibAlig.h"
#include "Cluster_Algo.h"
#include "plot_and_save.C"

using namespace std;

void ana_FEBv2::Loop()
{

	if (fChain == 0) return;
	TString         s("");
	s.Form("will save plots for RUN: , _HV_%d_SN_%d_MaxTrig_ ",hv_,sn_/*,mt_*/);
	std::cout<<s<<std::endl;

	double time_min = 4000; //To avoid the edges of the bc0
	double time_max = 87000; //To avoid the edges of the bc0

	// for retrig study
	s.Form("mkdir ScanId_%d/HV%d/plts_per_strip",sn_,hv_);
        gSystem->Exec(s);
        if(retrig_study){ 
		s.Form("mkdir ScanId_%d/HV%d/retriggering_study",sn_,hv_);
		gSystem->Exec(s);
	}
	std::vector<TH1F*> hLRT_V; 
        std::vector<TH1F*> hHRT_V;
        std::vector<TH1F*> hLRT_zoom_V;
        std::vector<TH1F*> hHRT_zoom_V;
        std::vector<TH1F*> hdeltaT_1D_V;

        std::vector<TH1F*> hLRT_retrig_V;
        std::vector<TH1F*> hHRT_retrig_V;

	for (int i=0;i<48;i++){ 
		TH1F *h1 = new TH1F(Form("hLRT_V_%d",i),"T (LR - Trig)",10000,-9000,1000); 
		hLRT_V.push_back(h1); 
                TH1F *h2 = new TH1F(Form("hHRT_V_%d",i),"T (HR - Trig)",10000,-9000,1000);
                hHRT_V.push_back(h2);
                TH1F *h3 = new TH1F(Form("hLRT_zoom_V_%d",i),"T (LR - Trig)",600,-2000,1000);
                hLRT_zoom_V.push_back(h3);
                TH1F *h4 = new TH1F(Form("hHRT_zoom_V_%d",i),"T (HR - Trig)",600,-2000,1000);
                hHRT_zoom_V.push_back(h4);
                TH1F *h5 = new TH1F(Form("deltaT_V_%d",i),"T (HR - LR)", 300, -30, 30);
                hdeltaT_1D_V.push_back(h5);
                TH1F *h6 = new TH1F(Form("hLRT_retrig_V_%d",i),"Tn = LRn - Trig; Tn-T1, n>1",150,0,150);
                hLRT_retrig_V.push_back(h6);
                TH1F *h7 = new TH1F(Form("hHRT_retrig_V_%d",i),"Tn = HRn - Trig; Tn-T1, n>1",150,0,150);
                hHRT_retrig_V.push_back(h7);
	}


	// Get time calibration for a given ScanId
	vector<vector<double>> time_corr_fine_vec;
	time_corr_fine_vec.clear();
	time_corr_fine_vec =  time_corr_fine_func(sn_);
	double time_corr_fine[3][32];
	for(int ii=0;ii<time_corr_fine_vec.at(0).size();ii++){
		time_corr_fine[0][ii] = time_corr_fine_vec.at(0).at(ii);				
	}
	for(int ii=0;ii<time_corr_fine_vec.at(1).size();ii++){
		time_corr_fine[1][ii] = time_corr_fine_vec.at(1).at(ii);
	}
	for(int ii=0;ii<time_corr_fine_vec.at(2).size();ii++){
		time_corr_fine[2][ii] = time_corr_fine_vec.at(2).at(ii);
	}

	// Get strip length for aligment
	vector<double> Strip_length_vec;
	Strip_length_vec.clear();
	vector<double> LR_to_conctor_vec;
	LR_to_conctor_vec.clear();
	vector<double> HR_to_conctor_vec;
	HR_to_conctor_vec.clear();
	HR_to_conctor_vec = Aligment_factor(sn_,1);
	LR_to_conctor_vec = Aligment_factor(sn_,2);
	Strip_length_vec = Aligment_factor(sn_,3);
	double HR_to_conctor[48];
	for(int ii=0;ii<HR_to_conctor_vec.size();ii++){
		HR_to_conctor[ii] = HR_to_conctor_vec.at(ii);
	}
	double LR_to_conctor[48];
	for(int ii=0;ii<LR_to_conctor_vec.size();ii++){
		LR_to_conctor[ii] = LR_to_conctor_vec.at(ii);
	}
	double Strip_length[48];
	for(int ii=0;ii<Strip_length_vec.size();ii++){
		Strip_length[ii] = Strip_length_vec.at(ii);
	}
	int reference_strip = 0;
	for(int iref = 0; iref < 48; iref++){
		if(Strip_length[iref] < Strip_length[reference_strip]) reference_strip = iref; //the smallest one
	}


	uint64_t nentries = fChain->GetEntriesFast();
	std::cout<<"N Events: " << nentries << std::endl;

	///Get Final bc0 to calculate overall noise
	fChain->Draw("bc0>>histbc0");
	TH1F *hbc0 = (TH1F*)gDirectory->Get("histbc0");
	std::cout << "============>>>>>>>>>>>>>>> hbc0->Integral(): " << hbc0->Integral() << std::endl;
	std::cout << "============>>>>>>>>>>>>>>> hbc0->GetXaxis()->GetXmax(): " << hbc0->GetXaxis()->GetXmax() << std::endl;
	last_bc0 = hbc0->GetXaxis()->GetXmax();//hbc0->GetEntries();//hbc0->GetXaxis()->GetXmax();
	std::cout<<"integral/3 : "<<hbc0->Integral()/3<<" " << last_bc0 <<" "<<hbc0->FindLastBinAbove(1)<<std::endl;
	std::cout<<"hbc0->GetMean(): "<<hbc0->GetMean()<<"; hbc0->GetRMS(): "<<hbc0->GetRMS()<<"; hbc0->GetEntries(): "<<hbc0->GetEntries()<<std::endl;
	TH1F *htrig= new TH1F("htrig", "trigger time",91000,0,91000);
	fChain->Draw("m_time(m_traw(frame))>>htrig","(m_channel(frame)==33)");
	std::cout<<"All triggers /3 : "<< (htrig->Integral())/3 << std::endl;

	// Open file to save the effs
	ofstream myfile;
	s.Form("outputs/_HV_%d_SN_%d_MaxTrig_.txt",hv_,sn_/*,mt_*/);
	myfile.open(s);

	// Get the Muon Window
	auto newtree = fChain->CloneTree();
	vector<double> mu_wind_vec = MuonWindow(sn_, newtree, time_min, time_max, nentries, isFEBv2r2);
	double muon_window_size = 50.;
	if(no_trig_case==true) muon_window_size = 1000.;
	double muW1_HR = mu_wind_vec.at(1)-muon_window_size/2;
	double muW2_HR = mu_wind_vec.at(1)+muon_window_size/2;
	double muW1_LR = mu_wind_vec.at(0)-muon_window_size/2;
	double muW2_LR = mu_wind_vec.at(0)+muon_window_size/2;
	double muW1_bkg = muW1_HR-5000.;//muW1_HR-1000.
	double muW2_bkg = muW2_HR-1000.;//muW2_HR-1000.
	TH2F* hHRT = new TH2F("hHRT", "T (HR - Trig)", 50,0,50,abs(muW1_HR-muW2_HR), muW1_HR, muW2_HR);
	TH2F* hLRT = new TH2F("hLRT", "T (LR - Trig)", 50,0,50,abs(muW1_LR-muW2_LR), muW1_LR, muW2_LR);

        for (int strip_numb=0;strip_numb<48;strip_numb++){
		hHRT_zoom_V.at(strip_numb)->GetXaxis()->SetRangeUser(mu_wind_vec.at(1)-150,mu_wind_vec.at(1)+150);
		hLRT_zoom_V.at(strip_numb)->GetXaxis()->SetRangeUser(mu_wind_vec.at(0)-150,mu_wind_vec.at(0)+150);
	}
	run_duration_in_sec = 10*(abs(muW1_HR)-abs(muW2_HR))*pow(10,-9);

	// Use it to print the final effs and cluster rate
	float or_eff_on = 0;
	float and_eff_on = 0;
	float HR_eff_on = 0;
	float LR_eff_on = 0;
	float cluster_eff_on = 0;
	float or_eff_out = 0;
	float and_eff_out = 0;
	float HR_eff_out = 0;
	float LR_eff_out = 0;
	float cluster_eff_out = 0;
	float cluster_rate_out = 0;

	double hHRn_gamma_AND  = 0.;
	double hHRn_gamma_HR  = 0.;
	double hHRn_gamma_LR  = 0.;
	double gamma_cluster_size = -1.;

	int side = Feb_Chamber(sn_).at(2);

	//for Y alignment
	double Y_align[48] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        double Y_align_n[48] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        double dT_align[48] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        double dT_align_n[48] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


	for (int jj=0;jj<2;jj++){ // jj=0->on muon window; jj=1->out muon window
		double sum_cluster_size = 0;
		if(jj==1){
			muW1_HR = muW1_bkg;
			muW2_HR = muW2_bkg;
			muW1_LR = muW1_bkg;
			muW2_LR = muW2_bkg;
		}

		// to evaluate the effs
		ntrig_allevent=0; //denominator
		nANDstrip_fired=0; //AND eff
		nANDstrip_noCut_fired=0; //OR eff
		nANDstripHR_noCut_fired=0; //HR eff
		nANDstripLR_noCut_fired=0; //LR eff
		N_cluster_good=0; // for cluster rate

		//for cluster study
		vector<double> ClusterN;
		vector<double> ClusterS;
		vector<double> ClusterL;
		ClusterN.clear();
		ClusterS.clear();
		ClusterL.clear();
		int ClusterE = 0;

		//start loop over all events
                int retrig_rate = 0; //for retrig study
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
			deltaT = -999;
			int trig[3];
			std::fill_n(trig, 3, 0);

			// struct to save the signals per strip
			stripTimes allStrips[48];

			trig_exists=0;
			bad_trig = 0;
			ntrig_event = 0;
			ntriggerHR_signal=0;
			ntriggerLR_signal=0;
			has_paired_strips=0;
			n_paired_srip=0;
			n_paired_srip_med=0;
			double trig_time = -99999.;

			//loop over all signals, no triggering, here we fill rate plots inportant at P5! here we also save the trigger info
			int LR_ns=0;
			int HR_ns=0;
			int signals_notrig=0;


        		double first_trig_HR[48] = {10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.};
                        double first_trig_LR[48] = {10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.,10000.};


			for (uint32_t i=0;i<nframe;i++){
				if (m_channel(frame[i])==32) {
					t_bc0[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
				}
				else if (m_channel(frame[i])==33){
					if (trig[m_fpga(frame[i])] != 0) bad_trig = 1; // prevents firing of one of the fpga's more than once
					trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
					if ((trig[m_fpga(frame[i])]) < time_max and trig[m_fpga(frame[i])] > time_min) {
						trig_exists=1; // trigger exits in this event
						ntrig_event++; // This value must reach 3 for all events.
						trig_time = m_time(m_traw(frame[i])); //trig[m_fpga(frame[i])];trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]
					}
				} else if(m_channel(frame[i])!=32 && m_channel(frame[i])!=33){
					std::vector<int> strip_and_side = m_strip_FEBv2r2(frame[i],side);
					int strip_side = -9999;
					int strip_numb = -9999;
					if(strip_and_side.size()>0){
						strip_side = strip_and_side.at(1);
						strip_numb = strip_and_side.at(0);	
					} else continue;
					signals_notrig++;
					if ((!isFEBv2r2 && c_side(m_channel(frame[i]))>0.5) || (isFEBv2r2 && strip_side==1)){
						double calibrated_time = m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])];
						double time_in_strip = 0;
						double time_in_strip_new = 0;
						double time_to_conector_new = 0;
						if(isFEBv2r2) {
							time_in_strip = calibrated_time - trig_time - LR_to_conctor[strip_numb]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[strip_numb];
							time_to_conector_new = LR_to_conctor[strip_numb]/V_conector-LR_to_conctor[reference_strip]/V_conector;
						} else {
							time_in_strip = calibrated_time - trig_time - LR_to_conctor[m_strip(frame[i])]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[m_strip(frame[i])];
							time_to_conector_new = LR_to_conctor[m_strip(frame[i])]/V_conector-LR_to_conctor[reference_strip]/V_conector;
						}
						double LR_trig = time_in_strip_new + time_to_conector_new;
						if ((LR_trig) < muW1_LR || (LR_trig) > muW2_LR ){
							if(isFEBv2r2) {
								LR_ns++;
								hLR->Fill(strip_numb,1/(last_bc0*(time_max-time_min-(muW2_LR-muW1_LR))*1e-9*0.9*Strip_length[strip_numb]*0.1)); //For Noise
							} else {
								hLR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min-(muW2_LR-muW1_LR))*1e-9*0.9*Strip_length[m_strip(frame[i])]*0.1)); //For Noise
							}
						}
					}
					if ((!isFEBv2r2 && c_side(m_channel(frame[i]))<0.5) || (isFEBv2r2 && strip_side==0)){
						double calibrated_time = m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])];
						double time_in_strip = 0;
						double time_in_strip_new = 0;
						double time_to_conector_new = 0;
						if(isFEBv2r2) {
							time_in_strip = calibrated_time - trig_time - HR_to_conctor[strip_numb]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[strip_numb];
							time_to_conector_new = HR_to_conctor[strip_numb]/V_conector-HR_to_conctor[reference_strip]/V_conector;
						} else {
							time_in_strip = calibrated_time - trig_time - HR_to_conctor[m_strip(frame[i])]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[m_strip(frame[i])];
							time_to_conector_new = HR_to_conctor[m_strip(frame[i])]/V_conector-HR_to_conctor[reference_strip]/V_conector;
						}
						double HR_trig = time_in_strip_new + time_to_conector_new;
						if ((HR_trig) < muW1_HR || (HR_trig) > muW2_HR ){
							if(isFEBv2r2) {
								HR_ns++;
								hHR->Fill(strip_numb,1/(last_bc0*(time_max-time_min-(muW2_HR-muW1_HR))*1e-9*0.9*Strip_length[strip_numb]*0.1)); //For Noise
							} else {
								hHR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min-(muW2_HR-muW1_HR))*1e-9*0.9*Strip_length[m_strip(frame[i])]*0.1)); //For Noise
							}
						} 
					}
				}
			}
			hLR_ns->Fill(LR_ns);
			hHR_ns->Fill(HR_ns);
			//Check Tigger
			if(no_trig_case==true) {
				trig_time=50000.;
				trig_exists=true;
				ntrig_event=3;
				bad_trig=false;
			}
			if(trig_time==-99999.) continue;
			if (!(trig_exists && ntrig_event ==3) ) continue;
			if (bad_trig) continue; //checar aqui

			int retrig_rate_temp = 0; //for retrig study
			
			//loop over all events signals in the triggered event
			for (uint32_t i=0;i<nframe;i++){
				std::vector<int> strip_and_side = m_strip_FEBv2r2(frame[i],side);
				int strip_side = -9999;
				int strip_numb = -9999;
				if(strip_and_side.size()>0){
					strip_side = strip_and_side.at(1);
					strip_numb = strip_and_side.at(0);
				}
				if (m_channel(frame[i])!=33 && m_channel(frame[i])!=32){
					if ((!isFEBv2r2 && c_side(m_channel(frame[i]))>0.5) || (isFEBv2r2 && strip_side==1)){ //LR 
						double calibrated_time = m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])];
						double time_in_strip = 0;
						double time_in_strip_new = 0;
						double time_to_conector_new = 0;
						// align and calibrate the signal
						if(isFEBv2r2) {
							time_in_strip = calibrated_time - trig_time - LR_to_conctor[strip_numb]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[strip_numb];
							time_to_conector_new = LR_to_conctor[strip_numb]/V_conector-LR_to_conctor[reference_strip]/V_conector;
						} else {
							time_in_strip = calibrated_time - trig_time - LR_to_conctor[m_strip(frame[i])]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[m_strip(frame[i])];
							time_to_conector_new = LR_to_conctor[m_strip(frame[i])]/V_conector-LR_to_conctor[reference_strip]/V_conector;
						} 

						double LR_trig = time_in_strip_new + time_to_conector_new;
						//std::cout << LR_trig << ";    " << calibrated_time << ";    " << trig_time << std::endl;
						if(jj==0) hLRT_->Fill(LR_trig); // Fill the signal - trig plot // checar se esta duplicado
						if(jj==0) hLRT_temp->Fill(LR_trig); // same but zoom
						if(jj==0) hLRT_V.at(strip_numb)->Fill(LR_trig);
                                                if(jj==0) hLRT_zoom_V.at(strip_numb)->Fill(LR_trig);
						if(isFEBv2r2 && strip_numb==23)hLRT_strip23->Fill(LR_trig);
                                                if(jj==0 && retrig_study){
                                                        if(first_trig_LR[strip_numb]==10000.) {first_trig_LR[strip_numb] = LR_trig;} else {
                                                                hLRT_retrig_V.at(strip_numb)->Fill(LR_trig-first_trig_LR[strip_numb]);
								retrig_rate_temp++;
                                                        }

                                                }
						// add signal to the strip struct
						if(isFEBv2r2) {
							allStrips[strip_numb].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
						} else {
							allStrips[m_strip(frame[i])].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
						}
					}
					if ((!isFEBv2r2 && c_side(m_channel(frame[i]))<0.5) || (isFEBv2r2 && strip_side==0)){ //HR
						double calibrated_time = m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])];
						double time_in_strip = 0;
						double time_in_strip_new = 0;
						double time_to_conector_new = 0;
						// align and calibrate the signal
						if(isFEBv2r2) {
							time_in_strip = calibrated_time - trig_time - HR_to_conctor[strip_numb]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[strip_numb];
							time_to_conector_new = HR_to_conctor[strip_numb]/V_conector-HR_to_conctor[reference_strip]/V_conector;
						} else {
							time_in_strip = calibrated_time - trig_time - HR_to_conctor[m_strip(frame[i])]/V_conector;
							time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[m_strip(frame[i])];
							time_to_conector_new = HR_to_conctor[m_strip(frame[i])]/V_conector-HR_to_conctor[reference_strip]/V_conector;
						}
						double HR_trig = time_in_strip_new + time_to_conector_new;
						if(jj==0) hHRT_->Fill(HR_trig); // Fill the signal - trig plot
						if(jj==0) hHRT_temp->Fill(HR_trig); // same but zoom // checar
                                                if(jj==0) hHRT_V.at(strip_numb)->Fill(HR_trig);
                                                if(jj==0) hHRT_zoom_V.at(strip_numb)->Fill(HR_trig);
						if(isFEBv2r2 && strip_numb==23)hHRT_strip23->Fill(HR_trig);
						if(jj==0 && retrig_study){
							if(first_trig_HR[strip_numb]==10000.) {first_trig_HR[strip_numb] = HR_trig;} else {
								hHRT_retrig_V.at(strip_numb)->Fill(HR_trig-first_trig_HR[strip_numb]);
								retrig_rate_temp++;
							}
				
						} 
						// add signal to the strip struct
						if(isFEBv2r2) {
							allStrips[strip_numb].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
						} else {
							allStrips[m_strip(frame[i])].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
						}
					}
				}
			} // end of first frame loop
			if(retrig_rate_temp>0 && retrig_study) retrig_rate++;

			int ntriggerHR_signal_strip = 0;
			int ntriggerLR_signal_strip = 0;

			//for clustering
			vector<int> strip_HR;
			vector<int> strip_LR;
			vector<double> time_HR;
			vector<double> time_LR;
			strip_HR.clear();
			strip_LR.clear();
			time_HR.clear();
			time_LR.clear();
			//loop over strips
			for (uint32_t i=0;i<48;i++) {
				int ntriggerHR_signal_per_strip = 0;
				for (uint32_t fhr=0; fhr<allStrips[i].HRframe.size(); fhr++){ // frame HR loop
					double time_in_strip = allStrips[i].HRtime[fhr] - trig_time/*trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]*/ - HR_to_conctor[i]/V_conector;
					double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
					double time_to_conector_new = HR_to_conctor[i]/V_conector-HR_to_conctor[reference_strip]/V_conector;
					double HR_trig = time_in_strip_new + time_to_conector_new;
					hHRT->Fill(i,HR_trig);
					if ( !((HR_trig) > muW1_HR and (HR_trig) < muW2_HR ))
					{
						allStrips[i].HRframeOK[fhr] = false;
						allStrips[i].HRframeBkg[fhr] = true;
						hHRn->Fill(i,1/((time_max-time_min-(muW2_HR-muW1_HR))*1e-9*0.9*Strip_length[i]*0.1)); //For Noise
					}
					if (allStrips[i].HRframeOK[fhr]){
						//strip_HR.push_back((int)i);
						//time_HR.push_back(time_in_strip_new); 
						ntriggerHR_signal++;
						ntriggerHR_signal_per_strip++;
					}
				}
				if(ntriggerHR_signal_per_strip>0) ntriggerHR_signal_strip++;
				nFiredHR_per_strip->Fill(i,ntriggerHR_signal_per_strip);

				int ntriggerLR_signal_per_strip = 0;
				for (uint32_t flr=0; flr<allStrips[i].LRframe.size(); flr++) { // frame LR loop
					double time_in_strip = allStrips[i].LRtime[flr] - trig_time/*trig[m_fpga(frame[allStrips[i].LRframe[flr]])]*/ - LR_to_conctor[i]/V_conector;
					double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
					double time_to_conector_new = LR_to_conctor[i]/V_conector-LR_to_conctor[reference_strip]/V_conector;
					double LR_trig = time_in_strip_new + time_to_conector_new;
					hLRT->Fill(i,LR_trig);
					if ( !((LR_trig) > muW1_LR and (LR_trig) < muW2_LR ))
					{
						allStrips[i].LRframeOK[flr] = false;
						allStrips[i].LRframeBkg[flr] = true;
						hLRn->Fill(i,1/((time_max-time_min-(muW2_LR-muW1_LR))*1e-9*0.9*Strip_length[i]*0.1)); //For Noise
					}
					if (allStrips[i].LRframeOK[flr]){
						//strip_LR.push_back((int)i);
						//time_LR.push_back(time_in_strip_new);
						ntriggerLR_signal++;
						ntriggerLR_signal_per_strip++;
					}

				}
				if(ntriggerLR_signal_per_strip>0) ntriggerLR_signal_strip++;
				nFiredLR_per_strip->Fill(i,ntriggerLR_signal_per_strip);

				allStrips[i].buildDeltas();
				if (allStrips[i].deltas.size() > 0) has_paired_strips = 1; //In the event there is at least 1 paired strip
				for (uint32_t d=0; d<allStrips[i].deltas.size(); d++){ //loop over deltas
					n_paired_srip++;
					hLHRn->Fill(i,1/(((muW2_HR-muW1_HR))*1e-9*0.9*Strip_length[i]*0.1)); 
					double time_corr_HR = time_corr_fine[m_fpga(frame[allStrips[i].deltas[d].second])][m_channel(frame[allStrips[i].deltas[d].second])];
					double time_corr_LR = time_corr_fine[m_fpga(frame[allStrips[i].deltas[d].first])][m_channel(frame[allStrips[i].deltas[d].first])];
					double time_in_strip = m_time(m_traw(frame[allStrips[i].deltas[d].second])) - time_corr_HR  - trig_time - HR_to_conctor[i]/V_conector;
					double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
					double time_to_conector_new = HR_to_conctor[i]/V_conector-HR_to_conctor[reference_strip]/V_conector;
					double HR_trig = time_in_strip_new + time_to_conector_new;
					double time_in_strip2 = m_time(m_traw(frame[allStrips[i].deltas[d].first])) - time_corr_LR - trig_time - LR_to_conctor[i]/V_conector;
					double time_in_strip2_new = time_in_strip2*Strip_length[reference_strip]/Strip_length[i];
					double time_to_conector_new2 = LR_to_conctor[i]/V_conector-LR_to_conctor[reference_strip]/V_conector;
					double LR_trig = time_in_strip2_new + time_to_conector_new2;

					strip_HR.push_back((int)i);
					time_HR.push_back(time_in_strip); //time_in_strip_new
					strip_LR.push_back((int)i);
					time_LR.push_back(time_in_strip2); //time_in_strip2_new

					deltaT=time_in_strip2_new-time_in_strip_new+Strip_length[reference_strip]/(V_conector);
					deltaT=deltaT/2.;//-i*0.3+7;//0.15 2
					if(deltaT>0) n_paired_srip_med++;

					vector<double> X_Y = convert_DeltaTAndStrip_To_XY(sn_,i,(time_in_strip-time_in_strip2-dT_align_factor[i]));//0.3 7.
                                        vector<double> X_Y_forAlign = convert_DeltaTAndStrip_To_XY(sn_,i,(time_in_strip-time_in_strip2));
					//std::cout << "test" << std::endl;
					if(align_XY && X_Y_forAlign.at(1)>-200 && X_Y_forAlign.at(1)<1500 && jj==1){
						Y_align[i] = Y_align[i]+X_Y_forAlign.at(1);
                                                Y_align_n[i] = Y_align_n[i]+1;
						dT_align[i] = dT_align[i]+(time_in_strip-time_in_strip2);
						dT_align_n[i] = dT_align_n[i]+1;
					}

                                        //vector<double> X_Y = convert_DeltaTAndStrip_To_XY(sn_,i,(m_time(m_traw(frame[allStrips[i].deltas[d].second]))-m_time(m_traw(frame[allStrips[i].deltas[d].first]))-time_corr_HR+time_corr_LR+10.));//0.3 7.

					h2_XY->Fill(X_Y.at(0),(X_Y.at(1)/*+Y_align_factor[i]*/));
					//std::cout << "X: " << X_Y.at(0) << "; Y: " << X_Y.at(1) << std::endl;

					sumT=time_in_strip2_new+time_in_strip_new;//LR_trig+HR_trig;
					//                                        hdeltaT->Fill(i,(m_time(m_traw(frame[allStrips[i].deltas[d].first]))-time_corr_LR+LR_to_conctor[i]/V_conector-(m_time(m_traw(frame[allStrips[i].deltas[d].second]))-time_corr_HR+HR_to_conctor[i]/V_conector)));
					hdeltaT->Fill(i,deltaT);
					hsumT->Fill(i,sumT);
					hdeltaT_1D->Fill(deltaT);
					if(jj==0) hdeltaT_1D_V.at(i)->Fill(time_in_strip-time_in_strip2/*deltaT*/);
					hsumT_1D->Fill(sumT);
				}
			} // strip loop ends

			bool cluster_good = false;

			//for cluster study
			if(plot_clust_par){
				double limit = 0.1;
				for(int ij=0;ij<100;ij++){
					vector<ClusterOneSide> cluster_HR = cluster_one_side(strip_HR, time_HR, limit);
					vector<ClusterOneSide> cluster_LR = cluster_one_side(strip_LR, time_LR, limit);

					std::vector<cluster> CLUSTER;
					CLUSTER.clear();
					if(cluster_HR.size()>0 && cluster_LR.size()>0) CLUSTER = final_clustering(cluster_LR,cluster_HR,Strip_length[reference_strip]/(V_conector),sn_);
					if(ClusterE==0) ClusterL.push_back(limit);
					if(ClusterE==0) ClusterN.push_back(CLUSTER.size()); else ClusterN.at(ij) = ClusterN.at(ij)+CLUSTER.size();
					if(ClusterE==0){
						double clustersizetemp=0.;
						for(int iii=0;iii<CLUSTER.size();iii++){
							clustersizetemp=clustersizetemp+CLUSTER.at(iii).size();
						}
						if(CLUSTER.size()>0) ClusterS.push_back(clustersizetemp/CLUSTER.size());
						else ClusterS.push_back(0.);
					} else {
						double clustersizetemp=0.;
						for(int iii=0;iii<CLUSTER.size();iii++){
							clustersizetemp=clustersizetemp+CLUSTER.at(iii).size();
						}
						if(CLUSTER.size()>0) ClusterS.at(ij) = ClusterS.at(ij)+clustersizetemp/CLUSTER.size();
					}
					limit = limit + 0.05;
				}
				ClusterE++;
			}

			// run cluster algos
			vector<ClusterOneSide> cluster_HR_new = cluster_one_side(strip_HR, time_HR, 3.);
			vector<ClusterOneSide> cluster_LR_new = cluster_one_side(strip_LR, time_LR, 3.);

			std::vector<cluster> CLUSTER;
			CLUSTER.clear();
			if(cluster_HR_new.size()>0 && cluster_LR_new.size()>0){
				CLUSTER = final_clustering(cluster_LR_new,cluster_HR_new,Strip_length[reference_strip]/(V_conector),sn_);
				// for cluster study
				if(plot_clust_par){
					double strip_diff = 999.;
					int LR_ind_smallest = -1;
					int HR_ind_smallest = -1;
					for(int ij=0;ij<cluster_HR_new.size();ij++){
						for(int ji=0;ji<cluster_LR_new.size();ji++){
							if(abs(cluster_HR_new.at(ij).strip() - cluster_LR_new.at(ji).strip())<strip_diff) {
								strip_diff = abs(cluster_HR_new.at(ij).strip() - cluster_LR_new.at(ji).strip());
								LR_ind_smallest = ji;
								HR_ind_smallest = ij;
							}
						}
					}
					strip_diff_h->Fill(cluster_HR_new.at(HR_ind_smallest).strip() - cluster_LR_new.at(LR_ind_smallest).strip());
				}
			}
			hNClusters->Fill(CLUSTER.size());
			sum_cluster_size = sum_cluster_size + CLUSTER.size();
			for(int ij=0;ij<CLUSTER.size();ij++){
				hsumTCluster_1D->Fill(CLUSTER.at(ij).sum_time());
				/*if(CLUSTER.at(ij).time()>0.)*/ cluster_good=true;
				if(CLUSTER.at(ij).size()>0.) hClusterSize->Fill(CLUSTER.at(ij).size());
				hdeltaTCluster_1D->Fill(CLUSTER.at(ij).time());
				hdeltaClusterT->Fill(CLUSTER.at(ij).strip(),(CLUSTER.at(ij).time()));
//                                int strip_cls_int = (int)CLUSTER.at(ij).strip();
//                                if((CLUSTER.at(ij).strip()-strip_cls_int)>0.5) strip_cls_int = strip_cls_int+1;
//				vector<double> X_Y_cls = convert_DeltaTAndStrip_To_XY(sn_,CLUSTER.at(ij).strip(),CLUSTER.at(ij).time_XY()/*+dT_align_factor[strip_cls_int]*/);
//				if(CLUSTER.at(ij).size()>=2) h2_XY_cls->Fill(X_Y_cls.at(0),X_Y_cls.at(1)/*+Y_align_factor[strip_cls_int]*/);
                                if(CLUSTER.at(ij).size()>=2) h2_XY_cls->Fill(CLUSTER.at(ij).X(),CLUSTER.at(ij).Y());
			}

			if(cluster_good) N_cluster_good++;
			if((ntriggerHR_signal>0) || (ntriggerLR_signal>0)) nANDstrip_noCut_fired++;
			if(ntriggerLR_signal>0) nANDstripLR_noCut_fired++;
			if(ntriggerHR_signal>0) nANDstripHR_noCut_fired++;

			ntrig_allevent++; // count denom of eff
			if (has_paired_strips) {
				nANDstrip_fired++;
				hHRT_2->Fill(ntriggerHR_signal_strip);
				hLRT_2->Fill(ntriggerLR_signal_strip);
			}
			hnPairs->Fill(n_paired_srip);
		} //Event loop ends
                std::cout << "===================================================================================" << std::endl;
		if(jj==0 && retrig_study) std::cout << "RETRIG RATE: " << (double)retrig_rate/(double)ntrig_allevent << std::endl;
                if(jj==0 && retrig_study) std::cout << "ntrig_allevent: " << ntrig_allevent << std::endl;
                std::cout << "===================================================================================" << std::endl;

		hHRT_temp->GetXaxis()->SetRangeUser(hHRT_temp->GetXaxis()->GetBinCenter(hHRT_temp->GetMaximumBin())-50,hHRT_temp->GetXaxis()->GetBinCenter(hHRT_temp->GetMaximumBin())+50);
		hLRT_temp->GetXaxis()->SetRangeUser(hLRT_temp->GetXaxis()->GetBinCenter(hLRT_temp->GetMaximumBin())-50,hLRT_temp->GetXaxis()->GetBinCenter(hLRT_temp->GetMaximumBin())+50);

		hHRT_strip23->GetXaxis()->SetRangeUser(hHRT_strip23->GetXaxis()->GetBinCenter(hHRT_strip23->GetMaximumBin())-50,hHRT_strip23->GetXaxis()->GetBinCenter(hHRT_strip23->GetMaximumBin())+50);
		hLRT_strip23->GetXaxis()->SetRangeUser(hLRT_strip23->GetXaxis()->GetBinCenter(hLRT_strip23->GetMaximumBin())-50,hLRT_strip23->GetXaxis()->GetBinCenter(hLRT_strip23->GetMaximumBin())+50);

		hsumT->GetYaxis()->SetRangeUser(hsumT_1D->GetXaxis()->GetBinCenter(hsumT_1D->GetMaximumBin())-50,hsumT_1D->GetXaxis()->GetBinCenter(hsumT_1D->GetMaximumBin())+50);
		hsumT_1D->GetXaxis()->SetRangeUser(hsumT_1D->GetXaxis()->GetBinCenter(hsumT_1D->GetMaximumBin())-50,hsumT_1D->GetXaxis()->GetBinCenter(hsumT_1D->GetMaximumBin())+50);
		hsumTCluster_1D->GetXaxis()->SetRangeUser(hsumTCluster_1D->GetXaxis()->GetBinCenter(hsumTCluster_1D->GetMaximumBin())-50,hsumTCluster_1D->GetXaxis()->GetBinCenter(hsumTCluster_1D->GetMaximumBin())+50);

		if(jj==0 && ClusterE!=0){
			for(int ij=0;ij<100;ij++){
				ClusterS.at(ij) = ClusterS.at(ij)/ClusterE;
				ClusterN.at(ij) = ClusterN.at(ij)/ClusterE;
				g_cluster_size->SetPoint(ij,ClusterL.at(ij),ClusterS.at(ij));
				g_cluster_number->SetPoint(ij,ClusterL.at(ij),ClusterN.at(ij));
			}
		}
		if(jj==0){
			or_eff_on = (nANDstrip_noCut_fired)/float(ntrig_allevent);
			and_eff_on = (nANDstrip_fired)/float(ntrig_allevent);
			HR_eff_on = (nANDstripHR_noCut_fired)/float(ntrig_allevent);
			LR_eff_on = (nANDstripLR_noCut_fired)/float(ntrig_allevent);
			cluster_eff_on = (N_cluster_good)/float(ntrig_allevent);
		}

		if(jj==1){
			or_eff_out = (nANDstrip_noCut_fired)/float(ntrig_allevent);
			and_eff_out = (nANDstrip_fired)/float(ntrig_allevent);
			HR_eff_out = (nANDstripHR_noCut_fired)/float(ntrig_allevent);
			LR_eff_out = (nANDstripLR_noCut_fired)/float(ntrig_allevent);
			cluster_eff_out = (N_cluster_good)/float(ntrig_allevent);
			cluster_rate_out = sum_cluster_size/(ntrig_allevent*muon_window_size*pow(10,-9)*active_area);
		}

		if(no_trig_case==true) and_eff_out = 0.;

		myfile << (nANDstrip_noCut_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstrip_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstripHR_noCut_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstripLR_noCut_fired)/float(ntrig_allevent) <<"\n";
		myfile << (N_cluster_good)/float(ntrig_allevent) <<"\n";
		myfile << sum_cluster_size/(ntrig_allevent*muon_window_size*pow(10,-9)*6000) <<"\n";
		if(jj==0) ntrig_allevent_muon_window = ntrig_allevent;
		if(jj==1) ntrig_allevent_gamma_window = ntrig_allevent;

		myfile<<""<<"\n";

		TString sss("");
		if(jj==0) sss.Form("muons_");
		if(jj==1) sss.Form("gammas_");
		if(jj==0 || (jj==1&&print_gamma_histos)){
			// Draw histos and save

			hLHRn_dif = (TH1F*)hHRn->Clone("hLHRn_dif");
			hLHRn_dif->Divide(hLRn);
			hLHRn_dif->GetYaxis()->SetTitle("HR/LR");			
			hLHRn_dif->SetTitle("Strips");
			hLHRn_dif->SetTitle("Rate ratio - HR/LR - with trig");
			plot_and_save_1D(hLHRn_dif, sss, "rate_ratio_with_trg", hv_, sn_);

			hLHR_dif = (TH1F*)hHR->Clone("hLHR_dif");
			hLHR_dif->Divide(hLR);
			hLHR_dif->GetYaxis()->SetTitle("HR/LR");
			hLHR_dif->SetTitle("Strips");
			hLHR_dif->SetTitle("Rate ratio - HR/LR - no trig");
			plot_and_save_1D(hLHR_dif, sss, "rate_ratio_no_trg", hv_, sn_);

			hLR_ns->GetYaxis()->SetTitle("events");
			hLR_ns->GetXaxis()->SetTitle("# all signals");
			hHR_ns->GetYaxis()->SetTitle("events");
			hHR_ns->GetXaxis()->SetTitle("# all signals");
			plot_and_save_2x_1D(hLR_ns, hHR_ns, sss, "all_signals", hv_,sn_, 1, false);

			hLR->GetYaxis()->SetTitle("Noise Hz/cm2");
			hLR->GetXaxis()->SetTitle("Strips");
			hHR->GetYaxis()->SetTitle("Noise Hz/cm2");
			hHR->GetXaxis()->SetTitle("Strips");
			plot_and_save_2x_1D(hLR, hHR, sss, "avarage_noise_side_no_trg", hv_,sn_, 0, true);

			if(jj==0) hLRn->Scale(1./ntrig_allevent_muon_window);
			if(jj==0) hHRn->Scale(1./ntrig_allevent_muon_window);
			if(jj==1) hLRn->Scale(1./ntrig_allevent_gamma_window);
			if(jj==1) hHRn->Scale(1./ntrig_allevent_gamma_window);
			if(jj==0) hLHRn->Scale(1./ntrig_allevent_muon_window);
			if(jj==1) hLHRn->Scale(1./ntrig_allevent_gamma_window);

			hLRn->GetYaxis()->SetTitle("Noise Hz/cm2");
			hLRn->GetXaxis()->SetTitle("Strips");
			hHRn->GetYaxis()->SetTitle("Noise Hz/cm2");
			hHRn->GetXaxis()->SetTitle("Strips");
			plot_and_save_2x_1D(hLRn, hHRn, sss, "avarage_noise_side_with_trg", hv_,sn_, 0, true);

			hLHRn->GetYaxis()->SetTitle("Noise Hz/cm2");
			hLHRn->GetXaxis()->SetTitle("Strips");
			plot_and_save_1D(hLHRn, sss, "avarage_noise_side_with_trg_AND", hv_,sn_,1111);


			hdeltaT->SetTitle("");
			hdeltaT->GetYaxis()->SetTitle("delta time [ns]");
			hdeltaT->GetXaxis()->SetTitle("Strips");
			plot_and_save_2D(hdeltaT, sss, "deltaT_srip", hv_, sn_, "P");


			h2_XY->SetTitle("");
			h2_XY->GetYaxis()->SetTitle("y [mm]");
			h2_XY->GetXaxis()->SetTitle("x [mm]");
			h2_XY->SetMarkerSize(4);
			plot_and_save_2D(h2_XY, sss, "X_Y", hv_, sn_, "P");

                        h2_XY_cls->SetTitle("");
                        h2_XY_cls->GetYaxis()->SetTitle("y [mm]");
                        h2_XY_cls->GetXaxis()->SetTitle("x [mm]");
                        h2_XY_cls->SetMarkerSize(4);
                        plot_and_save_2D(h2_XY_cls, sss, "X_Y_cls", hv_, sn_, "P");

			hdeltaClusterT->SetTitle("");
			hdeltaClusterT->GetYaxis()->SetTitle("Cluster delta time [ns]");
			hdeltaClusterT->GetXaxis()->SetTitle("Cluster Strips");
			plot_and_save_2D(hdeltaClusterT, sss, "Cluster_deltaT_srip", hv_, sn_, "P");

			hsumT->SetTitle("");
			hsumT->GetYaxis()->SetTitle("sum time [ns]");
			hsumT->GetXaxis()->SetTitle("Strips");
			//plot_and_save_2D(hsumT, sss, "sumT_srip", hv_, sn_, "P");

			hsumTCluster_1D->GetXaxis()->SetTitle("Cluster sum time [ns]");
			hsumTCluster_1D->SetTitle("");
			//plot_and_save_1D(hsumTCluster_1D, sss, "Cluster_sumT_srip", hv_, sn_,1111);

			hnPairs->GetYaxis()->SetTitle("Events");
			hnPairs->GetXaxis()->SetTitle("n paired strips");
			hnPairs->SetTitle("");
			plot_and_save_1D(hnPairs,sss, "n_paired_srip", hv_, sn_,1111);


			hClusterSize->SetTitle("");
			hClusterSize->GetYaxis()->SetTitle("Events");
			hClusterSize->GetXaxis()->SetTitle("Cluster size");
			plot_and_save_1D(hClusterSize,sss, "Cluster_size", hv_, sn_,1111);

			hNClusters->GetYaxis()->SetTitle("Events");
			hNClusters->GetXaxis()->SetTitle("Number of Clusters");
			hNClusters->SetTitle("");
			plot_and_save_1D(hNClusters,sss, "N_Clusters", hv_, sn_,1111);

			hdeltaT_1D->GetXaxis()->SetTitle("delta time [ns]");
			hdeltaT_1D->GetYaxis()->SetTitle("Events");
			hdeltaT_1D->SetTitle("");
			plot_and_save_1D(hdeltaT_1D, sss, "deltaT_1D_srip", hv_, sn_);

			hdeltaTCluster_1D->SetTitle("");
			hdeltaTCluster_1D->GetXaxis()->SetTitle("Cluster delta time [ns]");
			hdeltaTCluster_1D->GetYaxis()->SetTitle("Events");
			plot_and_save_1D(hdeltaTCluster_1D, sss, "Cluster_deltaT_1D_srip", hv_, sn_);

			hsumT_1D->SetTitle("");
			hsumT_1D->GetXaxis()->SetTitle("sum time [ns]");
			hsumT_1D->GetYaxis()->SetTitle("Events");
			//plot_and_save_1D(hsumT_1D, sss, "sumT_1D_srip", hv_, sn_);

			plot_and_save_2x_2D(hLRT, hHRT, sss, "DeltaTrig_side", hv_, sn_, "COLZ");
			//plot_and_save_2D(nFiredHR_per_strip, sss, "nFiredHR_per_strip", hv_, sn_, "COLZ");
			//plot_and_save_2D(nFiredLR_per_strip, sss, "nFiredLR_per_strip", hv_, sn_, "COLZ");

			plot_and_save_2x_1D(hLRT_temp, hHRT_temp, sss, "DeltaTrig_side_forTimeWindow", hv_,sn_, 0);
			plot_and_save_2x_1D(hLRT_, hHRT_, sss, "DeltaTrig_side_1d", hv_,sn_, 0);
			plot_and_save_2x_1D(hLRT_2, hHRT_2, sss, "nStrip_side_1d", hv_,sn_, 1111);

			plot_and_save_2x_1D(hLRT_strip23, hHRT_strip23, sss, "DeltaTrig_side_forTimeWindow_strip23", hv_,sn_, 0);

			if(jj==0){
				for (int i=0;i<48;i++){
					plot_and_save_2x_1D(hLRT_V.at(i), hHRT_V.at(i), sss, Form("DeltaTrig_side_1d_strip%d",i), hv_,sn_, 0, false, "plts_per_strip/");
                	                plot_and_save_2x_1D(hLRT_zoom_V.at(i), hHRT_zoom_V.at(i), sss, Form("DeltaTrig_zoom_side_1d_strip%d",i), hv_,sn_, 0, false, "plts_per_strip/");
                        		hdeltaT_1D_V.at(i)->GetXaxis()->SetTitle(Form("delta time [ns] strip%d",i));
                	        	hdeltaT_1D_V.at(i)->GetYaxis()->SetTitle("Events");
        	       			hdeltaT_1D_V.at(i)->SetTitle("");
	                        	plot_and_save_1D(hdeltaT_1D_V.at(i), sss, Form("deltaT_1D_srip%d",i), hv_, sn_, 0, "plts_per_strip/");
                        		//for retrig study
					if(retrig_study) plot_and_save_2x_1D(hLRT_retrig_V.at(i), hHRT_retrig_V.at(i), sss, Form("DeltaTrig_retrig_1d_strip%d",i), hv_,sn_, 0, false, "retriggering_study/");
//                                        plot_and_save_1D(hdeltaT_1D_V.at(i), sss, Form("deltaT_1D_srip%d",i), hv_, sn_, 0, "retriggering_study/");
					if(i!=0 && retrig_study) hLRT_retrig_V.at(0)->Add(hLRT_retrig_V.at(i));
                                        if(i!=0 && retrig_study) hHRT_retrig_V.at(0)->Add(hHRT_retrig_V.at(i));
				}
				if(retrig_study) plot_and_save_1D(hLRT_retrig_V.at(0), sss, Form("deltaT_1D_LR"), hv_, sn_, 0, "retriggering_study/");
                                if(retrig_study) plot_and_save_1D(hHRT_retrig_V.at(0), sss, Form("deltaT_1D_HR"), hv_, sn_, 0, "retriggering_study/");
				if(retrig_study) hLRT_retrig_V.at(0)->Add(hHRT_retrig_V.at(0));
				if(retrig_study) plot_and_save_1D(hLRT_retrig_V.at(0), sss, Form("deltaT_1D_ALL"), hv_, sn_, 0, "retriggering_study/");
			}

			//plot_and_save_1D(hbc0, sss, "hbc0", hv_, sn_);

			if(jj==0) hHRn_gamma_AND = hLHRn->Integral()/48.;
			if(jj==0) hHRn_gamma_HR = hHRn->Integral()/48.;
			if(jj==0) hHRn_gamma_LR = hLRn->Integral()/48.;
			if(jj==0) gamma_cluster_size = hClusterSize->GetMean(); 

			hLR->Reset(); hHR->Reset(); hLRn->Reset(); hHRn->Reset(); nFiredHR_per_strip->Reset(); nFiredLR_per_strip->Reset(); hdeltaT->Reset(); hsumT->Reset(); hdeltaT_1D->Reset(); hsumT_1D->Reset(); hLRT_->Reset(); hHRT_->Reset(); hLRT_temp->Reset(); hHRT_temp->Reset(); hLRT_2->Reset(); hHRT_2->Reset(); hnPairs->Reset(); hdeltaTCluster_1D->Reset(); hsumTCluster_1D->Reset(); hClusterSize->Reset(); hNClusters->Reset(); hdeltaClusterT->Reset(); hLRT_strip23->Reset(); hHRT_strip23->Reset(); hLHRn_dif->Reset(); hLHR_dif->Reset(); hLR_ns->Reset(); hHR_ns->Reset(); h2_XY->Reset(); h2_XY_cls->Reset(); hLHRn->Reset(); 
		}

	}
	myfile<<""<<"\n";
	myfile<<""<<"\n";
	myfile<<""<<"\n";
	myfile << "final OR eff:              " << (or_eff_on-or_eff_out)/(1-or_eff_out) <<"\n";
	myfile << "final AND eff:             " << (and_eff_on-and_eff_out)/(1-and_eff_out)  <<"\n";
	myfile << "final LR eff:              " << (HR_eff_on-HR_eff_out)/(1-HR_eff_out)  <<"\n";
	myfile << "final HR eff:              " << (LR_eff_on-LR_eff_out)/(1-LR_eff_out)  <<"\n";
	myfile << "final CLUSTER eff:         " << (cluster_eff_on-cluster_eff_out)/(1-cluster_eff_out)  <<"\n";
	myfile << "final gamma cluster rate:  " << cluster_rate_out/((and_eff_on-and_eff_out)/(1-and_eff_out))  <<"\n";
	myfile << "AND hit rate average divided by gamma cluster size: " << (hHRn_gamma_AND/gamma_cluster_size) << "\n";
	myfile << "HR hit rate average divided by gamma cluster size: " << (hHRn_gamma_HR/gamma_cluster_size) << "\n";
	myfile << "LR hit rate average divided by gamma cluster size: " << (hHRn_gamma_LR/gamma_cluster_size) << "\n";

	myfile.close();

	

	//for Y aligment
	if(align_XY){
                for(int i = 0; i<48; i++){
			std::cout << "-->> " << i << std::endl;
			std::cout << hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindFirstBinAbove(20))  << ";  " <<hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindLastBinAbove(20)) << std::endl;
                        std::cout << 0.5*(hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindLastBinAbove(20))-hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindFirstBinAbove(20)))*V_conector << "; " << Strip_length[i] << "  -->> " <<  (0.5*(hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindLastBinAbove(20))-hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindFirstBinAbove(20)))*V_conector-Strip_length[i])  << std::endl;
//        	        std::cout << std::endl;
//	                std::cout << "=============================" << std::endl;
//			std::cout << (abs(hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindLastBinAbove(20)))-abs(hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindFirstBinAbove(20))))/2  << ",";
//                	std::cout << std::endl;
//                	std::cout << "=============================" << std::endl;
		}
                std::cout << "=============================" << std::endl;
                for(int i = 0; i<48; i++){
                        std::cout << (abs(hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindLastBinAbove(20)))-abs(hdeltaT_1D_V.at(i)->GetXaxis()->GetBinCenter(hdeltaT_1D_V.at(i)->FindFirstBinAbove(20))))/2  << ",";
                }
                std::cout << std::endl;
                std::cout << "=============================" << std::endl;

/*		int best_strip = 0;
		double best_strip_diff = abs(Y_align[0]/Y_align_n[0]-Strip_length[0]);
		for(int i = 0; i<48; i++){
			std::cout << "Y_align[i]/Y_align_n[i]: " << Y_align[i]/Y_align_n[i] << ";  Y_align_n[i]: " << Y_align_n[i] << std::endl; 
			if(abs(Y_align[i]/Y_align_n[i]-Strip_length[i])<best_strip_diff){
				best_strip = i;
				best_strip_diff = abs(Y_align[i]/Y_align_n[i]-Strip_length[i]);
			}
		}
		std::cout << "Best strip: " << best_strip << "; Y-strip aligment factors:" << std::endl;
        	for(int i = 0; i<48; i++){
			std::cout << (Y_align[best_strip]/Y_align_n[best_strip]-Y_align[i]/Y_align_n[i])  << ",";
		}
		std::cout << std::endl;
		std::cout << "=============================" << std::endl;
		best_strip_diff=99999;
                for(int i = 0; i<48; i++){
                        if(abs(dT_align[i]/dT_align_n[i])<best_strip_diff){
                                best_strip = i;
				std::cout << "best strip in time: " << i << "; time:" << abs(dT_align[i]/dT_align_n[i]) << std::endl;
                                best_strip_diff = abs(dT_align[i]/dT_align_n[i]);
                        }
                }
                for(int i = 0; i<48; i++){
                        std::cout << (dT_align[best_strip]/dT_align_n[best_strip]-dT_align[i]/dT_align_n[i])  << ",";
                }
                std::cout << std::endl;
*/	}
	//for clustering study
	if(plot_clust_par){
		gStyle->SetOptFit(1111);
		TCanvas *cstrip_diff_h = new TCanvas("cstrip_diff_h","strip_diff_h",200,10,780,780);
		cstrip_diff_h->cd();
		cstrip_diff_h->SetGrid();
		strip_diff_h->Fit("gaus","","",-5,5);
		strip_diff_h->Draw();
		s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_strip_diff.root",sn_,hv_,hv_,sn_);
		cstrip_diff_h->SaveAs(s);
		s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_strip_diff.pdf",sn_,hv_,hv_,sn_);
		cstrip_diff_h->SaveAs(s);
		s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_strip_diff.png",sn_,hv_,hv_,sn_);
		cstrip_diff_h->SaveAs(s);
		gStyle->SetOptFit(0);


		g_cluster_size->SetMarkerStyle(7);
		g_cluster_size->SetMarkerColor(1);
		g_cluster_size->SetLineColor(1);
		g_cluster_size->SetFillColor(3);
		g_cluster_number->SetMarkerStyle(7);
		g_cluster_number->SetMarkerColor(2);
		g_cluster_number->SetLineColor(2);
		g_cluster_number->SetFillColor(3);
		TCanvas *cg = new TCanvas("cg", "comp", 200, 10, 700, 500);
		TPad *p1g = new TPad("p1g", "", 0, 0, 1, 1);
		TPad *p2g = new TPad("p2g", "", 0, 0, 1, 1);
		p2g->SetFillStyle(4000);
		p1g->Draw();
		p1g->cd();
		g_cluster_size->Draw("ALP");
		gPad->Update();
		Double_t xmin = p1g->GetUxmin();
		Double_t xmax = p1g->GetUxmax();
		Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
		Double_t ymin = g_cluster_number->GetHistogram()->GetMinimum();
		Double_t ymax = g_cluster_number->GetHistogram()->GetMaximum();
		Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
		p2g->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
		p2g->Draw();
		p2g->cd();
		g_cluster_number->Draw("LP");
		gPad->Update();
		TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
		axis->SetLineColor(kRed);
		axis->SetLabelColor(kRed);
		axis->Draw();
		gPad->Update();
		cg->cd();
		cg->Update();
		s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_cluster_number_size.root",sn_,hv_,hv_,sn_);
		cg->SaveAs(s);
		s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_cluster_number_size.pdf",sn_,hv_,hv_,sn_);
		cg->SaveAs(s);
		s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_cluster_number_size.png",sn_,hv_,hv_,sn_);
		cg->SaveAs(s);
	}

	gApplication->Terminate();

}


