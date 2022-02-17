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


using namespace std;

void ana_FEBv2::Loop()
{

	bool isFEBv2r2 = true;

	if (fChain == 0) return;
	TString         s("");
	s.Form("will save plots for RUN: , _HV_%d_SN_%d_MaxTrig_ ",hv_,sn_/*,mt_*/);

	std::cout<<s<<std::endl;
	double time_min = 4000;
	double time_max = 87000;

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

	vector<double> Strip_length_vec;
	Strip_length_vec.clear();
	vector<double> LR_to_conctor_vec;
	LR_to_conctor_vec.clear();
	vector<double> HR_to_conctor_vec;
	HR_to_conctor_vec.clear();

	HR_to_conctor_vec = Aligment_factor(sn_,1,1);
        LR_to_conctor_vec = Aligment_factor(sn_,2,1);
        Strip_length_vec = Aligment_factor(sn_,3,1);

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
        for(int iref = 1; iref < 48; iref++){
                if(Strip_length[iref] < Strip_length[reference_strip]) reference_strip = iref;
        }

	//reference_strip = 1;//

	double V_conector = 299.792*0.6;// mm/ns  - check it


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
	float sumT = -999;
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
	int nANDstripHR_noCut_fired=0;
	int nANDstripLR_noCut_fired=0;
        int N_cluster_good=0;
        int N_cluster_good_2=0;
	Long64_t trig_index = 0;
	float return_SR_int = 0;
	float return_BR_int = 0;
	float direct_SR_int = 0;
	float direct_BR_int = 0;
	float eff_medium = 0;
	int n_paired_srip = 0;
        int n_paired_srip_med = 0;


        TH1F* hLR_deltaTime_nextStrip = new TH1F("hLR_deltaTime_nextStrip", "LR neighbor strip delta time; ns", 100, -5., 5.);
        TH1F* hHR_deltaTime_nextStrip = new TH1F("hHR_deltaTime_nextStrip", "HR neighbor strip delta time; ns", 100, -5., 5.);


        TH1F* histo_diff_petiroc  = new TH1F("histo_diff_petiroc", "diff time for same petiroc", 60, -30, 30);

	TH1F* hLR = new TH1F("hLR", "LR", 50, 0.0, 50);
	TH1F* hHR = new TH1F("hHR", "HR", 50, 0.0, 50);
	TH1F* hLRn = new TH1F("hLRn", "hLRn", 50, 0.0, 50);
	TH1F* hHRn = new TH1F("hHRn", "hHRn", 50, 0.0, 50);
	TH2F* nPaired_per_strip = new TH2F("nPaired_per_strip", "# paired strip ; strip; #", 50,0,50,5, 0, 5);
	TH2F* nFiredHR_per_strip = new TH2F("nFiredHR_per_strip", "# fired strip ; strip; #", 50,0,50,10, 1, 11);
	TH2F* nFiredLR_per_strip = new TH2F("nFiredLR_per_strip", "# fired strip ; strip; #", 50,0,50,10, 1, 11);
	TH2F* hdeltaT = new TH2F("deltaT", "T (HR - LR)", 50,0,50,300, -30, 30);
	TH2F* hsumT = new TH2F("sumT", "T (HR + LR)", 50,0,50,100, -2750, -2650);
	TH1F* hdeltaT_1D = new TH1F("deltaT", "T (HR - LR)", 300, -30, 30);
	TH1F* hsumT_1D = new TH1F("sumT", "T (HR + LR)", 100, -1850, -1750);
	TH2F* hsumTdeltaT_2D = new TH2F("sumT_vs_deltaT", "sumT vs deltaT; T (HR + LR); T (HR - LR)", 1000, -2750, -2650, 300, -30, 30);
	TH1F* hLRT_ = new TH1F("hLRT_","T (LR - Trig)" ,10000,-9000,1000);
	TH1F* hHRT_ = new TH1F("hHRT_","T (HR - Trig)" ,10000,-9000,1000);
	TH1F* hLRT_temp = new TH1F("hLRT_muon_window","T (LR - Trig)" ,10000,-9000,1000);
	TH1F* hHRT_temp = new TH1F("hHRT_muon_window","T (HR - Trig)" ,10000,-9000,1000);
	TH1F* hLRT_2 = new TH1F("hLRT_2","n strip LR" ,20,0,20);
	TH1F* hHRT_2 = new TH1F("hHRT_2","n strip HR" ,20,0,20);
	TH1F* hLRT_3 = new TH1F("hLRT_3","n strip LR (LR - HR > 3 ns)" ,20,0,20);
	TH1F* hHRT_3 = new TH1F("hHRT_3","n strip HR (LR - HR > 3 ns)" ,20,0,20);
	TH1F* hnPairs = new TH1F("hnPairs","n paired strips" ,20,0,20);

	TH1F* hsumT_1D_med = new TH1F("sumT_med", "T (HR + LR)", 100, -2750, -2650);
	TH1F* hLRT__med = new TH1F("hLRT__med","T (LR - Trig)" ,10000,-9000,1000);
	TH1F* hHRT__med = new TH1F("hHRT__med","T (HR - Trig)" ,10000,-9000,1000);
	TH1F* hLRT_temp_med = new TH1F("hLRT_muon_window_med","T (LR - Trig)" ,10000,-9000,1000);
	TH1F* hHRT_temp_med = new TH1F("hHRT_muon_window_med","T (HR - Trig)" ,10000,-9000,1000);
	TH1F* hLRT_2_med = new TH1F("hLRT_2_med","n strip LR" ,20,0,20);
	TH1F* hHRT_2_med = new TH1F("hHRT_2_med","n strip HR" ,20,0,20);
	TH1F* hnPairs_med = new TH1F("hnPairs_med","n paired strips" ,20,0,20);

        TH1F* hdeltaTCluster_1D = new TH1F("hdeltaTCluster_1D", "Cluster T (HR - LR)", 300, -30, 30);
        TH1F* hsumTCluster_1D = new TH1F("hsumTCluster_1D", "Cluster T (HR + LR)", 100, -1850, -1750);
        TH1F* hClusterSize = new TH1F("hClusterSize","Cluster Size" ,10,0,10);
        TH1F* hNClusters = new TH1F("hNClusters","Number of Clusters" ,10,0,10);
        TH2F* hdeltaClusterT = new TH2F("hdeltaClusterT", "ClusterT (HR - LR)", 50,0,50,300, -30, 30);

        TH1F* hdeltaTCluster_2_1D = new TH1F("hdeltaTCluster_2_1D", "Cluster T (HR - LR)", 300, -30, 30);
        TH1F* hsumTCluster_2_1D = new TH1F("hsumTCluster_2_1D", "Cluster T (HR + LR)", 100, -2750, -2650);
        TH1F* hClusterSize_2 = new TH1F("hClusterSize_2","Cluster Size" ,10,0,10);
        TH1F* hNClusters_2 = new TH1F("hNClusters_2","Number of Clusters" ,10,0,10);
        TH2F* hdeltaClusterT_2 = new TH2F("hdeltaClusterT_2", "ClusterT (HR - LR)", 50,0,50,300, -30, 30);

        TH1F* strip_diff_h = new TH1F("strip_diff_h", " HR-LR strip cluster difference ", 19, -5, 5);

        TH1F* hHR_check_xtalk = new TH1F("hHR_check_xtalk", "HR signal in next strip", 50, 0.0, 50);
        TH1F* hLR_check_xtalk = new TH1F("hLR_check_xtalk", "LR signal in next strip", 50, 0.0, 50);

        TH1F* hHR_check_xtalk_temp = new TH1F("hHR_check_temp", "", 50, 0.0, 50);
        TH1F* hLR_check_xtalk_temp = new TH1F("hLR_check_temp", "", 50, 0.0, 50);

	
	TGraph *g_cluster_size = new TGraph();
	g_cluster_size->SetMarkerStyle(7);
	g_cluster_size->SetMarkerColor(1);
	g_cluster_size->SetLineColor(1);
	g_cluster_size->SetFillColor(3);
        TGraph *g_cluster_number = new TGraph();
	g_cluster_number->SetMarkerStyle(7);
	g_cluster_number->SetMarkerColor(2);
	g_cluster_number->SetLineColor(2);
	g_cluster_number->SetFillColor(3);

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
	TH2F *hHREvt[100]; //For histograms for each event
	TH2F *hLREvt[100];
	ofstream myfile;
	s.Form("outputs/_HV_%d_SN_%d_MaxTrig_.txt",hv_,sn_/*,mt_*/);
	myfile.open(s);


	auto newtree = fChain->CloneTree();


        vector<double> mu_wind_vec = MuonWindow(sn_, newtree, /*nframe, frame,*/ time_min, time_max, nentries/*, b_nframe, b_frame*/);
//	cout << "LR center: " << mu_wind_vec.at(0) << ";  HR center: " <<  mu_wind_vec.at(1) << endl;


        double muW1_HR = mu_wind_vec.at(1)-20.;
        double muW2_HR = mu_wind_vec.at(1)+20.;
        double muW1_LR = mu_wind_vec.at(0)-20.;
        double muW2_LR = mu_wind_vec.at(0)+20.;
        double muW1_bkg = muW1_HR-1000.;
        double muW2_bkg = muW2_HR-1000.;

	TH2F* hHRT = new TH2F("hHRT", "T (HR - Trig)", 50,0,50,abs(muW1_HR-muW2_HR), muW1_HR, muW2_HR);
	TH2F* hLRT = new TH2F("hLRT", "T (LR - Trig)", 50,0,50,abs(muW1_LR-muW2_LR), muW1_LR, muW2_LR);

	TH2F* hHRT_med = new TH2F("hHRT_med", "T (HR - Trig)", 50,0,50,abs(muW1_HR-muW2_HR), muW1_HR, muW2_HR);
	TH2F* hLRT_med = new TH2F("hLRT_med", "T (LR - Trig)", 50,0,50,abs(muW1_LR-muW2_LR), muW1_LR, muW2_LR);

	run_duration_in_sec = 10*(abs(muW1_HR)-abs(muW2_HR))*pow(10,-9);
	for (uint32_t i=0;i<100;i++) {
		hHREvt[i] = new TH2F(s, "T (return - Trig)", 50,0,50,abs(muW1_HR-muW2_HR),muW1_HR,muW2_HR);
		hLREvt[i] = new TH2F(s, "T (direct - Trig)", 50,0,50,abs(muW1_LR-muW2_LR),muW1_LR,muW2_LR);
	}


	double cluster_eff = 1;
	double cluster2_eff = 1;


	for (int jj=0;jj<2;jj++){
		double sum_cluster_size = 0;
                double sum_cluster2_size = 0;
		if(jj==1){
			muW1_HR = muW1_bkg;
			muW2_HR = muW2_bkg;
			muW1_LR = muW1_bkg;
			muW2_LR = muW2_bkg;
		}

		ntrig_allevent=0;
		nANDstrip_fired=0;
		nANDstrip_noCut_fired=0;
		nANDstripHR_noCut_fired=0;
		nANDstripLR_noCut_fired=0;
		nANDstrip_medium_fired=0;
		nANDstrip_tight_fired=0;
		N_cluster_good=0;
                N_cluster_good_2=0;


		vector<double> ClusterN;
                vector<double> ClusterS;
                vector<double> ClusterL;
		ClusterN.clear();
                ClusterS.clear();
                ClusterL.clear();
		int ClusterE = 0;

		int number_of_good_clusters = 0;
		int number_of_good_clusters_ratio1 = 0;
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry < 0) break;
			nb = fChain->GetEntry(jentry);   nbytes += nb;
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
                        n_paired_srip_med=0;
			medium_num=0;
			tight_num=0;
			double trig_time = -99999.;


                        for (uint32_t i=0;i<nframe;i++)
                        {

                        	if (m_channel(frame[i])==33){
                                        if (trig[m_fpga(frame[i])] != 0) bad_trig = 1; // prevents firing of one of the fpga's more than once
                                        trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
                                        if ((trig[m_fpga(frame[i])]) < time_max and trig[m_fpga(frame[i])] > time_min) {
                                                trig_exists=1; // trigger exits in this event
                                                ntrig_event++; // This value must reach 3 for all events.
                                                trig_time = m_time(m_traw(frame[i])); //trig[m_fpga(frame[i])];trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]
                                        }
				}


				std::vector<int> strip_and_side = m_strip_FEBv2r2(frame[i]);
        			int strip_side = -9999;
        			int strip_numb = -9999;
        			if(strip_and_side.size()>0){
			        	strip_side = strip_and_side.at(1);
			                strip_numb = strip_and_side.at(0);	
				} else continue;


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
                                                        if(jj==0) hLR->Fill(strip_numb,1/(last_bc0*(time_max-time_min-(muW2_LR-muW1_LR))*1e-9*1*Strip_length[strip_numb])); //For Noise
                                                } else {
                                                        if(jj==0) hLR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min-(muW2_LR-muW1_LR))*1e-9*1*Strip_length[m_strip(frame[i])])); //For Noise
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
			        			if(jj==0) hHR->Fill(strip_numb,1/(last_bc0*(time_max-time_min-(muW2_HR-muW1_HR))*1e-9*1*Strip_length[strip_numb])); //For Noise
			        		} else {
			        			if(jj==0) hHR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min-(muW2_HR-muW1_HR))*1e-9*1*Strip_length[m_strip(frame[i])])); //For Noise
			     	   		}
					} 
				}
			}



			if(trig_time==-99999.) continue;
			
			vector<int> petiroc;
			petiroc.clear();
			vector<double> time_petiroc;
			time_petiroc.clear();
                        for (uint32_t i=0;i<nframe;i++){
				std::vector<int> strip_and_side = m_strip_FEBv2r2(frame[i]);
	                        int strip_side = -9999;
        	                int strip_numb = -9999;
                	        if(strip_and_side.size()>0){
                        	        strip_side = strip_and_side.at(1);
                                	strip_numb = strip_and_side.at(0);
                        	}

				bool there_is_HR = true;
				if(m_channel(frame[i])!=32 && m_channel(frame[i])!=33){
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
                                                if ((HR_trig) > muW1_HR and (HR_trig) < muW2_HR ){
							for(int pr=0;pr<petiroc.size();pr++){
								if(petiroc.at(pr)==c_petiroc(m_channel(frame[i]))){ 
									there_is_HR=false;
									//histo_diff_petiroc->Fill((time_petiroc.at(pr)-m_time(m_traw(frame[i]))));
								}
							}
							if(there_is_HR){
								petiroc.push_back(c_petiroc(m_channel(frame[i])));
								//time_petiroc.push_back(m_time(m_traw(frame[i])));
							}
						}
					}					
				}
			}
			for(int ii = 0; ii<petiroc.size(); ii++){
				double time_temp = 9999999999999.;
	                        for (uint32_t i=0;i<nframe;i++){
					std::vector<int> strip_and_side = m_strip_FEBv2r2(frame[i]);
		                        int strip_side = -9999;
                		        int strip_numb = -9999;
                        		if(strip_and_side.size()>0){
                                		strip_side = strip_and_side.at(1);
                                		strip_numb = strip_and_side.at(0);
                        		}


                                        if(petiroc.at(ii)==c_petiroc(m_channel(frame[i]))){
	                                	if(m_channel(frame[i])!=32 && m_channel(frame[i])!=33){
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
                                                		if ((HR_trig) > muW1_HR and (HR_trig) < muW2_HR ){
									if(time_temp>m_time(m_traw(frame[i]))) time_temp = m_time(m_traw(frame[i]));
								}
                                			}
                        			}
					}
	                        }
				time_petiroc.push_back(time_temp);
			}

			for (uint32_t i=0;i<nframe;i++){
				std::vector<int> strip_and_side = m_strip_FEBv2r2(frame[i]);
	                        int strip_side = -9999;
        	                int strip_numb = -9999;
                	        if(strip_and_side.size()>0){
                        	        strip_side = strip_and_side.at(1);
                                	strip_numb = strip_and_side.at(0);
                        	}

				if (m_channel(frame[i])==32) {
					t_bc0[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
				}
				else if (m_channel(frame[i])==33){/*
					if (trig[m_fpga(frame[i])] != 0) bad_trig = 1; // prevents firing of one of the fpga's more than once
					trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
					if ((trig[m_fpga(frame[i])]) < time_max and trig[m_fpga(frame[i])] > time_min) {
						//std::cout<<"fpga: "<< m_fpga(frame[i]) << "Trig Exists !!!!! " << std::endl;
						trig_exists=1; // trigger exits in this event
						ntrig_event++; // This value must reach 3 for all events.
						trig_time = m_time(m_traw(frame[i])); //trig[m_fpga(frame[i])];trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]
					}
					if (debug) myfile<<"\t "<<"trig time "<< m_time(m_traw(frame[i])) <<"\n";
				*/}
				else {
					any_strip_fired = 1; // In this event there is at least one stripped fired. Without any requirement.
					if ((!isFEBv2r2 && c_side(m_channel(frame[i]))>0.5) || (isFEBv2r2 && strip_side==1)){
						//LR stips

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
						if(jj==0) hLRT_->Fill(LR_trig);
						if ((LR_trig) > muW1_LR and (LR_trig) < muW2_LR ){

							if(isFEBv2r2) {
								allStrips[strip_numb].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
                                                        	//if(jj==0) hLR->Fill(strip_numb,1/(last_bc0*(time_max-time_min)*1e-9*120)); //For Noise
							} else {
								allStrips[m_strip(frame[i])].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
								//if(jj==0) hLR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min)*1e-9*120)); //For Noise
							}
						} /*else {
                                                        if(isFEBv2r2) {
                                                                if(jj==0) hLR->Fill(strip_numb,1/(last_bc0*(time_max-time_min-(muW2_LR-muW1_LR))*1e-9*1*Strip_length[strip_numb])); //For Noise
                                                        } else {
                                                                if(jj==0) hLR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min-(muW2_LR-muW1_LR))*1e-9*1*Strip_length[m_strip(frame[i])])); //For Noise
                                                        }
						}*/
					}
					if ((!isFEBv2r2 && c_side(m_channel(frame[i]))<0.5) || (isFEBv2r2 && strip_side==0)){
						//HR strips
						bool pass = false;
                                                for(int pr=0;pr<petiroc.size();pr++){
                                                        if(petiroc.at(pr)==c_petiroc(m_channel(frame[i]))){ 
								histo_diff_petiroc->Fill((time_petiroc.at(pr)-m_time(m_traw(frame[i]))));
								if((time_petiroc.at(pr)-m_time(m_traw(frame[i])))>-1) pass = true; 
							}
						}
						//if(!pass) continue;

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
						if(jj==0) hHRT_->Fill(HR_trig);
                                                if ((HR_trig) > muW1_HR and (HR_trig) < muW2_HR ){
							if(isFEBv2r2) {
								allStrips[strip_numb].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
							        //if(jj==0) hHR->Fill(strip_numb,1/(last_bc0*(time_max-time_min)*1e-9*120)); //For Noise
							} else {
								allStrips[m_strip(frame[i])].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
							        //if(jj==0) hHR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min)*1e-9*120)); //For Noise
							}
                                                }/* else {
                                                        if(isFEBv2r2) {
                                                                if(jj==0) hHR->Fill(strip_numb,1/(last_bc0*(time_max-time_min-(muW2_HR-muW1_HR))*1e-9*1*Strip_length[strip_numb])); //For Noise
                                                        } else {
                                                                if(jj==0) hHR->Fill(m_strip(frame[i]),1/(last_bc0*(time_max-time_min-(muW2_HR-muW1_HR))*1e-9*1*Strip_length[m_strip(frame[i])])); //For Noise
                                                        }
                                                }*/
					}
				}
			} // end of first frame loop


			//if (ntrig_event!=3 and ntrig_event!=0) std::cout<< "nTRIG in this event : "<< ntrig_event  <<std::endl;
			if (!(trig_exists && ntrig_event ==3) ) continue;
			//if (!(trig_exists) ) continue;
			if (bad_trig) continue;
			if (jj==0) cout << "---------- starting event with trigger ------------ " << endl;

			int ntriggerHR_signal_strip = 0;
			int ntriggerLR_signal_strip = 0;

			vector<int> strip_HR;
                        vector<int> strip_LR;
                        vector<double> time_HR;
                        vector<double> time_LR;
			strip_HR.clear();
                        strip_LR.clear();
                        time_HR.clear();
                        time_LR.clear();

                        vector<double> cluster_new_try_strip;
                        vector<double> cluster_new_try_time;
			cluster_new_try_strip.clear();
			cluster_new_try_time.clear();

			for (uint32_t i=0;i<48;i++) {
//				if(i==22) continue;
/*				std::cout<< "LR size: "<< allStrips[i].LRframe.size()<<" HR size: "<<allStrips[i].HRframe.size() << std::endl;
				std::cout<< "-->> HR times: " << std::endl;
                                for (uint32_t fhr=0; fhr<allStrips[i].HRframe.size(); fhr++){
					cout << "allStrips[i].HRtime[fhr]: " << allStrips[i].HRtime[fhr] << endl;
				}
                                std::cout<< "-->> LR times: " << std::endl;
                                for (uint32_t flr=0; flr<allStrips[i].LRframe.size(); flr++){
                                        cout << "allStrips[i].LRtime[flr]: " << allStrips[i].LRtime[flr] << endl;
                                }
*/
				//cout << "muW1_HR: " << muW1_HR << endl;

				int ntriggerHR_signal_per_strip = 0;
				for (uint32_t fhr=0; fhr<allStrips[i].HRframe.size(); fhr++)
				{ // frame HR loop
					//if(i==31) cout << "HR" << endl;
					//cout << "allStrips[i].HRframe.size(): " << allStrips[i].HRframe.size() << endl; 
//					strip_HR.push_back((int)i);
//					time_HR.push_back(allStrips[i].HRtime[fhr]);
                                	double time_in_strip = allStrips[i].HRtime[fhr] - trig[m_fpga(frame[allStrips[i].HRframe[fhr]])] - HR_to_conctor[i]/V_conector;
                                	double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
                                	double time_to_conector_new = HR_to_conctor[i]/V_conector-HR_to_conctor[reference_strip]/V_conector;
                                	double HR_trig = time_in_strip_new + time_to_conector_new;
//					if((HR_trig) > muW1_HR and (HR_trig) < muW2_HR){ //aqui
//                                        	strip_HR.push_back((int)i);
//                                        	time_HR.push_back(HR_trig+trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]);
//                                                time_HR.push_back(time_in_strip_new);
//					}
//					if(jj==0) hHRT->Fill(i,allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]);
//					if(jj==0) hHRT_->Fill(allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]);
                                        if(jj==0) hHRT->Fill(i,HR_trig);
                                        //if(jj==0) hHRT_->Fill(HR_trig);

//                                        if ( !((allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) > muW1_HR and (allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) < muW2_HR ))
					if ( !((HR_trig) > muW1_HR and (HR_trig) < muW2_HR ))
					{
						// HR not in muon window
						allStrips[i].HRframeOK[fhr] = false;
//						if ( ((allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) > muW1_HR-1000-(9*(abs(muW1_HR)-abs(muW2_HR))) and (allStrips[i].HRtime[fhr]-trig[m_fpga(frame[allStrips[i].HRframe[fhr]])]) < muW2_HR-1000 ))
                                                if ( ((HR_trig) > muW1_HR-1000-(9*(abs(muW1_HR)-abs(muW2_HR))) and (HR_trig) < muW2_HR-1000 ))

							//trig = 2000 and 87000, signal is 1000,         
						{
							allStrips[i].HRframeBkg[fhr] = true;
							if(jj==0 && !isFEBv2r2) hLRn->Fill(i,1/(run_duration_in_sec*120)); //For Noise; newly added histogram
                                                        if(jj==0 && isFEBv2r2) hLRn->Fill(i,1/(run_duration_in_sec*120));
						}
					}
					if (allStrips[i].HRframeOK[fhr]){
                                                strip_HR.push_back((int)i);
                                                time_HR.push_back(time_in_strip_new); 
						ntriggerHR_signal++;
						ntriggerHR_signal_per_strip++;
					}
				}
				if(ntriggerHR_signal_per_strip>0) ntriggerHR_signal_strip++;
				if(jj==0) nFiredHR_per_strip->Fill(i,ntriggerHR_signal_per_strip);

				int ntriggerLR_signal_per_strip = 0;
				for (uint32_t flr=0; flr<allStrips[i].LRframe.size(); flr++) { // frame LR loop
					//if(i==31) cout << "LR" << endl;
//                                        strip_LR.push_back((int)i);
//                                        time_LR.push_back(allStrips[i].LRtime[flr]);

//					if(jj==0) hLRT->Fill(i,allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]);
//					if(jj==0) hLRT_->Fill(allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]);
                                        double time_in_strip = allStrips[i].LRtime[flr] - trig[m_fpga(frame[allStrips[i].LRframe[flr]])] - LR_to_conctor[i]/V_conector;
                                        double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
                                        double time_to_conector_new = LR_to_conctor[i]/V_conector-LR_to_conctor[reference_strip]/V_conector;
                                        double LR_trig = time_in_strip_new + time_to_conector_new;
//					if((LR_trig) > muW1_LR and (LR_trig) < muW2_LR){
//                                        	strip_LR.push_back((int)i);
//                                        	time_LR.push_back(LR_trig+trig[m_fpga(frame[allStrips[i].LRframe[flr]])]);
//                                                time_LR.push_back(time_in_strip_new);
//					}
                                        if(jj==0) hLRT->Fill(i,LR_trig);
                                        //if(jj==0) hLRT_->Fill(LR_trig);

//                                        if ( !((allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) > muW1_LR and (allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) < muW2_LR ))
					if ( !((LR_trig) > muW1_LR and (LR_trig) < muW2_LR ))
					{
						// LR not in muon window
						allStrips[i].LRframeOK[flr] = false;
//						if ( ((allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) > muW1_LR-1000-(9*(abs(muW1_LR)-abs(muW2_LR))) and (allStrips[i].LRtime[flr]-trig[m_fpga(frame[allStrips[i].LRframe[flr]])]) < muW2_LR-1000 ))
                                                if ( ((LR_trig) > muW1_LR-1000-(9*(abs(muW1_LR)-abs(muW2_LR))) and (LR_trig) < muW2_LR-1000 ))
						{ 
							allStrips[i].LRframeBkg[flr] = true;
                                                        if(jj==0 && !isFEBv2r2) hHRn->Fill(i,1/(run_duration_in_sec*120)); //For Noise; newly added histogram
							if(jj==0 && isFEBv2r2) hHRn->Fill(i,1/(run_duration_in_sec*120));

						}   //THIS AND THE PREVIOUS FOUR LINES ADDED BY ME
					}
					if (allStrips[i].LRframeOK[flr]){
                                                strip_LR.push_back((int)i);
                                                time_LR.push_back(time_in_strip_new);
						ntriggerLR_signal++;
						ntriggerLR_signal_per_strip++;
					}

				} // LR Frame loop ends
				if(ntriggerLR_signal_per_strip>0) ntriggerLR_signal_strip++;
				if(jj==0) nFiredLR_per_strip->Fill(i,ntriggerLR_signal_per_strip);

				allStrips[i].buildDeltas();
				//std::cout<< "how many delta computed : "<<allStrips[i].deltas.size() << std::endl;
				if (allStrips[i].deltas.size() > 0 )
				{
					has_paired_strips = 1; //In the event there is at least 1 paired strip
				}
				nPaired_per_strip->Fill(i,allStrips[i].deltas.size());
				for (uint32_t d=0; d<allStrips[i].deltas.size(); d++)
				{
                                        n_paired_srip++;

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

//					cout << "deltaT: " <<  (m_time(m_traw(frame[allStrips[i].deltas[d].first]))-m_time(m_traw(frame[allStrips[i].deltas[d].second]))) << endl;
					

					deltaT=time_in_strip2_new-time_in_strip_new+Strip_length[reference_strip]/(V_conector);
					deltaT=deltaT/2.;
					//m_time(m_traw(frame[allStrips[i].deltas[d].first])) - time_corr_LR  -  m_time(m_traw(frame[allStrips[i].deltas[d].second])) + time_corr_HR;//LR_trig-HR_trig;//time_in_strip2_new-time_in_strip_new+Strip_length[reference_strip]/(V_conector); //LR_trig-HR_trig; TRABALHAR AQUI
/*					if(i==24 || i==25 || i==26){
						cout << "i: " << i << "; deltaT: " << deltaT << "; HR: " << m_time(m_traw(frame[allStrips[i].deltas[d].second])) << "; LR: " << m_time(m_traw(frame[allStrips[i].deltas[d].first])) << endl;
					}
*/
					if(deltaT>0.) {cluster_new_try_time.push_back(deltaT); cluster_new_try_strip.push_back(i);} 

					//cout<<"                deltaT: " << deltaT << endl;
					if(deltaT>0/*minimum_time_diff[reference_strip]*/) n_paired_srip_med++;

					sumT=time_in_strip2_new+time_in_strip_new;//LR_trig+HR_trig;
					//if(jj==0) cout << "sumT: " << sumT << ";   trig time: " << trig_time << endl;
					if(jj==0) hdeltaT->Fill(i,deltaT);
					if(deltaT>0 && jj==0) hsumT->Fill(i,sumT);
					if(jj==0) hdeltaT_1D->Fill(deltaT);
					//cout << "strip: " << i << ";  sum time: " << sumT << endl;
					if(deltaT>0 && jj==0) hsumT_1D->Fill(sumT);
					if(deltaT>0 && jj==0) hsumTdeltaT_2D->Fill(sumT,deltaT);
					if (deltaT>0./*=minimum_time_diff[reference_strip]*/) {
						medium_num = 1; 
//						cout << "strip: " << i << endl;
					}
					if (deltaT > 8 ) {tight_num = 1; }
					if(deltaT > 0 && jj==0){
						hHRT_med->Fill(i,/*allStrips[i].HRtime[fhr]m_time(m_traw(frame[allStrips[i].deltas[d].second]))-trig_time*/HR_trig);
						hHRT__med->Fill(/*allStrips[i].HRtime[fhr]m_time(m_traw(frame[allStrips[i].deltas[d].second]))-trig_time*/HR_trig);
						hLRT_med->Fill(i,/*allStrips[i].LRtime[flr]m_time(m_traw(frame[allStrips[i].deltas[d].first]))-trig_time*/LR_trig);
						hLRT__med->Fill(/*allStrips[i].LRtime[flr]m_time(m_traw(frame[allStrips[i].deltas[d].first]))-trig_time*/LR_trig);
						hHRT_temp_med->Fill(/*allStrips_temp[i].HRtime[fhr]m_time(m_traw(frame[allStrips[i].deltas[d].second]))-trig_time*/HR_trig);
						hLRT_temp_med->Fill(/*allStrips_temp[i].LRtime[flr]m_time(m_traw(frame[allStrips[i].deltas[d].first]))-trig_time*/LR_trig);
						hsumT_1D_med->Fill(sumT);
					}

				}
			} // strip loop ends
			bool cluster_good = false;
                        bool cluster_good_2 = false;

			for(int iii=0;iii<strip_HR.size();iii++){
				if(iii==(strip_HR.size()-1)) continue;
                                if(iii==0) continue;
				hHR_check_xtalk_temp->Fill(strip_HR.at(iii));
				if(((strip_HR.at(iii)-strip_HR.at(iii+1))==-1 && abs(time_HR.at(iii)-time_HR.at(iii+1))<3) && ((strip_HR.at(iii)-strip_HR.at(iii-1))==1 && abs(time_HR.at(iii)-time_HR.at(iii-1))<3)){
					hHR_check_xtalk->Fill(strip_HR.at(iii));

				}
			}

                        for(int iii=0;iii<strip_LR.size();iii++){
                                if(iii==(strip_LR.size()-1)) break;
                                if(iii==0) continue;
                                hLR_check_xtalk_temp->Fill(strip_LR.at(iii));
                                if(((strip_LR.at(iii)-strip_LR.at(iii+1))==-1 && abs(time_LR.at(iii)-time_LR.at(iii+1))<3) && ((strip_LR.at(iii)-strip_LR.at(iii-1))==1 && abs(time_LR.at(iii)-time_LR.at(iii-1))<3)){
                                        hLR_check_xtalk->Fill(strip_LR.at(iii));
                                }
                        }


			if(medium_num){
                        	double limit = 0.1;
				for(int ij=0;ij<100;ij++){
//					vector<std::pair<vector<int>,vector<double>>> cluster_HR = cluster_one_side(strip_HR, time_HR, limit);
//		                	vector<std::pair<vector<int>,vector<double>>> cluster_LR = cluster_one_side(strip_LR, time_LR, limit);
                                        vector<ClusterOneSide> cluster_HR = cluster_one_side(strip_HR, time_HR, limit);
                                        vector<ClusterOneSide> cluster_LR = cluster_one_side(strip_LR, time_LR, limit);

					std::vector<cluster> CLUSTER;
					CLUSTER.clear();
					if(cluster_HR.size()>0 && cluster_LR.size()>0) CLUSTER = final_clustering(cluster_LR,cluster_HR,Strip_length[reference_strip]/(V_conector));
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


			// cluster new try

                        std::vector<cluster_new_try> CLUSTER_2;
                        CLUSTER_2.clear();
			if(cluster_new_try_strip.size()>0) CLUSTER_2 = cluster_new_try_func(cluster_new_try_strip,cluster_new_try_time);




//                        vector<std::pair<vector<int>,vector<double>>> cluster_HR_new = cluster_one_side(strip_HR, time_HR, 3.);
//                        vector<std::pair<vector<int>,vector<double>>> cluster_LR_new = cluster_one_side(strip_LR, time_LR, 3.);
			//cout << "HR: " << endl;
                        vector<ClusterOneSide> cluster_HR_new = cluster_one_side(strip_HR, time_HR, 3.);
                        //cout << "LR: " << endl;
                        vector<ClusterOneSide> cluster_LR_new = cluster_one_side(strip_LR, time_LR, 3.);
                        std::vector<cluster> CLUSTER;
                        CLUSTER.clear();
//                        if(cluster_HR_new.size()>0 && cluster_LR_new.size()>0) CLUSTER = final_clustering(cluster_LR_new,cluster_HR_new,Strip_length[reference_strip]/(V_conector));
                        if(cluster_HR_new.size()>0 && cluster_LR_new.size()>0){
				//cout << "cluster_HR_new.size(): " << cluster_LR_new.size() << ";  cluster_LR_new.size(): " << cluster_LR_new.size() << endl;
                                CLUSTER = final_clustering(cluster_LR_new,cluster_HR_new,Strip_length[reference_strip]/(V_conector));
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
                                //cout << "strip_diff abs: " << strip_diff << ";   strip_diff: " << cluster_HR_new.at(HR_ind_smallest).strip() - cluster_LR_new.at(LR_ind_smallest).strip()  << endl;
				strip_diff_h->Fill(cluster_HR_new.at(HR_ind_smallest).strip() - cluster_LR_new.at(LR_ind_smallest).strip());
/*                        	if(CLUSTER.size()==0){ 
					cout << "------------------>>>> strip_HR.size(): " << strip_HR.size() << ";   strip_LR.size(): " <<  strip_LR.size() << endl;
					cout << "cluster_HR_new.size(): " << cluster_HR_new.size() << ";  cluster_LR_new.size(): " << cluster_LR_new.size() << endl;

					cout << "     HR strips: ";
					for(int ij=0;ij<cluster_HR_new.size();ij++){ 
						for(int ji=0;ji<cluster_HR_new.at(ij).strips_v().size();ji++){
							cout <<  cluster_HR_new.at(ij).strips_v().at(ji) << ",";
						}	
						cout << ";  ";
					}
					cout << endl;
                                        cout << endl;
                                        cout << "     LR strips: ";
                                        for(int ij=0;ij<cluster_LR_new.size();ij++){
                                                for(int ji=0;ji<cluster_LR_new.at(ij).strips_v().size();ji++){
                                                        cout <<  cluster_LR_new.at(ij).strips_v().at(ji) << ",";
                                                }
                                                cout << ";  ";
                                        }
                                        cout << endl;
                                        cout << endl;

					cout << "smallest diff: " << cluster_HR_new.at(HR_ind_smallest).strip() - cluster_LR_new.at(LR_ind_smallest).strip() << endl;
				}*/
                        }

/*			if(medium_num!=0 && CLUSTER.size()==0){
				cout << "--->>> medium_num: " << medium_num << "; cluster_HR_new.size(): " <<  cluster_HR_new.size() << "; cluster_LR_new.size(): " << cluster_LR_new.size() << endl;
				std::vector<cluster> CLUSTER_new;
				CLUSTER_new = final_clustering(cluster_LR_new,cluster_HR_new,true);
			}
*/
			if(jj==0) hNClusters->Fill(CLUSTER.size());
			sum_cluster_size = sum_cluster_size + CLUSTER.size();
			for(int ij=0;ij<CLUSTER.size();ij++){
				if(jj==0) hsumTCluster_1D->Fill(CLUSTER.at(ij).sum_time());
                                if(CLUSTER.at(ij).time()>0.) cluster_good=true;
				if(jj==0) {hClusterSize->Fill(CLUSTER.at(ij).size());}//; cout << "CLUSTER.at(ij).size(): " << CLUSTER.at(ij).size() << endl;}
				if(jj==0) hdeltaTCluster_1D->Fill(CLUSTER.at(ij).time());
				if(jj==0) hdeltaClusterT->Fill(CLUSTER.at(ij).strip(),(CLUSTER.at(ij).time()));
				if(CLUSTER.at(ij).time()>0.&&jj==0) {
					number_of_good_clusters++;
					if(CLUSTER.at(ij).HR_LR_strips_ratio()==1) number_of_good_clusters_ratio1++;
				}
			}
			if(cluster_good) N_cluster_good++;

/*			if(!cluster_good){
				cout << "NOT EFFICIENT FOR CLUSTERING!!!" << endl;
				cout << "CLUSTER.size(): " << CLUSTER.size() << endl;
                        	for(int ij=0;ij<CLUSTER.size();ij++){
					cout << "CLUSTER.at(ij).time(): " <<CLUSTER.at(ij).time() << ";  CLUSTER.at(ij).strip(): " << CLUSTER.at(ij).strip() << endl;
                        	}
			}
*/

			if(jj==0) hNClusters_2->Fill(CLUSTER_2.size());
			sum_cluster2_size = sum_cluster2_size + CLUSTER_2.size();
			for(int ij=0;ij<CLUSTER_2.size();ij++){
                                if(CLUSTER_2.at(ij).time()>0.) cluster_good_2=true;
				if(jj==0) hClusterSize_2->Fill(CLUSTER_2.at(ij).size());
				if(jj==0) hdeltaTCluster_2_1D->Fill(CLUSTER_2.at(ij).time());
				if(jj==0) hdeltaClusterT_2->Fill(CLUSTER_2.at(ij).strip(),(CLUSTER_2.at(ij).time()));
			}
			if(cluster_good_2) N_cluster_good_2++;




			if(medium_num && jj==0) {
				hnPairs_med->Fill(n_paired_srip_med);
                                hHRT_2_med->Fill(ntriggerHR_signal_strip);
                                hLRT_2_med->Fill(ntriggerLR_signal_strip);
			}

			hHRT_temp_med->GetXaxis()->SetRangeUser(hHRT_temp_med->GetXaxis()->GetBinCenter(hHRT_temp_med->GetMaximumBin())-50,hHRT_temp_med->GetXaxis()->GetBinCenter(hHRT_temp_med->GetMaximumBin())+50);
			hLRT_temp_med->GetXaxis()->SetRangeUser(hLRT_temp_med->GetXaxis()->GetBinCenter(hLRT_temp_med->GetMaximumBin())-50,hLRT_temp_med->GetXaxis()->GetBinCenter(hLRT_temp_med->GetMaximumBin())+50);


			if((ntriggerHR_signal>0) || (ntriggerLR_signal>0)) nANDstrip_noCut_fired++;
			if(ntriggerLR_signal>0) nANDstripLR_noCut_fired++;
			if(ntriggerHR_signal>0) nANDstripHR_noCut_fired++;

			ntrig_allevent++; // count denom of eff
			if (has_paired_strips) {
				nANDstrip_fired++;
				if(jj==0) hHRT_2->Fill(ntriggerHR_signal_strip);
				if(jj==0) hLRT_2->Fill(ntriggerLR_signal_strip);
			}
			if (medium_num) {nANDstrip_medium_fired++;}// else { cout << "===================================================================" << endl;}
			if (tight_num) {nANDstrip_tight_fired++;}
			if (jj==0) hnPairs->Fill(n_paired_srip);
			//std::cout << "ntrig_allevent: " << ntrig_allevent << std::endl;
		} //Event loop ends


		cout << "---->> number_of_good_clusters: " << number_of_good_clusters << "; number_of_good_clusters_ratio1: " <<  number_of_good_clusters_ratio1 << endl;

		if(jj==0 && ClusterE!=0){
                	for(int ij=0;ij<100;ij++){
				ClusterS.at(ij) = ClusterS.at(ij)/ClusterE;
                        	ClusterN.at(ij) = ClusterN.at(ij)/ClusterE;
				//cout << "limit: " << ClusterL.at(ij) << ";  Number of Clusters: " << ClusterN.at(ij) << ";  Cluster size: " << ClusterS.at(ij) << endl;
				g_cluster_size->SetPoint(ij,ClusterL.at(ij),ClusterS.at(ij));
                                g_cluster_number->SetPoint(ij,ClusterL.at(ij),ClusterN.at(ij));
			}
		}

		myfile << (nANDstrip_noCut_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstrip_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstrip_medium_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstrip_tight_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstripHR_noCut_fired)/float(ntrig_allevent) <<"\n";
		myfile << (nANDstripLR_noCut_fired)/float(ntrig_allevent) <<"\n";
                myfile << (N_cluster_good)/float(ntrig_allevent) <<"\n";
                myfile << (N_cluster_good_2)/float(ntrig_allevent) <<"\n";
//		myfile << sum_cluster_size/(ntrig_allevent*40*pow(10,-9)*6000*cluster_eff) <<"\n";
//              myfile << sum_cluster2_size/(ntrig_allevent*40*pow(10,-9)*6000*cluster2_eff) <<"\n";
                myfile << sum_cluster_size/(ntrig_allevent) <<"\n";
                myfile << sum_cluster2_size/(ntrig_allevent) <<"\n";

		if(jj==0){
			cluster_eff = ((N_cluster_good)/float(ntrig_allevent));
			cluster2_eff = ((N_cluster_good_2)/float(ntrig_allevent));
		}

		myfile<<""<<"\n";

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
	hLR->GetYaxis()->SetTitle("Noise Hz/cm2");
	hLR->GetXaxis()->SetTitle("Strips");
	pad2->cd();
	pad2->SetGrid();
	pad2->SetLogy();
	hHR->Draw("Hist");
	hHR->GetYaxis()->SetTitle("Noise Hz/cm2");
	hHR->GetXaxis()->SetTitle("Strips");
	c1->cd();
	c1->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig__avarage_noise_side.root",hv_,sn_/*,mt_*/);
	c1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig__avarage_noise_side.pdf",hv_,sn_/*,mt_*/);
	c1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig__avarage_noise_side.png",hv_,sn_/*,mt_*/);
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
	newpad1->SetLogy();
	hLRn->Scale(1./ntrig_allevent);
	hLRn->Draw("Hist");
	hLRn->GetYaxis()->SetTitle("Noise Hz/cm");
	hLRn->GetXaxis()->SetTitle("Strips");
	newpad2->cd();
	newpad2->SetGrid();
	//newpad2->SetLogy();
	hHRn->Scale(1./ntrig_allevent);
	hHRn->Draw("Hist");
	hHRn->GetYaxis()->SetTitle("Noise Hz/cm");
	hHRn->GetXaxis()->SetTitle("Strips");
	can1->cd();
	can1->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_new_noise_histo.root",hv_,sn_/*,mt_*/);
	//can1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_new_noise_histo.pdf",hv_,sn_/*,mt_*/);
	//can1->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_new_noise_histo.png",hv_,sn_/*,mt_*/);
	//can1->SaveAs(s);

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
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_srip.root",hv_,sn_/*,mt_*/);
	c2->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_srip.pdf",hv_,sn_/*,mt_*/);
	c2->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_srip.png",hv_,sn_/*,mt_*/);
	c2->SaveAs(s);




	hdeltaT->GetYaxis()->SetRangeUser(2,30);
        hdeltaT->GetXaxis()->SetRangeUser(18,37);
        TCanvas *c2P = new TCanvas("c4P","The Ntuple canvas",200,10,900,780);
        TPad *pad4P = new TPad("pad4P","This is pad4",0.02,0.02,0.98,0.98,21);
        pad4P->UseCurrentStyle();
        pad4P->Draw();
        pad4P->cd();
        pad4P->SetRightMargin(0.12);
   	TProfile *prof = hdeltaT->ProfileX();
   	prof->Fit("pol1");
	prof->Draw();
	c2P->cd();
        c2P->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_srip_profile.png",hv_,sn_/*,mt_*/);
        //c2P->SaveAs(s);

	if(true){
        TCanvas *c2C = new TCanvas("c4C","The Ntuple canvas",200,10,900,780);
        TPad *pad4C = new TPad("pad4C","This is pad4C",0.02,0.02,0.98,0.98,21);
        pad4C->UseCurrentStyle();
        pad4C->Draw();
        pad4C->cd();
        pad4C->SetGrid();
        pad4C->SetRightMargin(0.12);
        hdeltaClusterT->GetYaxis()->SetTitle("Cluster delta time [ns]");
        hdeltaClusterT->GetXaxis()->SetTitle("Cluster Strips");
        hdeltaClusterT->GetYaxis()->SetLabelOffset(0.007);
        hdeltaClusterT->GetYaxis()->SetTitleOffset(1.5);
        hdeltaClusterT->SetLabelSize(0.02,"x");
        hdeltaClusterT->SetTitle("");
        hdeltaClusterT->Draw("P");
        c2C->cd();
        c2C->Update();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_deltaT_srip.root",hv_,sn_/*,mt_*/);
        c2C->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_deltaT_srip.pdf",hv_,sn_/*,mt_*/);
        c2C->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_deltaT_srip.png",hv_,sn_/*,mt_*/);
        c2C->SaveAs(s);
	}

	if(true){
        TCanvas *c2C = new TCanvas("c4C","The Ntuple canvas",200,10,900,780);
        TPad *pad4C = new TPad("pad4C","This is pad4C",0.02,0.02,0.98,0.98,21);
        pad4C->UseCurrentStyle();
        pad4C->Draw();
        pad4C->cd();
        pad4C->SetGrid();
        pad4C->SetRightMargin(0.12);
        hdeltaClusterT_2->GetYaxis()->SetTitle("Cluster2 delta time [ns]");
        hdeltaClusterT_2->GetXaxis()->SetTitle("Cluster2 Strips");
        hdeltaClusterT_2->GetYaxis()->SetLabelOffset(0.007);
        hdeltaClusterT_2->GetYaxis()->SetTitleOffset(1.5);
        hdeltaClusterT_2->SetLabelSize(0.02,"x");
        hdeltaClusterT_2->SetTitle("");
        hdeltaClusterT_2->Draw("P");
        c2C->cd();
        c2C->Update();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_deltaT_srip.root",hv_,sn_/*,mt_*/);
        //c2C->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_deltaT_srip.pdf",hv_,sn_/*,mt_*/);
        //c2C->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_deltaT_srip.png",hv_,sn_/*,mt_*/);
        //c2C->SaveAs(s);
	}



	TCanvas *c22 = new TCanvas("c42","The Ntuple canvas",200,10,900,780);
	TPad *pad42 = new TPad("pad42","This is pad42",0.02,0.02,0.98,0.98,21);
	pad42->UseCurrentStyle();
	pad42->Draw();
	pad42->cd();
	pad42->SetGrid();
	pad42->SetRightMargin(0.12);
	hsumT->GetYaxis()->SetTitle("sum time [ns]");
	hsumT->GetXaxis()->SetTitle("Strips");
	hsumT->GetYaxis()->SetLabelOffset(0.007);
	hsumT->GetYaxis()->SetTitleOffset(1.5);
	hsumT->SetLabelSize(0.02,"x");
	hsumT->SetTitle("");
	hsumT->Draw("P");
	c22->cd();
	c22->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_srip.root",hv_,sn_/*,mt_*/);
	//c22->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_srip.pdf",hv_,sn_/*,mt_*/);
	//c22->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_srip.png",hv_,sn_/*,mt_*/);
	//c22->SaveAs(s);



	if(true){
	TCanvas *c22C = new TCanvas("c42C","The Ntuple canvas",200,10,900,780);
	TPad *pad42C = new TPad("pad42C","This is pad42",0.02,0.02,0.98,0.98,21);
	pad42C->UseCurrentStyle();
	pad42C->Draw();
	pad42C->cd();
	pad42C->SetRightMargin(0.12);
	hsumTCluster_1D->GetXaxis()->SetTitle("Cluster sum time [ns]");
	hsumTCluster_1D->GetYaxis()->SetTitleOffset(1.5);
	hsumTCluster_1D->SetLabelSize(0.02,"x");
	hsumTCluster_1D->SetTitle("");
	hsumTCluster_1D->Draw();
	c22C->cd();
	c22C->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_sumT_srip.root",hv_,sn_/*,mt_*/);
	c22C->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_sumT_srip.pdf",hv_,sn_/*,mt_*/);
	c22C->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_sumT_srip.png",hv_,sn_/*,mt_*/);
	c22C->SaveAs(s);
	}


	if(true){
	TCanvas *c22C = new TCanvas("c42C","The Ntuple canvas",200,10,900,780);
	TPad *pad42C = new TPad("pad42C","This is pad42",0.02,0.02,0.98,0.98,21);
	pad42C->UseCurrentStyle();
	pad42C->Draw();
	pad42C->cd();
	pad42C->SetGrid();
	pad42C->SetRightMargin(0.12);
	hsumTCluster_2_1D->GetYaxis()->SetTitle("Cluster2 sum time [ns]");
	hsumTCluster_2_1D->GetXaxis()->SetTitle("Cluster2 Strips");
	hsumTCluster_2_1D->GetYaxis()->SetLabelOffset(0.007);
	hsumTCluster_2_1D->GetYaxis()->SetTitleOffset(1.5);
	hsumTCluster_2_1D->SetLabelSize(0.02,"x");
	hsumTCluster_2_1D->SetTitle("");
	hsumTCluster_2_1D->Draw("P");
	c22C->cd();
	c22C->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_sumT_srip.root",hv_,sn_/*,mt_*/);
	//c22C->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_sumT_srip.pdf",hv_,sn_/*,mt_*/);
	//c22C->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_sumT_srip.png",hv_,sn_/*,mt_*/);
	//c22C->SaveAs(s);
	}

	gStyle->SetOptStat(1111);
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
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_n_paired_srip.root",hv_,sn_/*,mt_*/);
	cPairs->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_n_paired_srip.pdf",hv_,sn_/*,mt_*/);
	cPairs->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_n_paired_srip.png",hv_,sn_/*,mt_*/);
	cPairs->SaveAs(s);    
        gStyle->SetOptStat(0);

	if(true){
        gStyle->SetOptStat(1111);
	TCanvas *cPairsC = new TCanvas("cPairsC","The Ntuple canvas",200,10,900,780);
	TPad *padPC = new TPad("padPC","This is padP",0.02,0.02,0.98,0.98,21);
	padPC->UseCurrentStyle();
	padPC->Draw();
	padPC->cd();
	padPC->SetGrid();
	padPC->SetRightMargin(0.12);
	hClusterSize->GetYaxis()->SetTitle("Events");
	hClusterSize->GetXaxis()->SetTitle("Cluster size");
	hClusterSize->GetYaxis()->SetLabelOffset(0.007);
	hClusterSize->GetYaxis()->SetTitleOffset(1.5);
	hClusterSize->SetLabelSize(0.02,"x");
	hClusterSize->SetTitle("");
	hClusterSize->Draw("Histo");
	cPairsC->cd();
	cPairsC->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_size.root",hv_,sn_/*,mt_*/);
	cPairsC->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_size.pdf",hv_,sn_/*,mt_*/);
	cPairsC->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_size.png",hv_,sn_/*,mt_*/);
	cPairsC->SaveAs(s);    
        gStyle->SetOptStat(0);
	}


        if(true){
        TCanvas *cPairsC = new TCanvas("cPairsC","The Ntuple canvas",200,10,900,780);
        TPad *padPC = new TPad("padPC","This is padP",0.02,0.02,0.98,0.98,21);
        padPC->UseCurrentStyle();
        padPC->Draw();
        padPC->cd();
        padPC->SetGrid();
        padPC->SetRightMargin(0.12);
        hClusterSize_2->GetYaxis()->SetTitle("Events");
        hClusterSize_2->GetXaxis()->SetTitle("Cluster2 size");
        hClusterSize_2->GetYaxis()->SetLabelOffset(0.007);
        hClusterSize_2->GetYaxis()->SetTitleOffset(1.5);
        hClusterSize_2->SetLabelSize(0.02,"x");
        hClusterSize_2->SetTitle("");
        hClusterSize_2->Draw("Histo");
        cPairsC->cd();
        cPairsC->Update();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_size.root",hv_,sn_/*,mt_*/);
        //cPairsC->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_size.pdf",hv_,sn_/*,mt_*/);
        //cPairsC->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_size.png",hv_,sn_/*,mt_*/);
        //cPairsC->SaveAs(s);
        }




	if(true){
        gStyle->SetOptStat(1111);
	TCanvas *cPairsCC = new TCanvas("cPairsCC","The Ntuple canvas",200,10,900,780);
	TPad *padPCC = new TPad("padPCC","This is padP",0.02,0.02,0.98,0.98,21);
	padPCC->UseCurrentStyle();
	padPCC->Draw();
	padPCC->cd();
	padPCC->SetGrid();
	padPCC->SetRightMargin(0.12);
	hNClusters->GetYaxis()->SetTitle("Events");
	hNClusters->GetXaxis()->SetTitle("Number of Clusters");
	hNClusters->GetYaxis()->SetLabelOffset(0.007);
	hNClusters->GetYaxis()->SetTitleOffset(1.5);
	hNClusters->SetLabelSize(0.02,"x");
	hNClusters->SetTitle("");
	hNClusters->Draw("Histo");
	cPairsCC->cd();
	cPairsCC->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters.root",hv_,sn_/*,mt_*/);
	cPairsCC->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters.pdf",hv_,sn_/*,mt_*/);
	cPairsCC->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters.png",hv_,sn_/*,mt_*/);
	cPairsCC->SaveAs(s);    
        gStyle->SetOptStat(0);
	}


	if(true){
	TCanvas *cPairsCC = new TCanvas("cPairsCC2","The Ntuple canvas",200,10,900,780);
	TPad *padPCC = new TPad("padPCC","This is padP",0.02,0.02,0.98,0.98,21);
	padPCC->UseCurrentStyle();
	padPCC->Draw();
	padPCC->cd();
	padPCC->SetGrid();
	padPCC->SetRightMargin(0.12);
	hNClusters_2->GetYaxis()->SetTitle("Events");
	hNClusters_2->GetXaxis()->SetTitle("Number of Clusters2");
	hNClusters_2->GetYaxis()->SetLabelOffset(0.007);
	hNClusters_2->GetYaxis()->SetTitleOffset(1.5);
	hNClusters_2->SetLabelSize(0.02,"x");
	hNClusters_2->SetTitle("");
	hNClusters_2->Draw("Histo");
	cPairsCC->cd();
	cPairsCC->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters2.root",hv_,sn_/*,mt_*/);
	//cPairsCC->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters2.pdf",hv_,sn_/*,mt_*/);
	//cPairsCC->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters2.png",hv_,sn_/*,mt_*/);
	//cPairsCC->SaveAs(s);    
	}



        TCanvas *cPairsm = new TCanvas("cPairsm","The Ntuple canvas",200,10,900,780);
        TPad *padPm = new TPad("padPm","This is padP",0.02,0.02,0.98,0.98,21);
        padPm->UseCurrentStyle();
        padPm->Draw();
        padPm->cd();
        padPm->SetGrid();
        padPm->SetRightMargin(0.12);
        hnPairs_med->GetYaxis()->SetTitle("Events");
        hnPairs_med->GetXaxis()->SetTitle("n paired strips");
        hnPairs_med->GetYaxis()->SetLabelOffset(0.007);
        hnPairs_med->GetYaxis()->SetTitleOffset(1.5);
        hnPairs_med->SetLabelSize(0.02,"x");
        hnPairs_med->SetTitle("");
        hnPairs_med->Draw("Histo");
        cPairsm->cd();
        cPairsm->Update();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_n_paired_srip_med.root",hv_,sn_/*,mt_*/);
        //cPairsm->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_n_paired_srip_med.pdf",hv_,sn_/*,mt_*/);
        //cPairsm->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_n_paired_srip_med.png",hv_,sn_/*,mt_*/);
        //cPairsm->SaveAs(s);


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
	hdeltaT_1D->Draw("Histo");
	chdeltaT_1D->cd();
	chdeltaT_1D->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_1D_srip.root",hv_,sn_/*,mt_*/);
	chdeltaT_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_1D_srip.pdf",hv_,sn_/*,mt_*/);
	chdeltaT_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_deltaT_1D_srip.png",hv_,sn_/*,mt_*/);
	chdeltaT_1D->SaveAs(s);


	if(true){
	TCanvas *chdeltaTC_1D = new TCanvas("chdeltaTC_1D","The Ntuple canvas",200,10,900,780);
	TPad *padhdeltaTC_1D = new TPad("padhdeltaTC_1D","This is hdeltaT_1D",0.02,0.02,0.98,0.98,21);
	padhdeltaTC_1D->UseCurrentStyle();
	padhdeltaTC_1D->Draw();
	padhdeltaTC_1D->cd();
	padhdeltaTC_1D->SetGrid();
	padhdeltaTC_1D->SetRightMargin(0.12);
	hdeltaTCluster_1D->GetXaxis()->SetTitle("Cluster delta time [ns]");
	hdeltaTCluster_1D->GetYaxis()->SetTitle("Events");
	hdeltaTCluster_1D->GetYaxis()->SetLabelOffset(0.007);
	hdeltaTCluster_1D->GetYaxis()->SetTitleOffset(1.5);
	hdeltaTCluster_1D->SetLabelSize(0.02,"x");
	hdeltaTCluster_1D->SetTitle("");
	hdeltaTCluster_1D->Draw("Histo");
	chdeltaTC_1D->cd();
	chdeltaTC_1D->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_deltaT_1D_srip.root",hv_,sn_/*,mt_*/);
	chdeltaTC_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_deltaT_1D_srip.pdf",hv_,sn_/*,mt_*/);
	chdeltaTC_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_deltaT_1D_srip.png",hv_,sn_/*,mt_*/);
	chdeltaTC_1D->SaveAs(s);
	}

	if(true){
	TCanvas *chdeltaTC_1D = new TCanvas("chdeltaTC_1D","The Ntuple canvas",200,10,900,780);
	TPad *padhdeltaTC_1D = new TPad("padhdeltaTC_1D","This is hdeltaT_1D",0.02,0.02,0.98,0.98,21);
	padhdeltaTC_1D->UseCurrentStyle();
	padhdeltaTC_1D->Draw();
	padhdeltaTC_1D->cd();
	padhdeltaTC_1D->SetGrid();
	padhdeltaTC_1D->SetRightMargin(0.12);
	hdeltaTCluster_2_1D->GetXaxis()->SetTitle("Cluster2 delta time [ns]");
	hdeltaTCluster_2_1D->GetYaxis()->SetTitle("Events");
	hdeltaTCluster_2_1D->GetYaxis()->SetLabelOffset(0.007);
	hdeltaTCluster_2_1D->GetYaxis()->SetTitleOffset(1.5);
	hdeltaTCluster_2_1D->SetLabelSize(0.02,"x");
	hdeltaTCluster_2_1D->SetTitle("");
	hdeltaTCluster_2_1D->Draw("Histo");
	chdeltaTC_1D->cd();
	chdeltaTC_1D->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_deltaT_1D_srip.root",hv_,sn_/*,mt_*/);
	//chdeltaTC_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_deltaT_1D_srip.pdf",hv_,sn_/*,mt_*/);
	//chdeltaTC_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_Cluster2_deltaT_1D_srip.png",hv_,sn_/*,mt_*/);
	//chdeltaTC_1D->SaveAs(s);
	}


	TCanvas *chsumT_1D = new TCanvas("chsumT_1D","",200,10,900,780);
	chsumT_1D->cd();
	hsumT_1D->GetXaxis()->SetTitle("sum time [ns]");
	hsumT_1D->GetYaxis()->SetTitle("Events");
	hsumT_1D->GetYaxis()->SetLabelOffset(0.007);
	hsumT_1D->GetYaxis()->SetTitleOffset(1.5);
	hsumT_1D->SetLabelSize(0.02,"x");
	hsumT_1D->SetTitle("");
	//hsumT_1D->Scale(1/(hsumT_1D->Integral(hsumT_1D->FindBin(9),hsumT_1D->FindBin(20))));
	hsumT_1D->Draw("Histo");
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_1D_srip.root",hv_,sn_/*,mt_*/);
	chsumT_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_1D_srip.pdf",hv_,sn_/*,mt_*/);
	chsumT_1D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_1D_srip.png",hv_,sn_/*,mt_*/);
	chsumT_1D->SaveAs(s);


        TCanvas *chsumT_1Dm = new TCanvas("chsumT_1Dm","",200,10,900,780);
        chsumT_1Dm->cd();
        hsumT_1D_med->GetXaxis()->SetTitle("sum time [ns]");
        hsumT_1D_med->GetYaxis()->SetTitle("Events");
        hsumT_1D_med->GetYaxis()->SetLabelOffset(0.007);
        hsumT_1D_med->GetYaxis()->SetTitleOffset(1.5);
        hsumT_1D_med->SetLabelSize(0.02,"x");
        hsumT_1D_med->SetTitle("");
        hsumT_1D_med->Draw("Histo");
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_1D_srip_med.root",hv_,sn_/*,mt_*/);
       //chsumT_1Dm->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_1D_srip_med.pdf",hv_,sn_/*,mt_*/);
        //chsumT_1Dm->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_1D_srip_med.png",hv_,sn_/*,mt_*/);
        //chsumT_1Dm->SaveAs(s);



	TCanvas *csumTdeltaT_2D = new TCanvas("csumTdeltaT_2D","",200,10,780,780);
	csumTdeltaT_2D->cd();
	hsumTdeltaT_2D->Draw("COLZ");
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_deltaT.root",hv_,sn_/*,mt_*/);
	//csumTdeltaT_2D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_deltaT.pdf",hv_,sn_/*,mt_*/);
	//csumTdeltaT_2D->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_sumT_deltaT.png",hv_,sn_/*,mt_*/);
	//csumTdeltaT_2D->SaveAs(s);



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
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side.root",hv_,sn_/*,mt_*/);
	c3->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side.pdf",hv_,sn_/*,mt_*/);
	c3->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side.png",hv_,sn_/*,mt_*/);
	c3->SaveAs(s);

	TCanvas *c3m = new TCanvas("c3m","The Ntuple canvas",200,10,780,780);
	TPad *pad3m = new TPad("pad3m","This is pad3",0.05,0.05,0.98,0.55,21);
	TPad *pad5m = new TPad("pad5m","This is pad5",0.05,0.55,0.98,0.98,21);
	pad3m->UseCurrentStyle(); pad5m->UseCurrentStyle();
	pad3m->Draw(); pad5m->Draw();
	pad3m->cd();
	pad3m->SetGrid();
	hLRT_med->Draw("COLZ");
	pad5m->cd();
	pad5m->SetGrid();
	hHRT_med->Draw("COLZ");
	c3m->cd();
	c3m->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_med.root",hv_,sn_/*,mt_*/);
	//c3m->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_med.pdf",hv_,sn_/*,mt_*/);
	//c3m->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_med.png",hv_,sn_/*,mt_*/);
	//c3m->SaveAs(s);


	TCanvas *c3n = new TCanvas("c3n","The Ntuple canvas",200,10,780,780);
	c3n->cd();
	c3n->SetGrid();
	nPaired_per_strip->Draw("COLZ");
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_npaired_per_strip.root",hv_,sn_/*,mt_*/);
	//c3n->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_npaired_per_strip.pdf",hv_,sn_/*,mt_*/);
	//c3n->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_npaired_per_strip.png",hv_,sn_/*,mt_*/);
	//c3n->SaveAs(s);

	TCanvas *c4n = new TCanvas("c4n","The Ntuple canvas",200,10,780,780);
	c4n->cd();
	c4n->SetGrid();
	nFiredHR_per_strip->Draw("COLZ");
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nFiredHR_per_strip.root",hv_,sn_/*,mt_*/);
	//c4n->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nFiredHR_per_strip.pdf",hv_,sn_/*,mt_*/);
	//c4n->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nFiredHR_per_strip.png",hv_,sn_/*,mt_*/);
	//c4n->SaveAs(s);

	TCanvas *c5n = new TCanvas("c5n","The Ntuple canvas",200,10,780,780);
	c5n->cd();
	c5n->SetGrid();
	nFiredLR_per_strip->Draw("COLZ");
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nFiredLR_per_strip.root",hv_,sn_/*,mt_*/);
	//c5n->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nFiredLR_per_strip.pdf",hv_,sn_/*,mt_*/);
	//c5n->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nFiredLR_per_strip.png",hv_,sn_/*,mt_*/);
	//c5n->SaveAs(s);


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
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_forTimeWindow.png",hv_,sn_/*,mt_*/);
	ct3->SaveAs(s);

        TCanvas *ct3m = new TCanvas("ct3m","The Ntuple canvas",200,10,780,780);
        TPad *padt3m = new TPad("padt3m","This is padt3",0.05,0.05,0.98,0.55,21);
        TPad *padt5m = new TPad("padt5m","This is padt5",0.05,0.55,0.98,0.98,21);
        padt3m->UseCurrentStyle(); padt5m->UseCurrentStyle();
        padt3m->Draw(); padt5m->Draw();
        padt3m->cd();
        hLRT_temp_med->Draw();
        padt5m->cd();
        hHRT_temp_med->Draw();
        ct3m->cd();
        ct3m->Update();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_forTimeWindow_med.png",hv_,sn_/*,mt_*/);
        ct3m->SaveAs(s);


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
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_1d.root",hv_,sn_/*,mt_*/);
	c10->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_1d.pdf",hv_,sn_/*,mt_*/);
	c10->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_1d.png",hv_,sn_/*,mt_*/);
	c10->SaveAs(s);

	TCanvas *c10m = new TCanvas("c10m","The Ntuple canvas",200,10,780,780);
	TPad *pad10m = new TPad("pad10m","This is pad10",0.05,0.05,0.98,0.55,21);
	TPad *pad11m = new TPad("pad11m","This is pad11",0.05,0.55,0.98,0.98,21);
	pad10m->UseCurrentStyle(); pad11m->UseCurrentStyle();
	pad10m->Draw(); pad11m->Draw();
	pad10m->cd();
	pad10m->SetGrid();
	hLRT__med->Draw("Hist");
	pad11m->cd();
	pad11m->SetGrid();
	hHRT__med->Draw("Hist");
	c10m->cd();
	c10m->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_1d_med.root",hv_,sn_/*,mt_*/);
	//c10m->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_1d_med.pdf",hv_,sn_/*,mt_*/);
	//c10m->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_DeltaTrig_side_1d_med.png",hv_,sn_/*,mt_*/);
	//c10m->SaveAs(s);

        gStyle->SetOptStat(1111);
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
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nStrip_side_1d.root",hv_,sn_/*,mt_*/);
	c11->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nStrip_side_1d.pdf",hv_,sn_/*,mt_*/);
	c11->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nStrip_side_1d.png",hv_,sn_/*,mt_*/);
	c11->SaveAs(s);
        gStyle->SetOptStat(0);

	TCanvas *c11m = new TCanvas("c11m","The Ntuple canvas",200,10,780,780);
	TPad *pad20m = new TPad("pad20m","This is pad10",0.05,0.05,0.98,0.55,21);
	TPad *pad21m = new TPad("pad21m","This is pad11",0.05,0.55,0.98,0.98,21);
	pad20m->UseCurrentStyle(); pad21m->UseCurrentStyle();
	pad20m->Draw(); pad21m->Draw();
	pad20m->cd();
	pad20m->SetGrid();
	hLRT_2_med->Draw("Hist");
	pad21m->cd();
	pad21m->SetGrid();
	hHRT_2_med->Draw("Hist");
	c11m->cd();
	c11m->Update();
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nStrip_side_1d_med.root",hv_,sn_/*,mt_*/);
	//c11m->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nStrip_side_1d_med.pdf",hv_,sn_/*,mt_*/);
	//c11m->SaveAs(s);
	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_nStrip_side_1d_med.png",hv_,sn_/*,mt_*/);
	//c11m->SaveAs(s);


        TCanvas *c5np = new TCanvas("c5np","The Ntuple canvas",200,10,780,780);
        c5np->cd();
        c5np->SetGrid();
        histo_diff_petiroc->Draw();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_time_diff_same_pr.root",hv_,sn_/*,mt_*/);
        c5np->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_time_diff_same_pr.pdf",hv_,sn_/*,mt_*/);
        c5np->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_time_diff_same_pr.png",hv_,sn_/*,mt_*/);
        c5np->SaveAs(s);


	gStyle->SetOptFit(1111);
        TCanvas *cstrip_diff_h = new TCanvas("cstrip_diff_h","strip_diff_h",200,10,780,780);
        cstrip_diff_h->cd();
        cstrip_diff_h->SetGrid();
	strip_diff_h->Fit("gaus","","",-5,5);
        strip_diff_h->Draw();
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_strip_diff.root",hv_,sn_/*,mt_*/);
        c5np->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_strip_diff.pdf",hv_,sn_/*,mt_*/);
        c5np->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_strip_diff.png",hv_,sn_/*,mt_*/);
        cstrip_diff_h->SaveAs(s);
	gStyle->SetOptFit(0);

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
        s.Form("plots_904/_HV_%d_SN_%d_cluster_number_size.root",hv_,sn_/*,mt_*/);
        cg->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_cluster_number_size.pdf",hv_,sn_/*,mt_*/);
        cg->SaveAs(s);
        s.Form("plots_904/_HV_%d_SN_%d_cluster_number_size.png",hv_,sn_/*,mt_*/);
        cg->SaveAs(s);

	if(true){
		TCanvas *cx = new TCanvas("cx","The Ntuple canvas",200,10,780,780);
        	cx->cd();
        	cx->SetGrid();
                hHR_check_xtalk->Divide(hHR_check_xtalk_temp);
        	hHR_check_xtalk->Draw("HIST");
        	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_check_xtalk_HR.root",hv_,sn_/*,mt_*/);
        	cx->SaveAs(s);
        	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_check_xtalk_HR.pdf",hv_,sn_/*,mt_*/);
        	cx->SaveAs(s);
        	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_check_xtalk_HR.png",hv_,sn_/*,mt_*/);
        	cx->SaveAs(s);
	}

	if(true){
		TCanvas *cx = new TCanvas("cx","The Ntuple canvas",200,10,780,780);
        	cx->cd();
        	cx->SetGrid();
                hLR_check_xtalk->Divide(hLR_check_xtalk_temp);
        	hLR_check_xtalk->Draw("HIST");
        	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_check_xtalk_LR.root",hv_,sn_/*,mt_*/);
        	cx->SaveAs(s);
        	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_check_xtalk_LR.pdf",hv_,sn_/*,mt_*/);
        	cx->SaveAs(s);
        	s.Form("plots_904/_HV_%d_SN_%d_MaxTrig_check_xtalk_LR.png",hv_,sn_/*,mt_*/);
        	cx->SaveAs(s);
	}



	gApplication->Terminate();

}

