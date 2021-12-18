#include <TH2.h>
#include <TChain.h>
#include "stripTimes.h"
#include "CalibAlig.h"


Int_t           fCurrent; 
double V_conector = 299.792*0.6;

Long64_t LoadTree(Long64_t entry, TTree *fChain)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		return kTRUE;//      Notify();
	}
	return centry;
}




std::vector<double>  MuonWindow(int sn_, TTree *fChain, /*UInt_t nframe, ULong64_t frame[65535],*/ double time_min, double time_max, uint64_t nentries/*, TBranch *b_nframe, TBranch *b_frame*/){
   TBranch *b_nframe;
   TBranch *b_frame;

   UInt_t          nframe;
   ULong64_t       frame[65535];   //[nframe]
   fChain->SetBranchAddress("nframe", &nframe, &b_nframe);
   fChain->SetBranchAddress("frame", frame, &b_frame);

	uint64_t t_bc0[3];
	Long64_t nbytes = 0, nb = 0;

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
        double V_conector = 299.792*0.6;// mm/ns  - check it


/*



	double HR_to_conctor_n[48] = {90.543564, 80.410461, 70.289078, 60.173553, 50.069733, 40.903408, 36.713425, 32.523438, 31.166557, 36.356541, 40.546543, 49.202293, 59.317810, 69.433334, 79.548859, 89.664398, 56.088192, 66.203720, 76.319244, 86.428909, 96.550293, 106.665825, 116.775482, 126.896873, 137.006546, 147.122055, 157.243454, 167.358978, 177.474487, 187.590027, 197.699677, 207.821075, 175.118668, 185.379196, 195.481796, 205.597305, 215.848694, 225.834229, 235.949753, 246.065292, 256.180817, 266.302185, 276.671875, 286.533234, 301.558472, 314.650421, 328.493317, 343.113739};
	double LR_to_conctor_n[48] = {2313.407715, 2303.935059, 2294.468262, 2285.000977, 2275.528320, 2266.067383, 2256.606445, 2247.133789, 2238.666504, 2228.200195, 2218.733398, 2209.266113, 2199.787598, 2190.332520, 2180.877441, 2171.387207, 2120.251953, 2110.785156, 2101.318359, 2091.851074, 2082.384277, 2072.917480, 2063.444824, 2053.983887, 2044.516846, 2035.049927, 2025.583008, 2016.116333, 2006.649414, 1997.182495, 1987.721558, 1978.243042, 1926.101929, 1916.629150, 1907.162231, 1897.701416, 1888.234497, 1878.773438, 1869.300659, 1855.762939, 1844.198853, 1832.120117, 1819.977539, 1807.770386, 1795.510132, 1783.159668, 1770.754272, 1757.014038};
	double Strip_length_n[48]  = {1468.002441, 1468.022339, 1468.062012, 1468.121460, 1468.200928, 1468.300049, 1468.419067, 1468.557861, 1467.717041, 1468.895020, 1469.093262, 1469.311279, 1469.549194, 1469.806763, 1470.084229, 1470.381348, 1470.698242, 1471.034790, 1471.391235, 1471.767212, 1472.162964, 1472.578247, 1473.013306, 1473.468018, 1473.942383, 1474.436279, 1474.949707, 1475.482788, 1476.035400, 1476.607544, 1477.199097, 1477.810181, 1477.132324, 1479.090576, 1479.760010, 1480.448730, 1480.165649, 1481.884155, 1482.630737, 1481.655518, 1477.538940, 1472.481323, 1466.406738, 1462.278442, 1452.549316, 1436.809814, 1420.182617, 1403.590942};

	double time_corr_fine_n[3][32] = {
		{12.5625,12.5625,13.5609,18.5609,13.5592,16.5592,12.5621,15.5622,14.5592,13.5592,15.5384,16.5592,21.5592,15.5592,16.5656,19.5617,16.5653,20.5582,20.5924,15.5701,20.5952,16.5721,17.567,15.5655,17.5594,15.5589,21.5639,19.564,12.5617,14.5618,15.598,12.5632},
		{12.5592,11.5598,13.5592,18.5592,13.5053,16.5592,12.5592,15.5592,14.5641,13.5592,14.571,16.5641,21.571,14.5712,16.5645,19.5647,15.5625,20.5625,20.5625,15.5625,20.5625,16.5592,17.5592,15.5634,17.5592,14.6062,21.5639,19.564,12.5686,14.5477,15.5673,12.5668},
		{12.5633,12.5648,13.5697,18.5647,12.752,16.5642,12.5691,15.5692,14.5691,13.5642,15.563,16.5642,21.5592,15.5591,16.5592,19.5592,15.8607,20.5592,20.5592,15.5592,20.439,16.5578,17.5592,15.5592,17.3047,15.5592,21.5592,19.5592,12.5565,14.5575,15.5573,12.5574}

	};
*/

	TH1F* hLRT_temp = new TH1F("hLRT_muon_window","T (LR - Trig)" ,10000,-9000,1000);
	TH1F* hHRT_temp = new TH1F("hHRT_muon_window","T (HR - Trig)" ,10000,-9000,1000);

	vector<double> _vec;
	_vec.clear();

        int bad_trig=0;
        int trig_exists=0;
        int ntrig_event = 0;



	for (Long64_t jentry=0; jentry<nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry, fChain);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

	        stripTimes allStrips_temp[48];

		bad_trig=0;
		trig_exists=0;
		ntrig_event = 0;
		int trig[3];
		std::fill_n(trig, 3, 0);


		for (uint32_t i=0;i<nframe;i++)
		{
			if (m_channel(frame[i])==32) {
				t_bc0[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
			}
			else if (m_channel(frame[i])==33){
				if (trig[m_fpga(frame[i])] != 0) bad_trig = 1; 
				trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
				if ((trig[m_fpga(frame[i])]) < time_max and trig[m_fpga(frame[i])] > time_min) {
					trig_exists=1; 
					ntrig_event++; 
				}
			}
			else {
				//			any_strip_fired = 1; 
				if (c_side(m_channel(frame[i])) > 0.5){
					allStrips_temp[m_strip(frame[i])].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
				}
				if (c_side(m_channel(frame[i])) < 0.5){
					allStrips_temp[m_strip(frame[i])].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
				}
			}
		} 
		if (!(trig_exists && ntrig_event ==3) ) continue;
		if (bad_trig) continue;

		for (uint32_t i=0;i<48;i++) {
			for (uint32_t fhr=0; fhr<allStrips_temp[i].HRframe.size(); fhr++){
				double time_in_strip = allStrips_temp[i].HRtime[fhr] - trig[m_fpga(frame[allStrips_temp[i].HRframe[fhr]])] - HR_to_conctor[i]/V_conector;
				double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
				double time_to_conector_new = HR_to_conctor[i]/V_conector-HR_to_conctor[reference_strip]/V_conector;
				double HR_trig = time_in_strip_new + time_to_conector_new;
				hHRT_temp->Fill(HR_trig);
			}
			for (uint32_t flr=0; flr<allStrips_temp[i].LRframe.size(); flr++) { 
				double time_in_strip = allStrips_temp[i].LRtime[flr] - trig[m_fpga(frame[allStrips_temp[i].LRframe[flr]])] - LR_to_conctor[i]/V_conector;
				double time_in_strip_new = time_in_strip*Strip_length[reference_strip]/Strip_length[i];
				double time_to_conector_new = LR_to_conctor[i]/V_conector-LR_to_conctor[reference_strip]/V_conector;
				double LR_trig = time_in_strip_new + time_to_conector_new;
				hLRT_temp->Fill(LR_trig);
			}
		}
	}




        double muW_HR =0.;
        double muW_LR =0.;

        if(!(hHRT_temp->Integral()==0 || hLRT_temp->Integral()==0)){
		hHRT_temp->GetXaxis()->SetRangeUser(hHRT_temp->GetXaxis()->GetBinCenter(hHRT_temp->GetMaximumBin())-50,hHRT_temp->GetXaxis()->GetBinCenter(hHRT_temp->GetMaximumBin())+50);
		hLRT_temp->GetXaxis()->SetRangeUser(hLRT_temp->GetXaxis()->GetBinCenter(hLRT_temp->GetMaximumBin())-50,hLRT_temp->GetXaxis()->GetBinCenter(hLRT_temp->GetMaximumBin())+50);
		hHRT_temp->Fit("gaus");
		hLRT_temp->Fit("gaus");	
		TF1* gauss_plus_contan = new TF1("gauss_plus_contant","[0]*exp(-0.5*pow(((x-[1])/[2]),2))+x*[4]+[5]");
		gauss_plus_contan->SetParameters(hHRT_temp->GetFunction("gaus")->GetParameter(0),hHRT_temp->GetFunction("gaus")->GetParameter(1),hHRT_temp->GetFunction("gaus")->GetParameter(2),0,0);
		hLRT_temp->Fit(gauss_plus_contan);
	        muW_LR = gauss_plus_contan->GetParameter(1);

		hHRT_temp->Fit(gauss_plus_contan);
                muW_HR = gauss_plus_contan->GetParameter(1);
	}

        _vec.push_back(muW_LR);
        _vec.push_back(muW_HR);



	return _vec;

}
