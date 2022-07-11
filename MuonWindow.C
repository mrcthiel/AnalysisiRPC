#include <TH2.h>
#include <TChain.h>
#include "stripTimes.h"
#include "CalibAlig.h"


Int_t           fCurrent; 

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




std::vector<double>  MuonWindow(int sn_, TTree *fChain, double time_min, double time_max, uint64_t nentries, bool isFEBv2r2){


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
        for(int iref = 1; iref < 48; iref++){
                if(Strip_length[iref] < Strip_length[reference_strip]) reference_strip = iref;
        }
        double V_conector = 299.792*0.568;// mm/ns  - check it


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
				if ((!isFEBv2r2 && c_side(m_channel(frame[i]))>0.5) || (isFEBv2r2 && strip_side==1)){
					if(isFEBv2r2){
						allStrips_temp[strip_numb].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
					} else {
						allStrips_temp[m_strip(frame[i])].addLRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
					}
				}
				if ((!isFEBv2r2 && c_side(m_channel(frame[i]))<0.5) || (isFEBv2r2 && strip_side==0)){
					if(isFEBv2r2){
                                        	allStrips_temp[strip_numb].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
                                        } else {
						allStrips_temp[m_strip(frame[i])].addHRframe(i,m_time(m_traw(frame[i]))-time_corr_fine[m_fpga(frame[i])][m_channel(frame[i])]);
					}
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




        double muW_HR =-1200.;
        double muW_LR =-1200.;

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
