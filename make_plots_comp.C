
void make_plots_comp(char* label, int numb_scans, ...){


	int allscans[40] =  {606,608,603,604,609,605,607,601,592,593,594,597,598,599,227,231,228,229,230,200,199,195,196,201,219,213,214,217,218,619,620,623,621,622,617,618,612,613,611,614};

	char* alllabels[40] = {"RE31_186_Feb1_Source_0ff_606_50Dac","RE31_186_Feb1_Source_0ff_608_40Dac_HR","RE31_186_Feb1_Source_10_603_40Dac","RE31_186_Feb1_Source_10_604_40Dac_LR","RE31_186_Feb1_Source_10_609_40Dac_HR","RE31_186_Feb1_Source_off_605_40Dac","RE31_186_Feb1_Source_off_607_40Dac_LR","RE31_186_Feb2_Source_10_601_40Dac","RE31_186_Feb2_Source_off_592_40Dac","RE31_186_Feb2_Source_off_593_30Dac","RE31_186_Feb2_Source_off_594_50Dac","RE31_186_Feb2_Source_off_597_20Dac","RE31_186_Feb2_Source_off_598_40Dac_LR","RE31_186_Feb2_Source_off_599_40Dac_HR","RE31_186_Feb6_Source_46_227_40Dac","RE31_186_Feb6_Source_4.6_231_40Dac","RE31_186_Feb6_Source_off_228_40Dac","RE31_186_Feb6_Source_off_229_30Dac","RE31_186_Feb6_Source_off_230_20Dac","RE41_188_Feb1_Source_0ff_200_30Dac","RE41_188_Feb1_Source_10_199_30Dac","RE41_188_Feb1_Source_4.6_195_40Dac","RE41_188_Feb1_Source_4.6_196_30Dac","RE41_188_Feb1_Source_off_201_40Dac","RE41_188_Feb2_Source_10_219_40Dac","RE41_188_Feb2_Source_4.6_213_40Dac","RE41_188_Feb2_Source_4.6_214_30Dac","RE41_188_Feb2_Source_off_217_40Dac","RE41_188_Feb2_Source_off_218_30Dac","RE41_189_Feb7_Source_100_619_30Dac","RE41_189_Feb7_Source_10_620_40Dac","RE41_189_Feb7_Source_10_623_30Dac","RE41_189_Feb7_Source_4.6_621_40Dac","RE41_189_Feb7_Source_6.9_622_40Dac","RE41_189_Feb7_Source_off_617_40Dac","RE41_189_Feb7_Source_off_618_30Dac","RE41_189_Feb9_Source_10_612_40Dac","RE41_189_Feb9_Source_4.6_613_30Dac","RE41_189_Feb9_Source_off_611_40Dac","RE41_189_Feb9_Source_off_614_30Dac"};

	double leg_x_v[6] = {6220,6220,6220,7280,7280,7280};
	double leg_y_v[6] = {0.6,0.39,0.18,0.15,0.36,0.57};
	double leb_x_v[6] = {6220,6220,6220,6220,6220,6220};
	double leb_y_v[6] = {1.15,1.09,1.03,0.97,0.91,0.85};


        va_list list;
        va_start(list, numb_scans);
        int val=0;
	int sn[numb_scans];
        for(int j=0; j<numb_scans; j++){
                val = va_arg(list, int);
                sn[j] = val;
        }
        va_end(list);


        std::vector<char *> wp;
        wp.push_back("nocut");
        wp.push_back("loose");
        wp.push_back("medium");
        wp.push_back("cluster2");
        wp.push_back("hr");
        wp.push_back("lr");
        for(int ii=0;ii<wp.size();ii++){

		TCanvas *c = new TCanvas("c", "c",1200,900);

		vector<TGraphErrors *> scurves_h;
		scurves_h.clear();

		vector<double> p1_2_vec;
                vector<double> eff_vec;
                vector<double> WP_2_vec;
		p1_2_vec.clear();
                eff_vec.clear();
                WP_2_vec.clear();

        	for(int jj=0;jj<numb_scans;jj++){
	
			double HVmax = 6400;
			double HVmin = 7800;
			double uLimit = HVmax, lLimit = HVmin;

		        string fileName("");
		        fileName = Form("ScanId_%i/Efficiency_SN%i_%s.root",sn[jj],sn[jj],wp.at(ii));
			TFile file_2(fileName.c_str(),"read");
			TGraphErrors *efficiency_2 = (TGraphErrors*)file_2.Get("Efficiency");
			TF1* sigmoid_2 = new TF1("sigmoid2","(1-sqrt((1-[0])*(1-[0])))/(1+exp([1]*([2]-x)))",6400,7200);
			if((1+jj)==5) {sigmoid_2->SetLineColor(46);} else sigmoid_2->SetLineColor(1+jj);
			Int_t n_2;
			n_2=efficiency_2->GetN(); //get ploted array dimention
			TGraphErrors* efficiency2_2 = new TGraphErrors();
			string title("");
			title = Form(" %s - %s ; HVeff; Eff",label,wp.at(ii));
		        efficiency2_2->SetTitle(title.c_str());
			for(Int_t i=0; i<n_2; i++) {
				efficiency2_2->SetPoint(i,efficiency_2->GetPointX(i),efficiency_2->GetPointY(i));
				efficiency2_2->SetPointError(i,efficiency_2->GetErrorX(i),efficiency_2->GetErrorY(i));
			}
			sigmoid_2->SetParName(0,"#epsilon_{max}");
			sigmoid_2->SetParName(1,"#lambda");
			sigmoid_2->SetParName(2,"HV_{50%}");
			sigmoid_2->SetParameter(0,0.98);
			sigmoid_2->SetParameter(1,0.01);
			sigmoid_2->SetParameter(2,7000);
		        efficiency2_2->Fit(sigmoid_2);
			double p1_2 = sigmoid_2->GetParameter(0);
			p1_2 = 1-sqrt((1-p1_2)*(1-p1_2));
			double p2_2 = sigmoid_2->GetParameter(1);
			double p3_2 = sigmoid_2->GetParameter(2);
			efficiency2_2->SetMarkerStyle(22);
			efficiency2_2->SetMarkerSize(2);
			//efficiency2_2->SetMarkerColor(1+jj);
                        if((1+jj)==5) {efficiency2_2->SetMarkerColor(46);} else efficiency2_2->SetMarkerColor(1+jj);

			efficiency2_2->GetXaxis()->SetLimits(6200,7600);  //GetXaxis()->SetRangeUser(2000, 7500);
                        efficiency2_2->GetYaxis()->SetRangeUser(0, 1.2);

			scurves_h.push_back(efficiency2_2);

			double knee_2 = p3_2 - log(1/0.95-1)/p2_2;
			double WP_2 = knee_2+120;

			p1_2_vec.push_back(p1_2);
			eff_vec.push_back(sigmoid_2->Eval(WP_2));
			WP_2_vec.push_back(WP_2);

		}



	scurves_h.at(0)->Draw("AP");
        TLatex* ltx_2= new TLatex();
        ltx_2->SetTextSize(0.04);
        ltx_2->SetTextColor(1);
        ltx_2->DrawLatex(leg_x_v[0], leg_y_v[0], Form("plateau = %.2f", p1_2_vec.at(0)));
        ltx_2->DrawLatex(leg_x_v[0], leg_y_v[0]-0.07, Form("Eff(WP) = %.2f", eff_vec.at(0)));
        ltx_2->DrawLatex(leg_x_v[0], leg_y_v[0]-2*0.07, Form("WP = %.0f V", WP_2_vec.at(0)));
        TLatex* ltxLeg2= new TLatex();
        ltxLeg2->SetTextSize(0.03);
        ltxLeg2->SetTextColor(1);
	for(int isn=0;isn<40;isn++){
		if(allscans[isn]==sn[0]) ltxLeg2->DrawLatex(leb_x_v[0], leb_y_v[0], alllabels[isn]);
	}



	for(int iii = 1; iii<scurves_h.size();iii++){
		scurves_h.at(iii)->Draw("P SAME");

	        TLatex* ltx_3= new TLatex();
	        ltx_3->SetTextSize(0.04);
	        //ltx_3->SetTextColor(iii);
                if((iii+1)==5) {ltx_3->SetTextColor(46);} else ltx_3->SetTextColor(iii+1);

	        ltx_3->DrawLatex(leg_x_v[iii], leg_y_v[iii], Form("plateau = %.2f", p1_2_vec.at(iii)));
	        ltx_3->DrawLatex(leg_x_v[iii], leg_y_v[iii]-0.07, Form("Eff(WP) = %.2f", eff_vec.at(iii)));
	        ltx_3->DrawLatex(leg_x_v[iii], leg_y_v[iii]-2*0.07, Form("WP = %.0f V", WP_2_vec.at(iii)));
	        TLatex* ltxLeg3= new TLatex();
	        ltxLeg3->SetTextSize(0.03);
	        //ltxLeg3->SetTextColor(iii);
                if((iii+1)==5) {ltxLeg3->SetTextColor(46);} else ltxLeg3->SetTextColor(iii+1);
	        for(int isn=0;isn<40;isn++){
	                if(allscans[isn]==sn[iii]) ltxLeg3->DrawLatex(leb_x_v[iii], leb_y_v[iii], alllabels[isn]);
	        }
	}		



        string fileName9("");
        fileName9 = Form("%s_%s.png",label,wp.at(ii));

        c->SaveAs(fileName9.c_str());

	}
}
