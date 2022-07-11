#include<stdlib.h>
#include<stdio.h>


void eff_curve(int sn, int count, ...){

	std::vector<char *> wp;
	wp.push_back("OR");
	wp.push_back("AND");
	wp.push_back("HR");
	wp.push_back("LR");
	wp.push_back("Cluster");

	std::vector<double> HVs;
	std::vector<float> effnocut_on, effloose_on, effhr_on, efflr_on, effnocut_out, effloose_out, effhr_out, efflr_out, effcluster_on, effcluster_out;
	HVs.clear();
	effnocut_on.clear();
	effloose_on.clear();
	effhr_on.clear();
	efflr_on.clear();
	effnocut_out.clear();
	effloose_out.clear();
	effhr_out.clear();
	efflr_out.clear();
	effcluster_on.clear();
	effcluster_out.clear();

	va_list list;
	va_start(list, count);
	double HV_val=0.;
	for(int j=0; j<count; j++){
		HV_val = va_arg(list, double);
		HVs.push_back(HV_val);
	}
	va_end(list);

	for(int j=0; j<count; j++){
		string fileName("");
		fileName = Form("outputs/_HV_%d_SN_%d_MaxTrig_.txt", (j+1), sn);
		FILE *fp = fopen(fileName.c_str(),"r");
		Int_t ncols;
		float effnocuton, efflooseon, effmediumon, efftighton, effhron, efflron, effclusteron, effcluster2on;
		float effnocutout, efflooseout, effmediumout, efftightout, effhrout, efflrout, effclusterout, effcluster2out;
		float temp;
		ncols=0;
		ncols=fscanf(fp,"%f", &effnocuton);
		effnocut_on.push_back(effnocuton);
		ncols=1;
		ncols=fscanf(fp,"%f", &efflooseon);
		effloose_on.push_back(efflooseon);
		ncols=2;
		ncols=fscanf(fp,"%f", &effhron);
		effhr_on.push_back(effhron);
		ncols=3;
		ncols=fscanf(fp,"%f", &efflron);
		efflr_on.push_back(efflron);
		ncols=4;
		ncols=fscanf(fp,"%f", &effclusteron);
		effcluster_on.push_back(effclusteron);
		ncols=5;
		ncols=fscanf(fp, "%f", &temp);
		ncols=6;
		ncols=fscanf(fp, "%*[^\n]\n");
		ncols=7;
		ncols=fscanf(fp,"%f", &effnocutout);
		effnocut_out.push_back(effnocutout);
		ncols=8;
		ncols=fscanf(fp,"%f", &efflooseout);
		effloose_out.push_back(efflooseout);
		ncols=9;
		ncols=fscanf(fp,"%f", &effhrout);
		effhr_out.push_back(effhrout);
		ncols=10;
		ncols=fscanf(fp,"%f", &efflrout);
		efflr_out.push_back(efflrout);
		ncols=11;
		ncols=fscanf(fp,"%f", &effclusterout);
		effcluster_out.push_back(effclusterout);
		fclose(fp);
	}


	//loop no wp
	for(int i=0;i<wp.size();i++){
		std::vector<float> eff_on, eff_out;
		eff_on.clear();
		eff_out.clear();
		if(i==0){eff_on=effnocut_on; eff_out=effnocut_out;}
		if(i==1){eff_on=effloose_on; eff_out=effloose_out;}
		if(i==2){eff_on=effhr_on; eff_out=effhr_out;}
		if(i==3){eff_on=efflr_on; eff_out=efflr_out;}
		if(i==4){eff_on=effcluster_on; eff_out=effcluster_out;}

		TGraphErrors* efficiency = new TGraphErrors();
		TGraphErrors* efficiency_temp = new TGraphErrors();

		//loop HV
		float HVmin = 6000;
		float HVmax = 7700;
		for(int j=0; j<count; j++){


			float eff_new = (eff_on.at(j)-eff_out.at(j))/(1.-eff_out.at(j));
			efficiency->SetPoint(j,(HVs.at(j)*1000),eff_new);
			efficiency->SetPointError(j,10,sqrt(eff_new*(1-eff_new))/sqrt(1000));

			efficiency_temp->SetPoint(j,(HVs.at(j)*1000),eff_new);
			efficiency_temp->SetPointError(j,10,sqrt(eff_new*(1-eff_new))/sqrt(1000));

			if(j==0) HVmin = HVs.at(j)*1000;
			if(j==(count-1)) HVmax = HVs.at(j)*1000;

		}



		TF1* sigmoid = new TF1("sigmoid","(1-sqrt((1-[0])*(1-[0])))/(1+exp([1]*([2]-x)))",6000,7400);
		sigmoid->SetParName(0,"#epsilon_{max}");
		sigmoid->SetParName(1,"#lambda");
		sigmoid->SetParName(2,"HV_{50%}");
		sigmoid->SetParameter(0,0.98);
		sigmoid->SetParameter(1,0.01);
		sigmoid->SetParameter(2,7000);
		efficiency->Fit(sigmoid);
		efficiency->SetMarkerStyle(22);
		efficiency->SetMarkerSize(2);
		string sTitle = "; HVeff; Eff";
		TH1D* Plotter = new TH1D("Plotter", sTitle.c_str(), 1, HVmin-100,HVmax+100);
		Plotter->SetStats(0);
		Plotter->SetMinimum(0.);
		Plotter->SetMaximum(1.08);
		Plotter->Draw();
		efficiency->Draw("P");
		double p1 = sigmoid->GetParameter(0);
		p1 = 1-sqrt((1-p1)*(1-p1));
		double p2 = sigmoid->GetParameter(1);
		double p3 = sigmoid->GetParameter(2);
		double uLimit = HVmax, lLimit = HVmin;
		TLatex* ltx = new TLatex();
		ltx->SetTextSize(0.04);
		double knee = p3 - log(1/0.95-1)/p2;
		double WP = knee+120;
		TLine* lWP = new TLine(WP, 0., WP, 1);
		lWP->SetLineStyle(2);
		lWP->Draw("SAME");
		cout << "WP: " << WP << endl;
		double shift = 0.5;
		double add = (uLimit-lLimit)/11., up = 0.3; 
		ltx->DrawLatex(HVmin, 0.48, Form("Eff(WP) = %.2f", sigmoid->Eval(WP)));
		ltx->DrawLatex(HVmin, 0.41, Form("WP = %.0f V", WP));
		ltx->DrawLatex(HVmin, 0.34, Form("knee = %.0f V", knee));
		ltx->DrawLatex(HVmin, 0.27, Form("HV 50% = %.0f V", p3));
		TLine* plateau = new TLine(lLimit-50, p1, uLimit+50, p1);
		plateau->SetLineStyle(2);
		plateau->Draw();
		add = (uLimit-lLimit)/11.;
		if ((knee - lLimit) < (uLimit-lLimit)*(3/11.)) add = knee + add;
		else add = lLimit+add;
		ltx->DrawLatex(add, p1+0.02, Form("plateau = %.2f", p1));
		string outfile = "Efficiency";
		string outfile_png = Form("ScanId_%d/Efficiency_SN%d_%s.png", sn, sn, wp.at(i));
		string outfile_pdf = Form("ScanId_%d/Efficiency_SN%d_%s.pdf", sn, sn, wp.at(i));
		string outfile_root =Form("ScanId_%d/Efficiency_SN%d_%s.root", sn, sn, wp.at(i));
		gPad->SaveAs(outfile_png.c_str());
		gPad->SaveAs(outfile_pdf.c_str());
		TFile* effout = new TFile(outfile_root.c_str(), "RECREATE");
		effout->cd();
		efficiency->Write("Efficiency");
		effout->Close();  


	}

}

