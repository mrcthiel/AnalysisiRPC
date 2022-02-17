#include<stdlib.h>
#include<stdio.h>


void eff_clust_curve(int sn, int count, ...){

	std::vector<char *> wp;
	wp.push_back("Or");
	wp.push_back("And");
	wp.push_back("medium");
	wp.push_back("tight");
	wp.push_back("hr");
	wp.push_back("lr");
        wp.push_back("cluster");
        wp.push_back("cluster2");



	std::vector<double> HVs;
	std::vector<float> effnocut_on, effloose_on, effmedium_on, efftight_on, effhr_on, efflr_on, effnocut_out, effloose_out, effmedium_out, efftight_out, effhr_out, efflr_out, effcluster_on, effcluster_out, effcluster2_on, effcluster2_out;
	HVs.clear();
	effnocut_on.clear();
	effloose_on.clear();
	effmedium_on.clear();
	efftight_on.clear();
	effhr_on.clear();
	efflr_on.clear();
	effnocut_out.clear();
	effloose_out.clear();
	effmedium_out.clear();
	efftight_out.clear();
	effhr_out.clear();
	efflr_out.clear();
        effcluster_on.clear();
        effcluster_out.clear();
        effcluster2_on.clear();
        effcluster2_out.clear();

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
		ncols=fscanf(fp,"%f", &effmediumon);
		effmedium_on.push_back(effmediumon);
		ncols=3;
		ncols=fscanf(fp,"%f", &efftighton);
		efftight_on.push_back(efftighton);
		ncols=4;
		ncols=fscanf(fp,"%f", &effhron);
		effhr_on.push_back(effhron);
		ncols=5;
		ncols=fscanf(fp,"%f", &efflron);
		efflr_on.push_back(efflron);
                ncols=6;
                ncols=fscanf(fp,"%f", &effclusteron);
                effcluster_on.push_back(effclusteron);
                ncols=7;
                ncols=fscanf(fp,"%f", &effcluster2on);
                effcluster2_on.push_back(effcluster2on);
		ncols=8;
		ncols=fscanf(fp, "%f", &temp);
                ncols=9;
                ncols=fscanf(fp, "%f", &temp);
                ncols=10;
                ncols=fscanf(fp, "%*[^\n]\n");
                ncols=11;
		ncols=fscanf(fp,"%f", &effnocutout);
		effnocut_out.push_back(effnocutout);
		ncols=12;
		ncols=fscanf(fp,"%f", &efflooseout);
		effloose_out.push_back(efflooseout);
		ncols=13;
		ncols=fscanf(fp,"%f", &effmediumout);
		effmedium_out.push_back(effmediumout);
		ncols=14;
		ncols=fscanf(fp,"%f", &efftightout);
		efftight_out.push_back(efftightout);
		ncols=15;
		ncols=fscanf(fp,"%f", &effhrout);
		effhr_out.push_back(effhrout);
		ncols=16;
		ncols=fscanf(fp,"%f", &efflrout);
		efflr_out.push_back(efflrout);
                ncols=17;
                ncols=fscanf(fp,"%f", &effclusterout);
                effcluster_out.push_back(effclusterout);
                ncols=18;
                ncols=fscanf(fp,"%f", &effcluster2out);
                effcluster2_out.push_back(effcluster2out);
		//fp->fclose();
		fclose(fp);
	}


	vector<double> clust_size_v;
        vector<double> clust_mult_v;
	clust_size_v.clear();
	clust_mult_v.clear();



        for(int j=0; j<count; j++){
                string fileName("");
                fileName = Form("ScanId_%d/HV%d/plots_904/_HV_%d_SN_%d_MaxTrig_Cluster_size.root",sn,(j+1),(j+1),sn);
                TFile _file0(fileName.c_str(),"read");
                TCanvas* n =(TCanvas*)_file0.Get("cPairsC");
                TPad* nn =(TPad*)n->GetPrimitive("padPC");
                TH1F *mm = (TH1F*)nn->GetPrimitive("hClusterSize");
                cout << mm->GetMean() << endl;
                clust_size_v.push_back(mm->GetMean());
        }

        for(int j=0; j<count; j++){
		string fileName("");
		fileName = Form("ScanId_%d/HV%d/plots_904/_HV_%d_SN_%d_MaxTrig_N_Clusters.root",sn,(j+1),(j+1),sn);
		TFile _file0(fileName.c_str(),"read");
		TCanvas* n =(TCanvas*)_file0.Get("cPairsCC");
		TPad* nn =(TPad*)n->GetPrimitive("padPCC");
		TH1F *mm = (TH1F*)nn->GetPrimitive("hNClusters");
		cout << mm->GetMean() << endl;
		clust_mult_v.push_back(mm->GetMean());
	}

	//loop no wp
	for(int i=0;i<wp.size();i++){
		if(i==3) continue;
		std::vector<float> eff_on, eff_out;
		eff_on.clear();
		eff_out.clear();
		if(i==0){eff_on=effnocut_on; eff_out=effnocut_out;}
		if(i==1){eff_on=effloose_on; eff_out=effloose_out;}
		if(i==2){eff_on=effmedium_on; eff_out=effmedium_out;}
		if(i==3){eff_on=efftight_on; eff_out=efftight_out;}
		if(i==4){eff_on=effhr_on; eff_out=effhr_out;}
		if(i==5){eff_on=efflr_on; eff_out=efflr_out;}
                if(i==6){eff_on=effcluster_on; eff_out=effcluster_out;}
                if(i==7){eff_on=effcluster2_on; eff_out=effcluster2_out;}

		TGraphErrors* efficiency = new TGraphErrors();
                TGraphErrors* efficiency_temp = new TGraphErrors();

                TGraphErrors* g_clust_mult = new TGraphErrors();
                TGraphErrors* g_clust_size = new TGraphErrors();


		//loop HV
		float HVmin = 6000;
                float HVmax = 7400;
		for(int j=0; j<count; j++){

std::cout << "eff_on.at(" << j << "): " << eff_on.at(j) << "; " << "eff_out.at(" << j << "): " << eff_out.at(j) << "; HV: " << HVs.at(j)*1000 << std::endl;


			float eff_new = (eff_on.at(j)-eff_out.at(j))/(1.-eff_out.at(j));
			efficiency->SetPoint(j,(HVs.at(j)*1000),eff_new);
			efficiency->SetPointError(j,10,sqrt(eff_new*(1-eff_new))/sqrt(1000));

                        efficiency_temp->SetPoint(j,(HVs.at(j)*1000),eff_new);
                        efficiency_temp->SetPointError(j,10,sqrt(eff_new*(1-eff_new))/sqrt(1000));

                        g_clust_mult->SetPoint(j,(HVs.at(j)*1000),clust_mult_v.at(j));
                        g_clust_size->SetPoint(j,(HVs.at(j)*1000),clust_size_v.at(j));


			if(j==0) HVmin = HVs.at(j)*1000;
                        if(j==(count-1)) HVmax = HVs.at(j)*1000;

		}

                /*efficiency->GetXaxis()->SetLimits(6000,7500);
                efficiency_temp->GetXaxis()->SetLimits(6000,7500);
                g_clust_mult->GetXaxis()->SetLimits(6000,7500);
                g_clust_size->GetXaxis()->SetLimits(6000,7500);
		*/

        TCanvas *cg = new TCanvas("cg", "comp", 200, 10, 700, 500);
        TPad *p1g = new TPad("p1g", "", 0, 0, 1, 1);
        TPad *p2g = new TPad("p2g", "", 0, 0, 1, 1);
        p2g->SetFillStyle(4000);
        p1g->Draw();
        p1g->cd();



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
		TH1D* Plotter = new TH1D("Plotter", sTitle.c_str(), 1, HVmin-400,HVmax);
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
		double shift = 0.5;
		double add = (uLimit-lLimit)/11., up = 0.3; 
		ltx->DrawLatex(HVmin-300, 0.88, Form("Eff(WP) = %.2f", sigmoid->Eval(WP)));
		ltx->DrawLatex(HVmin-300, 0.81, Form("WP = %.0f V", WP));
		ltx->DrawLatex(HVmin-300, 0.74, Form("knee = %.0f V", knee));
		ltx->DrawLatex(HVmin-300, 0.67, Form("HV 50% = %.0f V", p3));
		TLine* plateau = new TLine(lLimit-50, p1, uLimit+50, p1);
		plateau->SetLineStyle(2);
		plateau->Draw();
		add = (uLimit-lLimit)/11.;
		if ((knee - lLimit) < (uLimit-lLimit)*(3/11.)) add = knee + add;
		else add = lLimit+add;
		ltx->DrawLatex(add, p1+0.02, Form("plateau = %.2f", p1));

   TLegend* legend = new TLegend(0.14,0.4,0.45,0.55);
   legend->SetTextFont(42);
   legend->SetBorderSize(0);  // no border
   legend->SetFillStyle(4000);
   legend->SetFillColor(0);   // Legend background should be white
   legend->SetTextSize(0.04); // Increase entry font size! 
   legend->AddEntry(efficiency,"Efficiency","P");
   legend->AddEntry(g_clust_mult,"Cluster Multiplicity","P");
   legend->AddEntry(g_clust_size,"Cluster Size","P");
   legend->Draw();


        gPad->Update();

        Double_t xmin = p1g->GetUxmin();
        Double_t xmax = p1g->GetUxmax();
        Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
        Double_t ymin = 0.;// g_clust_size->GetHistogram()->GetMinimum();
        Double_t ymax = 10.;//g_clust_size->GetHistogram()->GetMaximum();
        Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
        p2g->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
    

        p2g->Draw();
        p2g->cd();

        g_clust_mult->SetMarkerStyle(29);
        g_clust_mult->SetMarkerSize(2);
        g_clust_mult->SetMarkerColor(4);

        g_clust_size->SetMarkerStyle(34);
        g_clust_size->SetMarkerSize(2);
        g_clust_size->SetMarkerColor(8);

                
	g_clust_size->Draw("P");
        g_clust_mult->Draw("P SAME");

        gPad->Update();
        TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
        axis->SetTitle("Cluster Multiplicity and Size");
	//axis->SetLineColor(kRed);
        //axis->SetLabelColor(kRed);
        axis->Draw();
        gPad->Update();
        cg->cd();
        cg->Update();


		string outfile = "Efficiency";

		string outfile_png = Form("ScanId_%d/Efficiency_Clust_SN%d_%s.png", sn, sn, wp.at(i));
		string outfile_pdf = Form("ScanId_%d/Efficiency_Clust_SN%d_%s.pdf", sn, sn, wp.at(i));
		string outfile_root =Form("ScanId_%d/Efficiency_Clust_SN%d_%s.root", sn, sn, wp.at(i));

		cg->SaveAs(outfile_png.c_str());

		gPad->SaveAs(outfile_pdf.c_str());
		TFile* effout = new TFile(outfile_root.c_str(), "RECREATE");
		effout->cd();
		efficiency->Write("Efficiency");
		g_clust_mult->Write("g_clust_mult");
		g_clust_size->Write("g_clust_size");
		effout->Close();  


	}

}

