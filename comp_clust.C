
void comp_clust(){

        std::vector<char *> wp;
        wp.push_back("And");
        wp.push_back("Or");
        wp.push_back("medium");
        wp.push_back("cluster");
        wp.push_back("hr");
        wp.push_back("lr");
        for(int ii=0;ii<wp.size();ii++){




	TCanvas *c = new TCanvas("c", "c",200, 10, 700, 500);
        TPad *p1g = new TPad("p1g", "", 0, 0, 1, 1);
        TPad *p2g = new TPad("p2g", "", 0, 0, 1, 1);
        p2g->SetFillStyle(4000);
        p1g->SetGrid();
        p1g->Draw();
        p1g->cd();
        


	double HVmax = 6400;
	double HVmin = 7800;
	double uLimit = HVmax, lLimit = HVmin;
        string fileName("");
        fileName = Form("ScanId_756/Efficiency_Clust_SN756_%s.root",wp.at(ii));
	TFile file_2(fileName.c_str(),"read");
	TGraphErrors *efficiency_2 = (TGraphErrors*)file_2.Get("Efficiency");
	TGraphErrors *efficiency_temp = (TGraphErrors*)file_2.Get("Efficiency");
	TF1* sigmoid_2 = new TF1("sigmoid2","(1-sqrt((1-[0])*(1-[0])))/(1+exp([1]*([2]-x)))",6400,7200);
	sigmoid_2->SetLineColor(4);
	Int_t n_2;
	n_2=efficiency_2->GetN(); //get ploted array dimention
	TGraphErrors* efficiency2_2 = new TGraphErrors();
	string title("");
	title = Form(" RE4/1 190 , FEB 18 , Petiroc 2C, FEBv2_2, 40fC - %s ; HVeff; Eff",wp.at(ii));
        efficiency2_2->SetTitle(title.c_str());
//        efficiency2_2->SetTitle(" SourceOFF, 40Dac ; HVeff; Eff");
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
	efficiency_temp->SetMarkerStyle(22);
	efficiency2_2->SetMarkerStyle(22);
	efficiency2_2->SetMarkerSize(2);
	efficiency2_2->SetMarkerColor(4);
    	efficiency2_2->GetYaxis()->SetRangeUser(0., 1.1);
	TAxis *axis1 = efficiency2_2->GetXaxis();
	axis1->SetLimits(5900,7600); 
	efficiency2_2->Draw("AP");
	TLatex* ltx_2= new TLatex();
	ltx_2->SetTextSize(0.035);
	ltx_2->SetTextColor(4);
	double knee_2 = p3_2 - log(1/0.95-1)/p2_2;
	double WP_2 = knee_2+120;
	ltx_2->DrawLatex(5950, 0.36, Form("plateau = %.2f", p1_2));
	ltx_2->DrawLatex(5950, 0.29, Form("Eff(WP) = %.2f", sigmoid_2->Eval(WP_2)));
	ltx_2->DrawLatex(5950, 0.22, Form("WP = %.0f V", WP_2));
//	ltx_2->DrawLatex(7650, 0.14, Form("knee = %.0f V", knee_2));
//	ltx_2->DrawLatex(7650, 0.07, Form("HV 50% = %.0f V", p3_2));




        string fileName2("");
        fileName2 = Form("ScanId_757/Efficiency_Clust_SN757_%s.root",wp.at(ii));
        TFile file_3(fileName2.c_str(),"read");
	TGraphErrors *efficiency_3 = (TGraphErrors*)file_3.Get("Efficiency");
	TF1* sigmoid_3 = new TF1("sigmoid3","(1-sqrt((1-[0])*(1-[0])))/(1+exp([1]*([2]-x)))",6400,7200);
	sigmoid_3->SetLineColor(8);
	Int_t n_3;
	n_3=efficiency_3->GetN(); //get ploted array dimention
	TGraphErrors* efficiency2_3 = new TGraphErrors();
	for(Int_t i=0; i<n_3; i++) {
		efficiency2_3->SetPoint(i,efficiency_3->GetPointX(i),efficiency_3->GetPointY(i));
		efficiency2_3->SetPointError(i,efficiency_3->GetErrorX(i),efficiency_3->GetErrorY(i));
	}
	sigmoid_3->SetParName(0,"#epsilon_{max}");
	sigmoid_3->SetParName(1,"#lambda");
	sigmoid_3->SetParName(2,"HV_{50%}");
	sigmoid_3->SetParameter(0,0.98);
	sigmoid_3->SetParameter(1,0.01);
	sigmoid_3->SetParameter(2,7000);
        efficiency2_3->Fit(sigmoid_3);
	double p1_3 = sigmoid_3->GetParameter(0);
	p1_3 = 1-sqrt((1-p1_3)*(1-p1_3));
	double p2_3 = sigmoid_3->GetParameter(1);
	double p3_3 = sigmoid_3->GetParameter(2);
	efficiency2_3->SetMarkerStyle(22);
	efficiency2_3->SetMarkerSize(2);
	efficiency2_3->SetMarkerColor(8);
	efficiency2_3->Draw("P SAME");
	TLatex* ltx_3= new TLatex();
	ltx_3->SetTextSize(0.035);
	ltx_3->SetTextColor(8);
	double knee_3 = p3_3 - log(1/0.95-1)/p2_3;
	double WP_3 = knee_3+120;
	ltx_3->DrawLatex(5950, 0.15, Form("plateau = %.2f", p1_3));
	ltx_3->DrawLatex(5950, 0.08, Form("Eff(WP) = %.2f", sigmoid_3->Eval(WP_3)));
	ltx_3->DrawLatex(5950, 0.01, Form("WP = %.0f V", WP_3));
//	ltx_3->DrawLatex(7300, 0.14, Form("knee = %.0f V", knee_3));
//	ltx_3->DrawLatex(7300, 0.07, Form("HV 50% = %.0f V", p3_3));




        string fileName3("");
        fileName3 = Form("ScanId_758/Efficiency_Clust_SN758_%s.root",wp.at(ii));
        TFile file_4(fileName3.c_str(),"read");
	TGraphErrors *efficiency_4 = (TGraphErrors*)file_4.Get("Efficiency");
	TF1* sigmoid_4 = new TF1("sigmoid4","(1-sqrt((1-[0])*(1-[0])))/(1+exp([1]*([2]-x)))",6400,7200);
	sigmoid_4->SetLineColor(6);
	Int_t n_4;
	n_4=efficiency_4->GetN(); //get ploted array dimention
	TGraphErrors* efficiency2_4 = new TGraphErrors();
	for(Int_t i=0; i<n_4; i++) {
		efficiency2_4->SetPoint(i,efficiency_4->GetPointX(i),efficiency_4->GetPointY(i));
		efficiency2_4->SetPointError(i,efficiency_4->GetErrorX(i),efficiency_4->GetErrorY(i));
	}
	sigmoid_4->SetParName(0,"#epsilon_{max}");
	sigmoid_4->SetParName(1,"#lambda");
	sigmoid_4->SetParName(2,"HV_{50%}");
	sigmoid_4->SetParameter(0,0.98);
	sigmoid_4->SetParameter(1,0.01);
	sigmoid_4->SetParameter(2,7000);
        efficiency2_4->Fit(sigmoid_4);
	double p1_4 = sigmoid_4->GetParameter(0);
	p1_4 = 1-sqrt((1-p1_4)*(1-p1_4));
	double p2_4 = sigmoid_4->GetParameter(1);
	double p3_4 = sigmoid_4->GetParameter(2);
	efficiency2_4->SetMarkerStyle(22);
	efficiency2_4->SetMarkerSize(2);
	efficiency2_4->SetMarkerColor(6);
	efficiency2_4->Draw("P SAME");
	TLatex* ltx_4= new TLatex();
	ltx_4->SetTextSize(0.035);
	ltx_4->SetTextColor(6);
	double knee_4 = p3_4 - log(1/0.95-1)/p2_4;
	double WP_4 = knee_4+120;
	ltx_4->DrawLatex(6300, 0.57, Form("plateau = %.2f", p1_4));
	ltx_4->DrawLatex(6300, 0.5, Form("Eff(WP) = %.2f", sigmoid_4->Eval(WP_4)));
	ltx_4->DrawLatex(6300, 0.43, Form("WP = %.0f V", WP_4));
//	ltx_4->DrawLatex(6850, 0.14, Form("knee = %.0f V", knee_4));
//	ltx_4->DrawLatex(6850, 0.07, Form("HV 50% = %.0f V", p3_4));



        string fileName4("");
        fileName4 = Form("ScanId_759/Efficiency_Clust_SN759_%s.root",wp.at(ii));
        TFile file_5(fileName4.c_str(),"read");
	TGraphErrors *efficiency_5 = (TGraphErrors*)file_5.Get("Efficiency");
	TF1* sigmoid_5 = new TF1("sigmoid5","(1-sqrt((1-[0])*(1-[0])))/(1+exp([1]*([2]-x)))",6400,7200);
	sigmoid_5->SetLineColor(2);
	Int_t n_5;
	n_5=efficiency_5->GetN(); //get ploted array dimention
	TGraphErrors* efficiency2_5 = new TGraphErrors();
	for(Int_t i=0; i<n_4; i++) {
		efficiency2_5->SetPoint(i,efficiency_5->GetPointX(i),efficiency_5->GetPointY(i));
		efficiency2_5->SetPointError(i,efficiency_5->GetErrorX(i),efficiency_5->GetErrorY(i));
	}
	sigmoid_5->SetParName(0,"#epsilon_{max}");
	sigmoid_5->SetParName(1,"#lambda");
	sigmoid_5->SetParName(2,"HV_{50%}");
	sigmoid_5->SetParameter(0,0.98);
	sigmoid_5->SetParameter(1,0.01);
	sigmoid_5->SetParameter(2,7000);
        efficiency2_5->Fit(sigmoid_5);
	double p1_5 = sigmoid_5->GetParameter(0);
	p1_5 = 1-sqrt((1-p1_5)*(1-p1_5));
	double p2_5 = sigmoid_5->GetParameter(1);
	double p3_5 = sigmoid_5->GetParameter(2);
	efficiency2_5->SetMarkerStyle(22);
	efficiency2_5->SetMarkerSize(2);
	efficiency2_5->SetMarkerColor(2);
	efficiency2_5->Draw("P SAME");
	TLatex* ltx_5= new TLatex();
	ltx_5->SetTextSize(0.035);
	ltx_5->SetTextColor(2);
	double knee_5 = p3_5 - log(1/0.95-1)/p2_5;
	double WP_5 = knee_5+120;
	ltx_5->DrawLatex(5950, 0.57, Form("plateau = %.2f", p1_5));
	ltx_5->DrawLatex(5950, 0.50, Form("Eff(WP) = %.2f", sigmoid_5->Eval(WP_5)));
	ltx_5->DrawLatex(5950, 0.43, Form("WP = %.0f V", WP_5));
//	ltx_5->DrawLatex(6400, 0.54, Form("knee = %.0f V", knee_5));
//	ltx_5->DrawLatex(6400, 0.47, Form("HV 50% = %.0f V", p3_5));

        TLatex* ltxLeg5= new TLatex();
        ltxLeg5->SetTextSize(0.035);
        ltxLeg5->SetTextColor(2);
        ltxLeg5->DrawLatex(5950, 1.05, "Source 4.6; 40ns deadtime; 1.4 kHz/cm2");

        TLatex* ltxLeg2= new TLatex();
        ltxLeg2->SetTextSize(0.035);
        ltxLeg2->SetTextColor(4);
        ltxLeg2->DrawLatex(5950, 0.99, "Source OFF; 40ns deadtime");

        TLatex* ltxLeg3= new TLatex();
        ltxLeg3->SetTextSize(0.035);
        ltxLeg3->SetTextColor(8);
        ltxLeg3->DrawLatex(5950, 0.93, "Source OFF; AUTORESET ON");

        TLatex* ltxLeg4= new TLatex();
        ltxLeg4->SetTextSize(0.035);
        ltxLeg4->SetTextColor(6);
        ltxLeg4->DrawLatex(5950, 0.87, "Source 4.6; AUTORESET ON; 1.4 kHz/cm2");

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


	TGraphErrors *efficiency_2_size = (TGraphErrors*)file_2.Get("g_clust_size");
	TGraphErrors *efficiency_2_mult = (TGraphErrors*)file_2.Get("g_clust_mult");
	TGraphErrors *efficiency_temp_size = (TGraphErrors*)file_2.Get("g_clust_size");
	TGraphErrors *efficiency_temp_mult = (TGraphErrors*)file_2.Get("g_clust_mult");
	efficiency_temp_size->SetMarkerColor(1);
	efficiency_temp_mult->SetMarkerColor(1);
	efficiency_2_size->SetMarkerColor(4);
	efficiency_2_mult->SetMarkerColor(4);
	efficiency_2_size->SetMarkerSize(1);
	efficiency_2_size->Draw("P");
	efficiency_2_mult->Draw("P SAME");
	
	TGraphErrors *efficiency_3_size = (TGraphErrors*)file_3.Get("g_clust_size");
	TGraphErrors *efficiency_3_mult = (TGraphErrors*)file_3.Get("g_clust_mult");
	efficiency_3_size->SetMarkerColor(8);
	efficiency_3_mult->SetMarkerColor(8);
	efficiency_3_size->SetMarkerSize(1);
	efficiency_3_size->Draw("P SAME");
	efficiency_3_mult->Draw("P SAME");
	
	TGraphErrors *efficiency_4_size = (TGraphErrors*)file_4.Get("g_clust_size");
	TGraphErrors *efficiency_4_mult = (TGraphErrors*)file_4.Get("g_clust_mult");
	efficiency_4_size->SetMarkerColor(6);
	efficiency_4_mult->SetMarkerColor(6);
	efficiency_4_size->SetMarkerSize(1);
	efficiency_4_size->Draw("P SAME");
	efficiency_4_mult->Draw("P SAME");

        TGraphErrors *efficiency_5_size = (TGraphErrors*)file_5.Get("g_clust_size");
        TGraphErrors *efficiency_5_mult = (TGraphErrors*)file_5.Get("g_clust_mult");
        efficiency_5_size->SetMarkerColor(2);
        efficiency_5_mult->SetMarkerColor(2);
        efficiency_5_size->SetMarkerSize(1);
        efficiency_5_size->Draw("P SAME");
        efficiency_5_mult->Draw("P SAME");


   TLegend* legend = new TLegend(0.12,0.55,0.45,0.7);
   legend->SetTextFont(42);
   legend->SetBorderSize(0);  // no border
   legend->SetFillStyle(4000);
   legend->SetFillColor(0);   // Legend background should be white
   legend->SetTextSize(0.04); // Increase entry font size! 
   legend->AddEntry(efficiency_temp,"Efficiency","P");
   legend->AddEntry(efficiency_temp_mult,"Cluster Multiplicity","P");
   legend->AddEntry(efficiency_temp_size,"Cluster Size","P");
   legend->Draw();


        gPad->Update();
        TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
        axis->SetTitle("Cluster Multiplicity and Size");
        axis->Draw();
        gPad->Update();
        c->cd();
        c->Update();


        string fileName9("");
        fileName9 = Form("Comp_clust_%s.png",wp.at(ii));

	c->SaveAs(fileName9.c_str());

	}
}
