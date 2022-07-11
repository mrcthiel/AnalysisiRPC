void plot_chamber_approval(){
		string fileName_1("");
		string fileName_2("");
		string fileName_3("");
		string fileName_4("");
		string fileName_5("");
		//fileName_2 = Form("_HV_11_SN_703_muons_X_Y_cls.root");
		//fileName_3 = Form("_HV_11_SN_704_muons_X_Y_cls.root");
		//fileName_1 = Form("_HV_11_SN_705_muons_X_Y_cls.root");
		fileName_2 = Form("ScanId_950/HV7/_HV_7_SN_950_muons_X_Y_cls.root");
		fileName_1 = Form("ScanId_958/HV7/_HV_7_SN_958_muons_X_Y_cls.root");
		fileName_3 = Form("ScanId_972/HV7/_HV_7_SN_972_muons_X_Y_cls.root");
		fileName_4 = Form("ScanId_981/HV7/_HV_7_SN_981_muons_X_Y_cls.root");
		fileName_5 = Form("ScanId_999/HV7/_HV_7_SN_999_muons_X_Y_cls.root");

		TFile _file1m(fileName_1.c_str(),"read");
		TFile _file2m(fileName_2.c_str(),"read");
		TFile _file3m(fileName_3.c_str(),"read");
		TFile _file4m(fileName_4.c_str(),"read");
		TFile _file5m(fileName_5.c_str(),"read");

		string fileNameh("");
		fileNameh = Form("h2_XY_cls");

		TCanvas* nm_1 =(TCanvas*)_file1m.Get("c4");
		TPad* nnm_1 =(TPad*)nm_1->GetPrimitive("pad4");
		TH2F *mmm_1 = (TH2F*)nnm_1->GetPrimitive(fileNameh.c_str());


		TCanvas* nm_2 =(TCanvas*)_file2m.Get("c4");
		TPad* nnm_2 =(TPad*)nm_2->GetPrimitive("pad4");
		TH2F *mmm_2 = (TH2F*)nnm_2->GetPrimitive(fileNameh.c_str());

		TCanvas* nm_3 =(TCanvas*)_file3m.Get("c4");
		TPad* nnm_3 =(TPad*)nm_3->GetPrimitive("pad4");
		TH2F *mmm_3 = (TH2F*)nnm_3->GetPrimitive(fileNameh.c_str());

		TCanvas* nm_4 =(TCanvas*)_file4m.Get("c4");
		TPad* nnm_4 =(TPad*)nm_4->GetPrimitive("pad4");
		TH2F *mmm_4 = (TH2F*)nnm_4->GetPrimitive(fileNameh.c_str());

		TCanvas* nm_5 =(TCanvas*)_file5m.Get("c4");
		TPad* nnm_5 =(TPad*)nm_5->GetPrimitive("pad4");
		TH2F *mmm_5 = (TH2F*)nnm_5->GetPrimitive(fileNameh.c_str());

		TCanvas *c2 = new TCanvas("c4n","The Ntuple canvas",0,0,600,600);
		mmm_1->SetMarkerStyle(20);
		mmm_1->SetMarkerSize(0.25);
		mmm_1->GetXaxis()->SetLabelSize(0.04);
		mmm_1->SetMarkerColor(kGreen+3);
		mmm_1->SetStats(0);
		mmm_1->Draw("P");

		mmm_2->SetMarkerStyle(20);
		mmm_2->SetMarkerSize(0.25);
		mmm_2->SetMarkerColor(4);
		mmm_2->GetXaxis()->SetLabelSize(0.04);
		mmm_2->SetStats(0);
		mmm_2->Draw("PSAME");

		mmm_3->SetMarkerStyle(20);
		mmm_3->SetMarkerSize(0.25);
		mmm_3->SetMarkerColor(6);
		mmm_3->GetXaxis()->SetLabelSize(0.04);
		mmm_3->SetStats(0);
		mmm_3->Draw("PSAME");

		mmm_4->SetMarkerStyle(20);
		mmm_4->SetMarkerSize(0.25);
		mmm_4->SetMarkerColor(1);
		mmm_4->GetXaxis()->SetLabelSize(0.04);
		mmm_4->SetStats(0);
		mmm_4->Draw("PSAME");

		mmm_5->SetMarkerStyle(20);
		mmm_5->SetMarkerSize(0.25);
		mmm_5->SetMarkerColor(2);
		mmm_5->GetXaxis()->SetLabelSize(0.04);
		mmm_5->SetStats(0);
		mmm_5->Draw("PSAME");
		//              For RE31
       //         Double_t x[7] = {0.,234.45,290.2,530.33,501.71,0.,0.};
       //         Double_t y[7] = {1467.,1467.,1416.22,42.65,0.,0.,1467.};
//              For RE41
              Double_t x[7] = {0.,279.,340.,539.,503.,0.,0.};
              Double_t y[7] = {1220.,1220.,1168.,43.,0.,0.,1220.};

		TPolyLine *pline = new TPolyLine(7,x,y);
		pline->SetFillColor(0);
		pline->SetLineColor(2);
		pline->SetLineWidth(4);
		//	pline->Draw("f");
		pline->Draw();

		string fileNameout("");
		fileNameout = Form("forICHEP_muon_strips_chamber_cls_FEB18_corrAlign.root");
		c2->SaveAs(fileNameout.c_str());
}
