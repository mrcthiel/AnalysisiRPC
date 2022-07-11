void plot_chamber(int hv, int sn){
	for(int i=0;i<4;i++){
		string fileName("");
		if(i==0) fileName = Form("ScanId_%d/HV%d/_HV_%d_SN_%d_muons_X_Y_cls.root",sn,hv,hv,sn);
		if(i==1) fileName = Form("ScanId_%d/HV%d/_HV_%d_SN_%d_muons_X_Y.root",sn,hv,hv,sn);
		if(i==2) fileName = Form("ScanId_%d/HV%d/_HV_%d_SN_%d_gammas_X_Y_cls.root",sn,hv,hv,sn);
		if(i==3) fileName = Form("ScanId_%d/HV%d/_HV_%d_SN_%d_gammas_X_Y.root",sn,hv,hv,sn);

		TFile _file0m(fileName.c_str(),"read");

		TCanvas* nm =(TCanvas*)_file0m.Get("c4");
		TPad* nnm =(TPad*)nm->GetPrimitive("pad4");
		string fileNameh("");
		if(i==0 || i==2) fileNameh = Form("h2_XY_cls");
		if(i==1 || i==3) fileNameh = Form("h2_XY");
		TH2F *mmm = (TH2F*)nnm->GetPrimitive(fileNameh.c_str());

		TCanvas *c2 = new TCanvas("c4n","The Ntuple canvas",0,0,1000,1500);
		mmm->SetMarkerStyle(20);
		mmm->SetMarkerSize(0.25);
		mmm->SetStats(0);
		//	mmm->SetMarkerColor(4);

		mmm->Draw("P");
//		For RE41
		Double_t x[7] = {0.,279.,340.,539.,503.,0.,0.};
		Double_t y[7] = {1220.,1220.,1168.,43.,0.,0.,1220.};
//		For RE31
//                Double_t x[7] = {0.,234.45,290.2,530.33,501.71,0.,0.};
//                Double_t y[7] = {1467.,1467.,1416.22,42.65,0.,0.,1467.};

		TPolyLine *pline = new TPolyLine(7,x,y);
		pline->SetFillColor(0);
		pline->SetLineColor(2);
		pline->SetLineWidth(4);
		//	pline->Draw("f");
		pline->Draw();

		string fileNameout("");
		if(i==0) fileNameout = Form("ScanId_%d/HV%d/muon_strips_chamber_cls.png",sn,hv);
                if(i==1) fileNameout = Form("ScanId_%d/HV%d/muon_strips_chamber.png",sn,hv);
                if(i==2) fileNameout = Form("ScanId_%d/HV%d/gammas_strips_chamber_cls.png",sn,hv);
                if(i==3) fileNameout = Form("ScanId_%d/HV%d/gammas_strips_chamber.png",sn,hv);
		c2->SaveAs(fileNameout.c_str());
	}
}
