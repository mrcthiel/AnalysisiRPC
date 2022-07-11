void plot_and_save_2D(TH2F *h, TString sss, string ss, int hv_,int sn_, TString style){
	gStyle->SetOptStat(0);
	TCanvas *c2 = new TCanvas("c4","The Ntuple canvas",200,10,900,780);
	TPad *pad4 = new TPad("pad4","This is pad4",0.02,0.02,0.98,0.98,21);
	pad4->UseCurrentStyle();
	pad4->Draw();
	pad4->cd();
	pad4->SetGrid();
	pad4->SetRightMargin(0.12);
	h->GetYaxis()->SetLabelOffset(0.007);
	h->GetYaxis()->SetTitleOffset(1.5);
	h->SetLabelSize(0.02,"x");
	h->Draw(style);
	c2->cd();
	c2->Update();
	TString s("");
	s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.root",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
	c2->SaveAs(s);
	s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.pdf",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
	c2->SaveAs(s);
	s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.png",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
	c2->SaveAs(s);
}

void plot_and_save_1D(TH1F *h, TString sss, string ss, int hv_,int sn_, int gst=0, TString ssss=""){
	gStyle->SetOptStat(gst);
	TCanvas *c22 = new TCanvas("c42C","The Ntuple canvas",200,10,900,780);
	TPad *pad42C = new TPad("pad42C","This is pad42",0.02,0.02,0.98,0.98,21);
	pad42C->UseCurrentStyle();
	pad42C->Draw();
	pad42C->cd();
        pad42C->SetGrid();
	pad42C->SetRightMargin(0.12);
	h->GetYaxis()->SetTitleOffset(1.5);
	h->SetLabelSize(0.02,"x");
	h->Draw();
	c22->cd();
	c22->Update();
	TString s("");
	/*s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.root",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        if(plts_per_strip)*/ s.Form("ScanId_%d/HV%d/%s_HV_%d_SN_%d_%s%s.root",sn_,hv_,ssss.Data(),hv_,sn_,sss.Data(),ss.c_str());
	c22->SaveAs(s);
	/*s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.pdf",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        if(plts_per_strip)*/ s.Form("ScanId_%d/HV%d/%s_HV_%d_SN_%d_%s%s.pdf",sn_,hv_,ssss.Data(),hv_,sn_,sss.Data(),ss.c_str());
	c22->SaveAs(s);
	/*s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.png",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        if(plts_per_strip)*/ s.Form("ScanId_%d/HV%d/%s_HV_%d_SN_%d_%s%s.png",sn_,hv_,ssss.Data(),hv_,sn_,sss.Data(),ss.c_str());
	c22->SaveAs(s);
}

void plot_and_save_2x_1D(TH1F *h1, TH1F *h2, TString sss, string ss, int hv_,int sn_, int gst=0, bool log=false, TString ssss=""){
        gStyle->SetOptStat(gst);
        TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",200,10,780,780);
        TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.05,0.98,0.55,21);
        TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.55,0.98,0.98,21);
        pad1->UseCurrentStyle(); pad2->UseCurrentStyle();
        pad1->Draw(); pad2->Draw();
        pad1->cd();
        pad1->SetGrid();
        if(log) pad1->SetLogy();
        h1->Draw("Hist");
        pad2->cd();
        pad2->SetGrid();
        if(log) pad2->SetLogy();
        h2->Draw("Hist");
        c1->cd();
        c1->Update();
        TString s("");
        /*s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.root",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
	if(plts_per_strip)*/ s.Form("ScanId_%d/HV%d/%s_HV_%d_SN_%d_%s%s.root",sn_,hv_,ssss.Data(),hv_,sn_,sss.Data(),ss.c_str());
        c1->SaveAs(s);
        /*s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.pdf",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        if(plts_per_strip)*/ s.Form("ScanId_%d/HV%d/%s_HV_%d_SN_%d_%s%s.pdf",sn_,hv_,ssss.Data(),hv_,sn_,sss.Data(),ss.c_str());
        c1->SaveAs(s);
        /*s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.png",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        if(plts_per_strip)*/ s.Form("ScanId_%d/HV%d/%s_HV_%d_SN_%d_%s%s.png",sn_,hv_,ssss.Data(),hv_,sn_,sss.Data(),ss.c_str());
        c1->SaveAs(s);

}

void plot_and_save_2x_2D(TH2F *h1, TH2F *h2, TString sss, string ss, int hv_,int sn_, TString style){
        TCanvas *c3 = new TCanvas("c3","The Ntuple canvas",200,10,780,780);
        TPad *pad3 = new TPad("pad3","This is pad3",0.05,0.05,0.98,0.55,21);
        TPad *pad5 = new TPad("pad5","This is pad5",0.05,0.55,0.98,0.98,21);
        pad3->UseCurrentStyle(); pad5->UseCurrentStyle();
        pad3->Draw(); pad5->Draw();
        pad3->cd();
        pad3->SetGrid();
        h1->Draw(style);
        pad5->cd();
        pad5->SetGrid();
        h2->Draw(style);
        c3->cd();
        c3->Update();
        TString s("");
        s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.root",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        c3->SaveAs(s);
        s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.pdf",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        c3->SaveAs(s);
        s.Form("ScanId_%d/HV%d/_HV_%d_SN_%d_%s%s.png",sn_,hv_,hv_,sn_,sss.Data(),ss.c_str());
        c3->SaveAs(s);
}

