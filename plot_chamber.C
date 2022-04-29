void plot_chamber(){
/*	TFile _file0g("ScanId_768/HV6/_HV_6_SN_768_gammas_X_Y.root","read");
	TCanvas* ng =(TCanvas*)_file0g.Get("c4");
	TPad* nng =(TPad*)ng->GetPrimitive("pad4");
	TH2F *mmg = (TH2F*)nng->GetPrimitive("h2_XY");
        mmg->SetMarkerStyle(20);
        mmg->SetMarkerSize(0.25);
        mmg->SetStats(0);
        mmg->Draw("P");

*/

        TFile _file0m("ScanId_768/HV6/_HV_6_SN_768_muons_X_Y.root","read");
        TCanvas* nm =(TCanvas*)_file0m.Get("c4");
        TPad* nnm =(TPad*)nm->GetPrimitive("pad4");
        TH2F *mmm = (TH2F*)nnm->GetPrimitive("h2_XY");

        TCanvas *c2 = new TCanvas("c4n","The Ntuple canvas",0,0,1000,1500);
//        mmm->SetMarkerStyle(20);
//        mmm->SetMarkerSize(0.25);
        mmm->SetStats(0);
//	mmm->SetMarkerColor(4);

//        mmg->Draw("P");
        mmm->Draw("P ");


	Double_t x[7] = {0.,279.,340.,539.,503.,0.,0.};
	Double_t y[7] = {1220.,1220.,1168.,43.,0.,0.,1220.};
	TPolyLine *pline = new TPolyLine(7,x,y);
	pline->SetFillColor(0);
	pline->SetLineColor(2);
	pline->SetLineWidth(4);
//	pline->Draw("f");
	pline->Draw();

	c2->SaveAs("muon_strips_chamber.png");
}
