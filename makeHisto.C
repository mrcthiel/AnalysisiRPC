 
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"


void makeHisto(Int_t rn) {
	TString s("");
	//s.Form("/data/beamdump/root_trees/_HV_-999_SN_-999_MaxTrig_-999_Run_%d.root",rn);
	//s.Form("/data/beamdump/root_trees/tdc_time_calibration_FEB1.root");
	s.Form("/data/904/root_trees/_HV_-999_SN_-999_MaxTrig_-999_Run_1000.root");
	//s.Form("/data/root_trees/Run_1167.root");
	std::cout<<s<<std::endl;
	TFile *f1 = TFile::Open(s);
	TTree *evt = (TTree*)f1->Get("evt");
	gROOT->SetStyle("Plain");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	//Draw bc0 time 
	TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",200,10,780,780);
	TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.02,0.98,0.98,21);
	pad1->UseCurrentStyle();
	pad1->Draw();
	pad1->cd();
	pad1->SetGrid();
	pad1->SetLogy();
	TH1F *hist1= new TH1F("hist1", "bc0 time",100,89624.6,89625.4);
	evt->Draw("m_time(m_traw(frame))>>hist1","(m_channel(frame)==32)");
	std::cout<<"All events / 3: "<< (hist1->Integral())/3 << std::endl;
	hist1->SetLineColor(4);
	hist1->SetLineWidth(2);
	hist1->GetXaxis()->SetTitle("bc0 time [ns]");
	hist1->GetYaxis()->SetTitle("Events");
	hist1->GetYaxis()->SetLabelOffset(0.007);
	hist1->GetYaxis()->SetTitleOffset(1.3);
	hist1->SetLabelSize(0.02,"x");
	hist1->SetTitle("");
	hist1->Draw("Histo");
	c1->cd();
	c1->Update();
	s.Form("plots_904/run_%d_bc_0.root",rn);
	c1->SaveAs(s);
	s.Form("plots_904/run_%d_bc_0.pdf",rn);
	c1->SaveAs(s);
	//Draw trigger time
	TCanvas *c2 = new TCanvas("c2","The Ntuple canvas",200,10,780,780);
	TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.98,0.98,21);
	pad2->UseCurrentStyle();
	pad2->Draw();
	pad2->cd();
	pad2->SetGrid();
	pad2->SetLogy();
	TH1F *hist2= new TH1F("hist2", "trigger time",91000,0,91000);
	evt->Draw("m_time(m_traw(frame))>>hist2","(m_channel(frame)==33)");
	std::cout<<"All triggers /3 : "<< (hist2->Integral())/3 << std::endl;
	hist2->SetLineColor(4);
	hist2->SetLineWidth(2);
	hist2->GetXaxis()->SetTitle("trigger time [ns]");
	hist2->GetYaxis()->SetTitle("Events");
	hist2->GetYaxis()->SetLabelOffset(0.007);
	hist2->GetYaxis()->SetTitleOffset(1.3);
	hist2->SetLabelSize(0.02,"x");
	hist2->SetTitle("");
	hist2->Draw("Histo");
	c2->cd();
	c2->Update();
	s.Form("plots_904/run_%d_triggerTime.root",rn);
	c2->SaveAs(s);
	s.Form("plots_904/run_%d_triggerTime.pdf",rn);
	c2->SaveAs(s);
	//time vs strip	
	TCanvas *c3 = new TCanvas("c3","The Ntuple canvas",200,10,900,780);
	TPad *pad3 = new TPad("pad3","This is pad3",0.02,0.02,0.98,0.98,21);
	pad3->UseCurrentStyle();
	pad3->Draw();
	pad3->cd();
	pad3->SetGrid();
	pad3->SetLogz();
	TH2F *hist3= new TH2F("hist3", "direct strip time",91000,0,91000,50,0,50);
	evt->Draw("m_strip(frame):m_time(m_traw(frame))>>hist3","(m_channel(frame)<32)&&(c_side(m_channel(frame)))<0.5");
	hist3->SetLineColor(4);
	hist3->SetLineWidth(2);
	hist3->GetXaxis()->SetTitle("time stamp");
	hist3->GetYaxis()->SetTitle("Strips");
	hist3->GetYaxis()->SetLabelOffset(0.007);
	hist3->GetYaxis()->SetTitleOffset(1.3);
	hist3->SetLabelSize(0.02,"x");
	hist3->SetTitle("");
	hist3->Draw("COLZ");
	c3->cd();
	c3->Update();
	s.Form("plots_904/run_%d_2D_timeDirectvsStrip.root",rn);
	c3->SaveAs(s);
	s.Form("plots_904/run_%d_2D_timevsDirectStrip.pdf",rn);
	c3->SaveAs(s);
	// strip vs bc0
	TCanvas *c4 = new TCanvas("c4","The Ntuple canvas",200,10,900,780);
	TPad *pad4 = new TPad("pad4","This is pad4",0.02,0.02,0.98,0.98,21);
	pad4->UseCurrentStyle();
	pad4->Draw();
	pad4->cd();
	pad4->SetGrid();
	pad4->SetLogz();
	pad4->SetRightMargin(0.12);
	TH2F *hist4= new TH2F("hist4", "trigger time",100,0,15,50,0,50);
	evt->Draw("m_strip(frame):(bc0*90*10^(-6))/3600>>hist4","m_channel(frame)<32&&m_time(m_traw(frame))>0");
	hist4->SetLineColor(4);
	hist4->SetLineWidth(2);
	hist4->GetXaxis()->SetTitle("time [hours]");
	hist4->GetYaxis()->SetTitle("Strips");
	hist4->GetYaxis()->SetLabelOffset(0.007);
	hist4->GetYaxis()->SetTitleOffset(1.3);
	hist4->SetLabelSize(0.02,"x");
	hist4->SetTitle("");
	hist4->Draw("COLZ");
	c4->cd();
   	c4->Update();
	s.Form("plots_904/run_%d_2D_bc0timevsStrips.root",rn);
