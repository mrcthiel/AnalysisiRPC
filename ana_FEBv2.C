#define ana_FEBv2_cxx
#include "ana_FEBv2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
using namespace std;

void ana_FEBv2::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ana_FEBv2.C
//      root> ana_FEBv2 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   TString         s("");
   s.Form("will save plots for RUN: , _HV_%d_SN_%d_MaxTrig_%d ",hv_,sn_,mt_);

   std::cout<<s<<std::endl;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout<<"N Events: " << nentries << std::endl;
   //Long64_t nentries = 100;

   Long64_t nbytes = 0, nb = 0;
   uint64_t t_bc0[3],trig[3];
   bool debug=false;
   float y = 0;
   float last_bc0 = 0;
   float run_duration_in_sec = 0;
   int any_strip_fired=0;
   int trig_exists=0;
   int ntrig_all=0;
   int ntrig_allevent=0;
   int ntrig_signal=0;
   Bool_t histo_exists=0;
   int nstrip_fired=0; 
   int Rstrip_fired=0; 
   int Dstrip_fired=0; 
   int Rstrip_bkg_fired=0; 
   int Dstrip_bkg_fired=0; 
   int nRstrip_fired=0; 
   int nDstrip_fired=0; 
   int nRstrip_bkg_fired=0; 
   int nDstrip_bkg_fired=0; 
   Long64_t trig_index = 0;
   float return_SR_int = 0;
   float return_BR_int = 0;
   float direct_SR_int = 0;
   float direct_BR_int = 0;


   TH1F* hLR = new TH1F("hLR", "LR", 50, 0.0, 50);
   TH1F* hHR = new TH1F("hHR", "HR", 50, 0.0, 50);
   TH2F* hdeltaT = new TH2F("deltaT", "T (HR - LR)", 50,0,50,1000, -500, 500);
   TH2F* hsumT = new TH2F("sumT", "T (HR + LR)", 50,0,50,20000, 0, 200000);
   TH2F* hHRT = new TH2F("hHRT", "T (HR - Trig)", 50,0,50,20000, -10000, 10000);
   TH2F* hLRT = new TH2F("hLRT", "T (LR - Trig)", 50,0,50,20000, -10000, 10000);
   TH1F* hLRT_ = new TH1F("hLRT_","T (LR - Trig)" ,20000,-10000,10000);
   TH1F* hHRT_ = new TH1F("hHRT_","T (HR - Trig)" ,20000,-10000,10000);
   TH1F* hLRT_2 = new TH1F("hLRT_2","n any strip" ,2,0,2);
   TH1F* hHRT_2 = new TH1F("hHRT_2","n strip" ,200,0,200);
   s.Form("plots_gif/events/_HV_%d_SN_%d_MaxTrig_%d_.root",hv_,sn_,mt_);

   TFile *outputFile = new TFile(s,"RECREATE");
   ///Get Final bc0 to calculate overall noise
   fChain->Draw("bc0>>histbc0");
   TH1F *hbc0 = (TH1F*)gDirectory->Get("histbc0");
   last_bc0 = hbc0->GetXaxis()->GetXmax();
   std::cout<<"integral/3 : "<<hbc0->Integral()/3<<" " << last_bc0 <<" "<<hbc0->FindLastBinAbove(1)<<std::endl;
   run_duration_in_sec = (last_bc0*90*pow(10,-6));
   std::cout<<"The data taking time is : "<< run_duration_in_sec <<std::endl;
   //fChain->Draw("(1)>>histTrigger","(m_channel(frame)==33)");
   //TH1F *histTrig = (TH1F*)gDirectory->Get("histTrigger");
   TH2F *hHREvt[100]; //For histograms for each event
   TH2F *hLREvt[100];
   for (uint32_t i=0;i<100;i++) {
   hHREvt[i] = new TH2F(s, "T (return - Trig)", 50,0,50,80,-940,-860);
   hLREvt[i] = new TH2F(s, "T (direct - Trig)", 50,0,50,80,-940,-860);
   }
   ofstream myfile;
   s.Form("outputs/_HV_%d_SN_%d_MaxTrig_%d_.txt",hv_,sn_,mt_);
   myfile.open(s);
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //std::cout<<"Hello "<<jentry<<" bc0 "<<bc0<<std::endl;
      // if (Cut(ientry) < 0) continue;
      float delta[3][100] = {};
      float sum[3][100] = {};
      float time_hr[3][100] = {};
      float time_lr[3][100] = {};
      any_strip_fired=0;
      trig_exists=0;
      nstrip_fired=0; 
	Rstrip_fired=0; 
	Dstrip_fired=0; 
	Rstrip_bkg_fired=0; 
	Dstrip_bkg_fired=0; 
      for (uint32_t i=0;i<nframe;i++)
        {
	  if (m_channel(frame[i])==32) {
    		t_bc0[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
		}
	else if (m_channel(frame[i])==33){
		trig[m_fpga(frame[i])]=m_time(m_traw(frame[i]));
		//std::cout<< "trig time : " << trig[m_fpga(frame[i])] << "  bc0 time : "<< t_bc0[m_fpga(frame[i])] <<std::endl;
		//std::cout<<  int(trig[m_fpga(frame[i])] - t_bc0[m_fpga(frame[i])]) <<std::endl;
		if ((trig[m_fpga(frame[i])] - t_bc0[m_fpga(frame[i])]) > -85000 and (trig[m_fpga(frame[i])] - t_bc0[m_fpga(frame[i])]) < -500) {
			trig_exists=1;
		}
		if (debug) myfile<<"\t "<<"trig time "<< m_time(m_traw(frame[i])) <<"\n"; 
		}
	else {
		//if (m_time(m_traw(frame[i]))>220  and m_time(m_traw(frame[i]))<88900) {
		any_strip_fired = 1;
		if (c_side(m_channel(frame[i])) < 1){
			//direct stips
	                hLR->Fill(m_strip(frame[i]),1/(run_duration_in_sec*150));
			time_lr[m_fpga(frame[i])][m_channel(frame[i])] = m_time(m_traw(frame[i]));
		}
       		if (c_side(m_channel(frame[i])) > 0.1){
			//return strips
               		hHR->Fill(m_strip(frame[i]),1/(run_duration_in_sec*150));
			time_hr[m_fpga(frame[i])][m_channel(frame[i])] = m_time(m_traw(frame[i]));
			}
	//}
	}
        }

      //if (ntrig!=3) std::cout<< "nTRIG in this event : "<< ntrig  <<std::endl;
      if (any_strip_fired>0) hLRT_2->Fill(any_strip_fired);

      for (uint32_t i=0;i<nframe;i++)
      {
	//std::cout<< " bc0 "<<bc0<<" channel "<<m_channel(frame[i])<<std::endl;
	if (m_channel(frame[i])==32) continue;
	if (m_channel(frame[i])==33) continue;
	//if (m_time(m_traw(frame[i]))>220 ) {
        hHRT->Fill(m_strip(frame[i]),time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]);
        hLRT->Fill(m_strip(frame[i]),time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]);//}
        hHRT_->Fill(time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]);
        hLRT_->Fill(time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]);
	//Enter muon window
	if(trig_exists>0){
		ntrig_all++;
	if ((time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) > -930 and (time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) < -850 ){ 
	//std::cout<< " jentry "<<trig_index<<" filled "<<std::endl;
		ntrig_signal++;
		Rstrip_fired =1;
	if (trig_index < 100 ) hHREvt[int(trig_index)]->Fill(m_strip(frame[i]),time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]);
	}
	if ((time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) > -5000 and (time_hr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) < -1000 ) Rstrip_bkg_fired=1;
	if ((time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) > -930 and (time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) < -850 ){ 
		Dstrip_fired =1;
		nstrip_fired++;
	if (trig_index < 100 ) hLREvt[int(trig_index)]->Fill(m_strip(frame[i]),time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]);
	}
	if ((time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) > -5000 and (time_lr[m_fpga(frame[i])][m_channel(frame[i])]-trig[m_fpga(frame[i])]) < -1000 ) Dstrip_bkg_fired=1;
   			}
        
	if (m_channel(frame[i])<32&&m_channel(frame[i])>15&&c_side(m_channel(frame[i]))>0.1){
		//take channels connected to return channels only
		delta[m_fpga(frame[i])][m_channel(frame[i])] = time_hr[m_fpga(frame[i])][m_channel(frame[i])]-time_lr[m_fpga(frame[i])][m_channel(frame[i])-1];
		sum[m_fpga(frame[i])][m_channel(frame[i])] = time_hr[m_fpga(frame[i])][m_channel(frame[i])]+time_lr[m_fpga(frame[i])][m_channel(frame[i])-1];
		}
	if (m_channel(frame[i])<15&&c_side(m_channel(frame[i]))>0.1){
		//take channels connected to return channels only
		delta[m_fpga(frame[i])][m_channel(frame[i])] = time_hr[m_fpga(frame[i])][m_channel(frame[i])]-time_lr[m_fpga(frame[i])][m_channel(frame[i])+1];
		sum[m_fpga(frame[i])][m_channel(frame[i])] = time_hr[m_fpga(frame[i])][m_channel(frame[i])]+time_lr[m_fpga(frame[i])][m_channel(frame[i])+1];
		}
	if (c_side(m_channel(frame[i]))<1){delta[m_fpga(frame[i])][m_channel(frame[i])] =0;}

	hdeltaT->Fill(m_strip(frame[i]),delta[m_fpga(frame[i])][m_channel(frame[i])]);
	hsumT->Fill(m_strip(frame[i]),sum[m_fpga(frame[i])][m_channel(frame[i])]);

      }
	if(trig_exists>0){
		ntrig_allevent++;
		if (Rstrip_fired) nRstrip_fired++;
		if (Dstrip_fired) nDstrip_fired++;
		if (Rstrip_bkg_fired) nRstrip_bkg_fired++;
		if (Dstrip_bkg_fired) nDstrip_bkg_fired++;
      		hHRT_2->Fill(nstrip_fired);
		if (trig_index < 100){
		hHREvt[int(trig_index)]->Write();
		hLREvt[int(trig_index)]->Write();
		trig_index++;
		}
	}
   }
outputFile->Close();
std::cout<<"any strip fired: " << hLRT_2->Integral() << std::endl;
return_SR_int = hHRT_->Integral(hHRT_->FindBin(-940),hHRT_->FindBin(-880));
return_BR_int = ((hHRT_->Integral(hHRT_->FindBin(-3000),hHRT_->FindBin(-1000)))/100)+((hHRT_->Integral(hHRT_->FindBin(1000),hHRT_->FindBin(3000)))/100);
direct_SR_int = hLRT_->Integral(hLRT_->FindBin(-940),hLRT_->FindBin(-880));
direct_BR_int = ((hLRT_->Integral(hLRT_->FindBin(-3000),hLRT_->FindBin(-1000)))/100)+((hLRT_->Integral(hLRT_->FindBin(1000),hLRT_->FindBin(3000)))/100);
myfile<<"purity return: " << (return_SR_int-return_BR_int)/return_SR_int <<"\n";
myfile<<"purity direct: " << (direct_SR_int-direct_BR_int)/direct_SR_int <<"\n";
myfile<<"Count narrow return: " << return_SR_int << " bkg : " << return_BR_int <<"\n";
myfile<<"Count narrow direct: " << direct_SR_int << " bkg : " << direct_BR_int <<"\n";
myfile<<"Triggered events: " << ntrig_allevent <<"\n";
myfile<<"at least one return fired in muon window: " << nRstrip_fired <<"\n";
myfile<<"at least one direct fired in muon window: " << nDstrip_fired <<"\n";
myfile<<"at least one return fired in bkg window 50x(muon window) norm count: " << nRstrip_bkg_fired <<"\n";
myfile<<"at least one direct fired in bkg window 50x(muon window) norm count: " << nDstrip_bkg_fired <<"\n";
myfile<<"Efficiency direct: " << (nDstrip_fired-(nDstrip_bkg_fired/50))/float(ntrig_allevent-(nDstrip_bkg_fired/50)) <<"\n";
myfile<<"Efficiency return: " << (nRstrip_fired-(nRstrip_bkg_fired/50))/float(ntrig_allevent-(nRstrip_bkg_fired/50)) <<"\n";
myfile.close();
ofstream mydcsfile;
s.Form("outputs/_HV_%d_SN_%d_MaxTrig_%d__forwebDCS.txt",hv_,sn_,mt_);
mydcsfile.open(s);
mydcsfile<<(nDstrip_fired-(nDstrip_bkg_fired/50))/float(ntrig_allevent-(nDstrip_bkg_fired/50));
mydcsfile.close();



gStyle->SetOptStat(0);
TCanvas *c1 = new TCanvas("c1","The Ntuple canvas",200,10,780,780);
TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.05,0.98,0.55,21);
TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.55,0.98,0.98,21);
pad1->UseCurrentStyle(); pad2->UseCurrentStyle();
pad1->Draw(); pad2->Draw();
pad1->cd();
pad1->SetGrid();
pad1->SetLogy();
hLR->Draw("Hist");
hLR->GetYaxis()->SetTitle("Noise Hz/cm");
hLR->GetXaxis()->SetTitle("Strips");
pad2->cd();
pad2->SetGrid();
pad2->SetLogy();
hHR->Draw("Hist");
hHR->GetYaxis()->SetTitle("Noise Hz/cm");
hHR->GetXaxis()->SetTitle("Strips");
c1->cd();
c1->Update();
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__avarage_noise_side.root",hv_,sn_,mt_);
c1->SaveAs(s);
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__avarage_noise_side.pdf",hv_,sn_,mt_);
c1->SaveAs(s);
TCanvas *c2 = new TCanvas("c4","The Ntuple canvas",200,10,900,780);
TPad *pad4 = new TPad("pad4","This is pad4",0.02,0.02,0.98,0.98,21);
pad4->UseCurrentStyle();
pad4->Draw();
pad4->cd();
pad4->SetGrid();
pad4->SetLogz();
//pad4->SetLogy();
pad4->SetRightMargin(0.12);
hdeltaT->GetYaxis()->SetTitle("delta time [ns]");
hdeltaT->GetXaxis()->SetTitle("Strips");
hdeltaT->GetYaxis()->SetLabelOffset(0.007);
hdeltaT->GetYaxis()->SetTitleOffset(1.5);
hdeltaT->SetLabelSize(0.02,"x");
hdeltaT->SetTitle("");
hdeltaT->Draw("COLZ");
c2->cd();
c2->Update();
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__deltaT_srip.root",hv_,sn_,mt_);
c2->SaveAs(s);
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__deltaT_srip.pdf",hv_,sn_,mt_);
c2->SaveAs(s);
TCanvas *c3 = new TCanvas("c3","The Ntuple canvas",200,10,780,780);
TPad *pad3 = new TPad("pad3","This is pad3",0.05,0.05,0.98,0.55,21);
TPad *pad5 = new TPad("pad5","This is pad5",0.05,0.55,0.98,0.98,21);
pad3->UseCurrentStyle(); pad5->UseCurrentStyle();
pad3->Draw(); pad5->Draw();
pad3->cd();
pad3->SetGrid();
hLRT->Draw("COLZ");
pad5->cd();
pad5->SetGrid();
hHRT->Draw("COLZ");
c3->cd();
c3->Update();
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side.root",hv_,sn_,mt_);
c3->SaveAs(s);
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side.pdf",hv_,sn_,mt_);
c3->SaveAs(s);
TCanvas *c10 = new TCanvas("c10","The Ntuple canvas",200,10,780,780);
TPad *pad10 = new TPad("pad10","This is pad10",0.05,0.05,0.98,0.55,21);
TPad *pad11 = new TPad("pad11","This is pad11",0.05,0.55,0.98,0.98,21);
pad10->UseCurrentStyle(); pad11->UseCurrentStyle();
pad10->Draw(); pad11->Draw();
pad10->cd();
pad10->SetGrid();
hLRT_->Draw("Hist");
pad11->cd();
pad11->SetGrid();
hHRT_->Draw("Hist");
c10->cd();
c10->Update();
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side_1d.root",hv_,sn_,mt_);
c10->SaveAs(s);
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__DeltaTrig_side_1d.pdf",hv_,sn_,mt_);
c10->SaveAs(s);
TCanvas *c11 = new TCanvas("c11","The Ntuple canvas",200,10,780,780);
TPad *pad20 = new TPad("pad20","This is pad10",0.05,0.05,0.98,0.55,21);
TPad *pad21 = new TPad("pad21","This is pad11",0.05,0.55,0.98,0.98,21);
pad20->UseCurrentStyle(); pad21->UseCurrentStyle();
pad20->Draw(); pad21->Draw();
pad20->cd();
pad20->SetGrid();
hLRT_2->Draw("Hist");
pad21->cd();
pad21->SetGrid();
hHRT_2->Draw("Hist");
c11->cd();
c11->Update();
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__Time_side_1d.root",hv_,sn_,mt_);
c11->SaveAs(s);
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__Time_side_1d.pdf",hv_,sn_,mt_);
c11->SaveAs(s);

TCanvas *csum = new TCanvas("csum","The Ntuple canvas",200,10,900,780);
TPad *padsum = new TPad("padsum","This is padsum",0.02,0.02,0.98,0.98,21);
padsum->UseCurrentStyle();
padsum->Draw();
padsum->cd();
padsum->SetGrid();
padsum->SetLogz();
//pad4->SetLogy();
padsum->SetRightMargin(0.12);
hsumT->GetYaxis()->SetTitle("sum time [ns]");
hsumT->GetXaxis()->SetTitle("Strips");
hsumT->GetYaxis()->SetLabelOffset(0.007);
hsumT->GetYaxis()->SetTitleOffset(1.5);
hsumT->SetLabelSize(0.02,"x");
hsumT->SetTitle("");
hsumT->Draw("COLZ");
csum->cd();
csum->Update();
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__sumT_srip.root",hv_,sn_,mt_);
csum->SaveAs(s);
s.Form("plots_gif/_HV_%d_SN_%d_MaxTrig_%d__sumT_srip.pdf",hv_,sn_,mt_);
csum->SaveAs(s);



}
