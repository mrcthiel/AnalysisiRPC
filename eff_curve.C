void eff_curve(Int_t sn ){


  TGraphErrors* efficiency = new TGraphErrors();

  // open a file in read mode.

  string fileInName("");
  //fileInName = fileInName + "" + fileIn + ".txt";
  //fileInName = "eff_input.txt";
  //fileInName = Form("eff_input_SN%d.txt", sn);
  //fileInName = Form("eff_input_medium_SN%d", sn);
  fileInName = Form("eff_input_loose_SN156.txt");

  cout << fileInName.c_str() << endl;
  string sHV, sEff, sEff_Err, line, sTitle;
  float HV, Eff, Eff_Err, HVmin, HVmax;


  int i = -2;
  ifstream readFile(fileInName.c_str());


  
  while(getline(readFile,line))   {
    if (i == -2){
      i++;
      sTitle=line;
      continue;
    }
      i++;
      stringstream iss(line);
      getline(iss, sHV, '\t');
      getline(iss, sEff, '\t');
      getline(iss, sEff_Err, '\t');
      
      istringstream fsHV(sHV); fsHV >> HV;
      istringstream fsEff(sEff); fsEff >> Eff;
      istringstream fsEff_Err(sEff_Err); fsEff_Err >> Eff_Err;
      
      cout << "At " << sHV.c_str() << " eff = " <<  sEff.c_str() << " #pm " << sEff_Err.c_str() << endl;
      
      efficiency->SetPoint(i,HV, Eff/100.);
      efficiency->SetPointError(i,1e-3, sqrt(Eff/100*(1-Eff/100))/sqrt(1000));
      if (i == 0) HVmin = HV;

  }

  HVmax = HV;

  readFile.close();
  

 //TF1* sigmoid = new TF1("sigmoid","[0]/(1+exp([1]*([2]-x)))",HVmin-100,HVmax+100);
 TF1* sigmoid = new TF1("sigmoid","[0]/(1+exp([1]*([2]-x)))",6000,7200);
 sigmoid->SetParName(0,"#epsilon_{max}");
 sigmoid->SetParName(1,"#lambda");
 sigmoid->SetParName(2,"HV_{50%}");
 sigmoid->SetParameter(0,0.98);
 sigmoid->SetParameter(1,0.01);
 sigmoid->SetParameter(2,7000);

            //****************************************************


efficiency->Fit(sigmoid);

efficiency->SetMarkerStyle(22);
efficiency->SetMarkerSize(2);

 sTitle = sTitle + "; HVeff; Eff";
 TH1D* Plotter = new TH1D("Plotter", sTitle.c_str(), 1, HVmin-100,HVmax+100);

 Plotter->SetStats(0);

 Plotter->SetMinimum(0.);
 Plotter->SetMaximum(1.08);

 Plotter->Draw();

efficiency->Draw("P");


 //gPad->SetLogy(1);

 

 TLegend *leg = new TLegend(0.608739,0.1024784,0.8085924,0.2523343,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.04);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);

  
  leg->AddEntry(efficiency,"Efficiency","P");

  //  leg->Draw();


  	
  double p1 = sigmoid->GetParameter(0);
  double p2 = sigmoid->GetParameter(1);
  double p3 = sigmoid->GetParameter(2);
   	

  double uLimit = HVmax, lLimit = HVmin;

  TLatex* ltx = new TLatex();
  ltx->SetTextSize(0.04);

  double knee = p3 - log(1/0.95-1)/p2;


  
  double WP = knee+120;
  cout << "Val at WP " << WP << endl;

  
  TLine* lWP = new TLine(WP, 0., WP, 1);
  lWP->SetLineStyle(2);
  lWP->Draw();

  double shift = 0.5;
  
  double add = (uLimit-lLimit)/11., up = 0.3; 
  ltx->DrawLatex(WP+shift*add, 0.42+up, Form("Eff(WP) = %.2f", sigmoid->Eval(WP)));
  ltx->DrawLatex(WP+shift*add, 0.35+up, Form("WP = %.0f V", WP));
  ltx->DrawLatex(WP+shift*add, 0.27+up, Form("knee = %.0f V", knee));
  ltx->DrawLatex(WP+shift*add, 0.20+up, Form("HV 50% = %.0f V", p3));


  TLine* plateau = new TLine(lLimit-50, p1, uLimit+50, p1);
  plateau->SetLineStyle(2);
  plateau->Draw();

  add = (uLimit-lLimit)/11.;
  
  if ((knee - lLimit) < (uLimit-lLimit)*(3/11.)) add = knee + add;
  else add = lLimit+add;

  ltx->DrawLatex(add, p1+0.02, Form("plateau = %.2f", p1));



  cout << "knee = " << knee << endl;

  string outfile = "Efficiency";
  string outfile_png = Form("Efficiency_loose_50dac_Beamdump.png");	    //,sn);
  string outfile_pdf = Form("Efficiency_loose_50dac_Beamdump.pdf");	    //,sn);
  string outfile_root =Form("Efficiency_loose_50dac_Beamdump.root");	    //,sn);

  gPad->SaveAs(outfile_png.c_str());
  gPad->SaveAs(outfile_pdf.c_str());

  TFile* effout = new TFile(outfile_root.c_str(), "RECREATE");
  effout->cd();
  efficiency->Write("Efficiency");
  effout->Close();  



}

