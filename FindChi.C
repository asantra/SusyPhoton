{
  TFile f1("ReweightedOldHLTNewMETMap.root");
  int nBins  = 15 ;
  float chi2 = 0;

  //// rebinning
  h_met_eeSample_diEMPtRhoReweighted->Rebin(4);
  h_met_ffSample_RhoReweighted->Rebin(4);
  h_met_candidate->Rebin(4);
  

  //// adding two histogram
  TH1F *eeffMET = new TH1F("eeffMET", "ee and ff MET added together", 200, 0., 200.);
  

  ///// normalization
  Float_t t1=h_met_eeSample_diEMPtRhoReweighted->Integral(0,50);
  h_met_eeSample_diEMPtRhoReweighted->Scale(1./t1);
  Float_t t2=h_met_ffSample_RhoReweighted->Integral(0,50);
  h_met_ffSample_RhoReweighted->Scale(1./t2);
  Float_t t3=h_met_candidate->Integral(0,50);
  h_met_candidate->Scale(1./t3);
  Float_t t4=h_met_eeSample_PurityReweighted->Integral(0,200);
  h_met_eeSample_PurityReweighted->Scale(1./t4);
  Float_t t5=h_met_ffSample_PurityReweighted->Integral(0,200);
  h_met_ffSample_PurityReweighted->Scale(1./t5);
  

  eeffMET->Add(h_met_eeSample_PurityReweighted);
  eeffMET->Add(h_met_ffSample_PurityReweighted);
  eeffMET->Rebin(4);
  Float_t t6=eeffMET->Integral(0,50);
  eeffMET->Scale(1./t6);

  /// setting color
  h_met_candidate->SetLineColor(kGreen+3);
  h_met_eeSample_diEMPtRhoReweighted->SetLineColor(4);
  h_met_ffSample_RhoReweighted->SetLineColor(2);
  eeffMET->SetLineColor(2);
  
  
  //h_met_ffSample_RhoReweighted->Draw();
  
  h_met_eeSample_diEMPtRhoReweighted->Draw();
  h_met_ffSample_RhoReweighted->Draw("sames");
  //eeffMET->Draw("sames");

  for(int i=1; i<=nBins; ++i){
    float x   = h_met_ffSample_RhoReweighted->GetBinContent(i);
    //printf("x:%f for bin %d\n",x,i);
    float dx  = h_met_ffSample_RhoReweighted->GetBinError(i);
    //printf("dx:%f for bin %d\n",dx,i);
    float y   = h_met_eeSample_diEMPtRhoReweighted->GetBinContent(i);
    //printf("y:%f for bin %d\n",y,i);
    float dy  = h_met_eeSample_diEMPtRhoReweighted->GetBinError(i);
    //printf("dy:%f for bin %d\n",dy,i);
    if(sqrt(dx*dx+dy*dy)!=0){
      float chi = (x-y)/(sqrt(dx*dx+dy*dy));
      printf("%f,\n",chi);
      chi2     += chi*chi;
    }
  }
  printf("Chi2: %f for Bins upto: %d\n",chi2,nBins);
  printf("Chi2/dof %f for Bins upto: %d\n",chi2/nBins,nBins);
  

}
