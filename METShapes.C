{
  TFile f1("PhotonAll_MET.root");
  Float_t x1=h_met_candidate->Integral(0,200);
  h_met_candidate->Scale(1./x1);
  h_met_candidate->SetLineColor(1);
 
  h_met_eeSample->Sumw2();
  h_met_eeSample->Rebin(4);
  Float_t x2=h_met_eeSample->Integral(0,50);
  h_met_eeSample->Scale(1./x2);
  h_met_eeSample->SetLineColor(2);
  
  h_met_ffSample->Sumw2();
  h_met_ffSample->Rebin(4);
  Float_t x3=h_met_ffSample->Integral(0,50);
  h_met_ffSample->Scale(1./x3);
  h_met_ffSample->SetLineColor(4);

  //h_met_candidate->Draw();
  TH1F *eeOverffMET = new TH1F("eeOverffMET", "ee MET to ff MET ratio", 50, 0.0, 200.0);
  //h_met_eeSample->Draw();
  //h_met_ffSample->Draw("sames");
  eeOverffMET->Sumw2();
  eeOverffMET->Divide(h_met_eeSample, h_met_ffSample);
  eeOverffMET->Draw();
  
}
