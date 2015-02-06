{  
   TFile f2("ReRecoCheck.root");
   TH1F *h_metY_rereco = (TH1F*)h_metY_rereco->Clone("h_metY_rereco");
   Float_t Y1=h_metY_rereco->Integral(0,400);
   h_metY_rereco->Scale(1.0/Y1);
   h_metY_rereco->SetLineColor(4);

  
   
   
   TFile f3("PromptRecoCheck.root");
   TH1F *h_metY_prompt=(TH1F*)h_metY_prompt->Clone("h_metY_prompt");
   Float_t Y2=h_metY_prompt->Integral(0,400);
   h_metY_prompt->Scale(1.0/Y2);
   h_metY_prompt->SetLineColor(2);
   //h_nvtY_weight->Sumw2();
   h_metY_prompt->Draw();
   h_metY_rereco->Draw("sames");

}
