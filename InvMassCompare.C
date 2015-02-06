{
  TFile f1("PhotonAll_InvMassNewTry.root");
  Float_t x1=h_gg_invmass_all_passing_shower_cut->Integral(0,150);
  h_gg_invmass_all_passing_shower_cut->Scale(1./x1);
  h_gg_invmass_all_passing_shower_cut->SetLineColor(2);
  Float_t x2=h_gg_invmass_all_failing_shower_cut->Integral(0,150);
  h_gg_invmass_all_failing_shower_cut->Scale(1./x2);
  h_gg_invmass_all_failing_shower_cut->SetLineColor(4);
  
  h_gg_invmass_all_passing_shower_cut->Draw();
  h_gg_invmass_all_failing_shower_cut->Draw("sames");
  
}
