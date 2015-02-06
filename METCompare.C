{
   TFile f1("SymmetricCheckHLTDiff.root");
   Float_t x1=h_met_candidate_HLT->Integral(0,200);
   //h_met_candidate->Sumw2();
   h_met_candidate_HLT->Scale(1./x1);
   h_met_candidate_HLT->Rebin(4);
   
   h_met_candidate_HLT->SetMarkerStyle(kFullDotLarge);
   h_met_candidate_HLT->SetMarkerColor(kGreen+1);
   Float_t x2=h_met_ffSample_HLT->Integral(0,200);
   h_met_ffSample_HLT->Scale(1./x2);
   h_met_ffSample_HLT->Rebin(4);
   //h_met_ffSample_Rho25Reweighted->SetLineColor(2);
   h_met_ffSample_HLT->SetMarkerStyle(kFullSquare);
   h_met_ffSample_HLT->SetMarkerColor(kRed);
   Float_t x3=h_met_eeSample_HLT->Integral(0,200);
   h_met_eeSample_HLT->Scale(1./x3);
   h_met_eeSample_HLT->Rebin(4);
   //h_met_eeSample_Rho25Reweighted->SetLineColor(4);
   h_met_eeSample_HLT->SetMarkerStyle(kFullTriangleUp);
   h_met_ffSample_HLT->SetMarkerColor(kBlue);
   h_met_candidate_HLT->Draw();
   h_met_ffSample_HLT->Draw("sames");
   h_met_eeSample_HLT->Draw("sames");


}
