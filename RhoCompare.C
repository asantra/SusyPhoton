{
   TFile f1("SymmetricCheckHLTDiff.root");
   /*Float_t x1=h_rho_candidate->Integral(0,50);
   //h_rho_candidate->Sumw2();
   h_rho_candidate->Scale(1./x1);
   h_rho_candidate->Rebin(2);
   //h_met_candidate->SetLineColor(kGreen+3);
   h_rho_candidate->SetMarkerStyle(kFullDotLarge);
   h_rho_candidate->SetMarkerColor(kGreen+1);
   Float_t x2=h_rho_ffSample->Integral(0,50);
   h_rho_ffSample->Scale(1./x2);
   h_rho_ffSample->Rebin(2);
   //h_met_ffSample_Rho25Reweighted->SetLineColor(2);
   h_rho_ffSample->SetMarkerStyle(kFullSquare);
   h_rho_ffSample->SetMarkerColor(kRed);
   Float_t x3=h_rho_eeSample->Integral(0,50);
   h_rho_eeSample->Scale(1./x3);
   h_rho_eeSample->Rebin(2);
   //h_met_eeSample_Rho25Reweighted->SetLineColor(4);
   h_rho_eeSample->SetMarkerStyle(kFullTriangleUp);
   h_rho_ffSample->SetMarkerColor(kBlue);
   
   h_rho_ffSample->Draw();
   h_rho_candidate->Draw("sames");
   h_rho_eeSample->Draw("sames");*/
   Float_t y1=h_diEMPt_candidate_HLT->Integral(0,200);
   h_diEMPt_candidate_HLT->Scale(1./y1);
   h_diEMPt_candidate_HLT->SetMarkerStyle(kFullDotLarge);
   h_diEMPt_candidate_HLT->SetMarkerColor(kGreen+1);
   h_diEMPt_candidate_HLT->Rebin(2);
   Float_t y2=h_diEMPt_eeSample_HLT->Integral(0,200);
   h_diEMPt_eeSample_HLT->Scale(1./y2);
   h_diEMPt_eeSample_HLT->Rebin(2);
   h_diEMPt_eeSample_HLT->SetMarkerStyle(kFullTriangleUp);
   h_diEMPt_eeSample_HLT->SetMarkerColor(kBlue);
   h_diEMPt_eeSample_HLT->Draw();
   h_diEMPt_candidate_HLT->Draw("sames");
   Float_t y3=h_diEMPt_ffSample_HLT->Integral(0,200);
   h_diEMPt_ffSample_HLT->Scale(1./y3);
   h_diEMPt_ffSample_HLT->Rebin(2);
   h_diEMPt_ffSample_HLT->SetMarkerStyle(kFullSquare);
   h_diEMPt_ffSample_HLT->SetMarkerColor(kRed);
   h_diEMPt_ffSample_HLT->Draw("sames");
   

}
