{
   TFile f1("ReweightedOldHLT.root");
   Float_t fitchi2;
   Int_t NDF;
   TObjArray *component = new TObjArray(2);        // MC histograms are put in this array

   
   //Float_t x1=h_met_eeSample->Integral(0,200);
   //Float_t x2=h_met_ffSample_nVtxReweighted->Integral(0,200);
   //Float_t x3=h_met_candidate->Integral(0,200);
   //h_met_eeSample->Scale(1./x1);
   //h_met_ffSample->Scale(1./x2);
   //h_met_candidate->Scale(1./x3);

   component->Add(h_met_eeSample_diEMPtRhoReweighted);
   component->Add(h_met_ffSample_RhoReweighted);
   
   TFractionFitter* fit = new TFractionFitter(h_met_candidate, component); // initialise
   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(0,50);                    // use only the first 15 bins in the fit
   Int_t status = fit->Fit();              // perform the fit
   fitchi2=fit->GetChisquare();
   NDF=fit->GetNDF();
   std::cout << "fit status: " << status << std::endl;
   std::cout << "chi2:       " << fitchi2 << std::endl;
   std::cout << "chi2/ndf:   " << fitchi2/NDF << std::endl;
   TH1F *CandidateOverFit = new TH1F("CandidateOverFit", "Candidate diEMPt to fit diEMPt ratio", 200, 0.0, 200.0);
   if (status == 0) {                       // check on fit status
     TH1F* result = (TH1F*) fit->GetPlot();
     h_met_candidate->SetLineColor(4);
     h_met_candidate->Draw("Ep");
     result->SetLineColor(2);
     result->Draw("same");
     //TH1F *CandidateOverFit = new TH1F("CandidateOverFit", "Candidate MET to fit MET ratio", 200, 0.0, 200.0);
     CandidateOverFit->Sumw2();
     CandidateOverFit->Divide(h_met_candidate, result); 
     //CandidateOverFit->Draw();
   }
 }
