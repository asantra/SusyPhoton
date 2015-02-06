 {
   TFile f1("ImportantReReco.root");
   h_met_eeSample_diEMPtRhoReweighted->Sumw2();
   h_met_eeSample_diEMPtRhoReweighted->Rebin(4);
   Float_t x1=h_met_eeSample_diEMPtRhoReweighted->Integral(0,50);
   h_met_eeSample_diEMPtRhoReweighted->Scale(1.0/x1);
   h_met_eeSample_diEMPtRhoReweighted->SetLineColor(4);
   //h_met_eeSample_diEMPtRhoReweighted->GetXaxis()->SetRangeUser(0,80);



   //h_met_eeSample->Sumw2();
   //h_met_eeSample_PurityReweighted->Rebin(4);
   Float_t x11=h_met_eeSample_PurityReweighted->Integral(0,200);
   h_met_eeSample_PurityReweighted->Scale(1.0/x11);
   h_met_eeSample_PurityReweighted->SetLineColor(4);
   //h_met_eeSample->GetXaxis()->SetRangeUser(0,80);


   h_met_rereco_all->Sumw2();
   //h_met_rereco_all->Rebin(4);
   Float_t y11=h_met_rereco_all->Integral(0,200);
   h_met_rereco_all->Scale(1.0/y11);
   h_met_rereco_all->SetLineColor(kGreen+3);
   h_met_rereco_all->GetXaxis()->SetRangeUser(0,80);

   //TH1F *Fakerhoreweighted_rereco=(TH1F*)h_met_rereco_rhoreweighted->Clone("Fakerhoreweighted_rereco");

   //h_met_ffSample_RhoReweighted->Sumw2();
   //h_met_ffSample_RhoReweighted->Rebin(4);
   Float_t x2=h_met_ffSample_RhoReweighted->Integral(0,200);
   h_met_ffSample_RhoReweighted->Scale(1.0/x2);
   h_met_ffSample_RhoReweighted->SetLineColor(2);
   //h_met_ffSample_RhoReweighted->GetXaxis()->SetRangeUser(0,80);


   //h_met_ffSample->Sumw2();
   //h_met_ffSample_PurityReweighted->Rebin(4);
   Float_t x22=h_met_ffSample_PurityReweighted->Integral(0,200);
   //h_met_ffSample_PurityReweighted->Scale(1.0/x22);
   h_met_ffSample_PurityReweighted->SetLineColor(2);
   //h_met_ffSample->GetXaxis()->SetRangeUser(0,80);


   h_met_PhoFake_PurityReweighted->Sumw2();
   //h_met_PhoFake_PurityReweighted->Rebin(4);
   Float_t x3=h_met_PhoFake_PurityReweighted->Integral(0,200);
   //h_met_PhoFake_PurityReweighted->Scale(1.0/x3);
   h_met_PhoFake_PurityReweighted->SetLineColor(kCyan+3);
   //h_met_PhoFake_rereco->GetXaxis()->SetRangeUser(0,80);

   h_met_FakePho_PurityReweighted->Sumw2();
   //h_met_FakePho_PurityReweighted->Rebin(4);
   Float_t x4=h_met_FakePho_PurityReweighted->Integral(0,200);
   //h_met_FakePho_PurityReweighted->Scale(1.0/x4);
   h_met_FakePho_PurityReweighted->SetLineColor(kOrange-6);
   //h_met_FakePho_rereco->GetXaxis()->SetRangeUser(0,80);

  


   TH1F *h3 = new TH1F("h3","ee*0.466175+ff*0.182506+yf*0.260527+fy*0.09079",200,0,200);
   TH1F *h4 = new TH1F("h4","ee*0.466175+ff*0.533825",200,0,200);
   TH1F *h5 = new TH1F("h5","ee*0.466175+ff*0.182506+yf*0.260527",200,0,200);
   
   
   h3->Sumw2();
   h4->Sumw2();
   h5->Sumw2();
   
   h4->Add(h_met_eeSample_PurityReweighted, h_met_ffSample_RhoReweighted);
   //h4->Draw();
   h5->Add(h4, h_met_PhoFake_PurityReweighted);
   //h5->Draw();
   h3->Add(h5, h_met_FakePho_PurityReweighted);
   
   //h3->Draw();
   Float_t z1=h4->Integral(0,200);// changed from h3 to h4
   h4->Scale(1.0/z1);
   h4->SetLineColor(kPink-3);
   
   h4->Rebin(4);
   h4->GetXaxis()->SetRangeUser(0,80);
   h_met_rereco_all->Rebin(4);
   h_met_rereco_all->GetXaxis()->SetRangeUser(0,80);
   //h_met_rereco->Draw();
   

   
   TCanvas *Q = new TCanvas("Q","Comparison",900,900); //1200,900
   TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();

   gStyle->SetStatTextColor(kPink-3);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2); 
   h4->Draw();
   //Fakerhoreweighted_rereco->Draw();
   gStyle->SetOptStat(111111);
   


   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_met_rereco_all->Draw("sames");
   gStyle->SetOptStat(111111);
   
   
   
   

   
   
   
   
   Q->cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad2->Draw();
   pad2->cd();
   
   //Fake_rereco->SetStats(0);
   //Electron_rereco->SetStats(0);
   TH1F *h2 = new TH1F("h2","candidate over mixture",50,0,50);
   h2->Sumw2();
   
   //h2->SetStats(0);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.15);
   gStyle->SetStatH(0.15);
   h2->Divide(h_met_rereco_all, h4);
   //gStyle->SetStatTextColor(kGreen+3);
   
   h2->GetXaxis()->SetRangeUser(0,20);
   h2->GetYaxis()->SetRangeUser(0.2,2.0);
   //h2->SetMarkerStyle(21);
   h2->Draw("ep");
   
   Q->SaveAs("CandidateMixtureeeffMET.eps");

   TCanvas *R = new TCanvas("R","Ratio",600,450); //1200,900
   R->cd();
   //gStyle->SetStatTextColor(kCyan+3);
   
   h2->Draw();


   //

 
}
       
