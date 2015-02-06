{

  TFile *F2 = new TFile("ImportantReReco.root","READ");
  TH1F *DiEMPtRatioWeightedee[1000]; //= (TH1F*)Fnew->Get("candidateeediEMPt");
  TH1F *RhoRatioWeightedee[1000];
  TH1F *RhoRatioWeightedff[1000];
 // TH1F *h_met_ffSample = (TH1F*)F2->Get("h_met_ffSample");
  

  Float_t x1= h_met_ffSample->Integral();
  h_met_ffSample->Scale(1./x1);
  h_met_ffSample->SetLineColor(2);
  h_met_ffSample->GetXaxis()->SetRangeUser(0,80);

  for(int g=0; g<1000; ++g){
    char *nameGraph        = new char[35];
    char *nameGraph2        = new char[35];
    char *nameGraph3        = new char[35];
    sprintf(nameGraph,"DiEMPtReweightedee%d",g+1);
    DiEMPtRatioWeightedee[g] = (TH1F*)F2->Get(nameGraph);
    sprintf(nameGraph2,"RhoReweightedee%d",g+1);
    RhoRatioWeightedee[g] = (TH1F*)F2->Get(nameGraph2);
    sprintf(nameGraph3,"RhoReweightedff%d",g+1);
    RhoRatioWeightedff[g] = (TH1F*)F2->Get(nameGraph3);
    /*Float_t x2 = DiEMPtRatioWeightedee[g]->Integral(0,200);
    DiEMPtRatioWeightedee[g]->Scale(1./x2);
    DiEMPtRatioWeightedee[g]->SetLineColor(4);*/
    DiEMPtRatioWeightedee[g]->GetXaxis()->SetRangeUser(0,80);
    RhoRatioWeightedee[g]->GetXaxis()->SetRangeUser(0,80);
    RhoRatioWeightedff[g]->GetXaxis()->SetRangeUser(0,80);
  }

  
  TCanvas *Q = new TCanvas("Q","Comparison",600,450);
  Q->cd();

  /*gStyle->SetStatTextColor(4);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  h_met_ffSample->Draw();
  gStyle->SetOptStat(111111);*/


  gStyle->SetStatTextColor(2);
  gStyle->SetStatY(0.5);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  RhoRatioWeightedee[0]->Draw();
  
  for(int j=1;j<1000;++j){
    if(j%20==0)printf("processed:%d\n",j);

    
    RhoRatioWeightedee[j]->Draw("sames");
  }


  Q->SaveAs("RhoErrorPropagatedeeMETReRecoCor.eps");
  Q->SaveAs("RhoErrorPropagatedeeMETReRecoCor.pdf");


  TCanvas *R = new TCanvas("R","Comparison",600,450);
  R->cd();

  /*gStyle->SetStatTextColor(4);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  h_met_ffSample->Draw();
  gStyle->SetOptStat(111111);*/


  gStyle->SetStatTextColor(2);
  gStyle->SetStatY(0.5);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.15);
  RhoRatioWeightedff[0]->Draw();
  
  for(int j=1;j<1000;++j){
    if(j%20==0)printf("processed:%d\n",j);

    
    RhoRatioWeightedff[j]->Draw("sames");
  }


  R->SaveAs("ErrorPropagatedffMETReRecoCor.eps");
  R->SaveAs("ErrorPropagatedffMETReRecoCor.pdf");
  



}
