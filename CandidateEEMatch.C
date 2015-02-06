#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>
#include <TFractionFitter.h>
#include <TCanvas.h>
#include <TPad.h>


void CandidateEEMatch(){
  TFile *F2 = new TFile("NotimpSymmetricPt_eeRemoved.root","READ");
  TH1F *h_met_candidate=(TH1F*)F2->Get("h_met_candidate");
  TFile *F23 = new TFile("ErrorCorEWK_SymPt.root","READ");
  TH1F *h_met_eySample_Scaled_WithErr=(TH1F*)F23->Get("rerecoewkcorerrorsym");
  //h_met_candidate->Sumw2();
  //h_met_rereco->Rebin(4);
  
//   Float_t x1 = h_met_rereco->Integral(0,50);
//   printf("integral200:%f\n", x1);
  //h_met_rereco->Scale(1./x1);
  
  TFile *F3 = new TFile("NormalizedSymm60ErrorCorElectron.root", "READ");
  TH1F *g_electron_cor=(TH1F*)F3->Get("normalized60rerecoeecorerrorsym");// normalized to MET 40 of data
  TFile *F4 = new TFile("NormalizedSymm60ErrorCorFake.root", "READ");
  TH1F *g_fake_cor=(TH1F*)F4->Get("normalized60rerecoffcorerrorsym");
  
  TFile *F1 = new TFile("EWKMixtureSymmetricPt_eeRemoved.root", "RECREATE");
  
  Float_t bins[] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,150,300};
  Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1; 
  TH1F *h_electron_cor = new TH1F("h_electron_cor","Electron MET after reweighting",binnum, bins);
  TH1F *h_fake_cor = new TH1F("h_fake_cor","Fake MET after reweighting",binnum, bins);
  TH1F *h_mixture_cor = new TH1F("h_mixture_cor","Electron+Fake MET after purity reweighting",binnum, bins);
  TH1F *h_mixture_cor_unscaled = new TH1F("h_mixture_cor_unscaled","Electron+Fake MET after purity reweighting, not scaled to yy-ey",binnum, bins);
  Double_t Err;
  float sum   = 0, summ60 = 0, suma60 = 0, sumam60 = 0, toterr = 0;
  float X   = h_met_eySample_Scaled_WithErr->Integral(0,11);// - 
  float Y   = h_met_candidate->Integral(0,11);
  float see = Y - X  ;
  //cout << Y << "," << X << endl;
  //float see60 = h_met_rereco->Integral(15,50);
  float see60 = h_met_candidate->IntegralAndError(12,17,Err);
  const float frac1 = 0.466175 ;
  for(int i=0;i<18;++i){
    float y = g_electron_cor->GetBinContent(i);
    float z = g_fake_cor->GetBinContent(i);
    //float ey = h_met_eySample_Scaled_WithErr->GetBinContent(i);
    
    float r(0);
    
    r = frac1*y+(1-frac1)*z;
    if(i<=11)summ60 = summ60 + r;
    if(i>=12)sum = sum + r;
  }
  cout << see << endl;
  cout << summ60 << endl;
  cout << see/summ60 << endl;
  for(int i=0;i<18;++i){
    float y  = g_electron_cor->GetBinContent(i);
    float z  = g_fake_cor->GetBinContent(i);
    
    float r(0);
    r = frac1*y+(1-frac1)*z;
    h_mixture_cor_unscaled->SetBinContent(i,r);
    r = r*see/summ60;
    if(i<=11)sumam60  = sumam60 + r;
    if(i>=12)suma60 = suma60 + r;
    float y1 = g_electron_cor->GetBinError(i);
    float z1 = g_fake_cor->GetBinError(i);
    float err(0), err1(0);
    err1 = sqrt(frac1*frac1*y1*y1+(1-frac1)*(1-frac1)*z1*z1);
    err  = sqrt(frac1*frac1*y1*y1+(1-frac1)*(1-frac1)*z1*z1)*see/summ60;
    if(i>=12)toterr = toterr + err*err*err*err;
    h_electron_cor->SetBinContent(i,y);
    h_electron_cor->SetBinError(i,y1);
    h_fake_cor->SetBinContent(i,z);
    h_fake_cor->SetBinError(i,z1);


    h_mixture_cor->SetBinContent(i,r);
    h_mixture_cor->SetBinError(i,err);
    h_mixture_cor_unscaled->SetBinError(i, err1);
  }
  printf("unreweighted mixture integral after 60:%f\n",sum);
  printf("unreweighted mixture integral upto 60 :%f\n",summ60);
  printf("reweighted mixture integral upto 60   :%f\n",sumam60);
  printf("reweighted mixture integral above 60  :%f\n",suma60);
  printf("reweighted mixture integral error 60  :%f\n",sqrt(toterr));
  printf("candidate integral upto 60            :%f\n",see);
  printf("candidate integral after 60           :%f\n",see60);
  printf("candidate integral error after 60     :%f\n",Err);
  


  TCanvas *Q = new TCanvas("Q","Comparison",900,900); //1200,900
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();
   
   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   pad1->SetLogy();
   h_met_candidate->SetLineColor(kGreen+3);
   h_met_candidate->GetXaxis()->SetRangeUser(0,300);
   h_met_candidate->Draw();
   //Fakerhoreweighted_prompt->Draw();
   gStyle->SetOptStat(111111);
   
   

   gStyle->SetStatTextColor(2);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_mixture_cor->SetLineColor(2);
   h_mixture_cor->Draw("sames");
   gStyle->SetOptStat(111111);
   
   
   
   
   Q->cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad2->Draw();
   pad2->cd();
   TH1F *h2 = new TH1F("h2","candidate over mixture",binnum,bins);
   for(int i =0 ; i<18; ++i){
     /*float y = g_electron_cor->Eval(i*4);
     float ery= g_electron_cor->GetErrorY(i);
     float y2= g_fake_cor->Eval(i*4);
     float ery2= g_fake_cor->GetErrorY(i);
     if(y2!=0)h2->SetBinContent(i,y/y2);
     float erz=0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     printf("error%d:%f\n",i,erz);
     h2->SetBinError(i,erz);*/

     float y = h_met_candidate->GetBinContent(i);
     float ery= h_met_candidate->GetBinError(i);
     float y2= h_mixture_cor->GetBinContent(i);
     float ery2= h_mixture_cor->GetBinError(i);
     if(y2!=0)h2->SetBinContent(i,y/y2);
     float erz=0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     //printf("error%d:%f\n",i,erz);
     h2->SetBinError(i,erz);

   }
   
   //h2->Sumw2();
   
   //h2->SetStats(0);
   

   for(int i = 0; i< 18 ; ++i){
     float so = h2->GetBinContent(i);
     //printf("so:%f\n",so);
   }
   //gStyle->SetStatTextColor(kGreen+3);
   
   h2->GetXaxis()->SetRangeUser(0,300);
   h2->GetYaxis()->SetRangeUser(0.2,2.0);
   //h2->SetMarkerStyle(21);
   h2->Draw("ep");
   
   
   
   
   //g2->Draw("AP");
   
   Q->SaveAs("CandidateReweighted60Electron0_46Fake0_54MixtureCorPfMetReRecoEWKeeSubtractSymmetricNoStack.eps");
   Q->SaveAs("CandidateReweighted60Electron0_46Fake0_54MixtureCorPfMetReRecoEWKeeSubtractSymmetricNoStack.pdf");

   TCanvas *R = new TCanvas("R","Ratio",600,450); //1200,900
   R->cd();
   h2->Draw();
   F1->cd();
   h_mixture_cor->Write();
   h_mixture_cor_unscaled->Write();
   h_met_eySample_Scaled_WithErr->Write();
   h_met_candidate->Write();

  /*for(int i =0 ;i<50; ++i){
    Float_t y1 = g_electron_cor->Eval(i);
    Float_t y2 = g_electron_cor->GetErrorY(i);
    h_electron_cor->SetBinContent(i,y1);
    h_electron_cor->SetBinError(i,y2);
  }*/

  //h_electron_cor->Draw();
  return;
  
 

}
