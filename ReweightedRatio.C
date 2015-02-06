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
#include <TLegend.h>
void ReweightedRatio(){
  TFile *F2 = new TFile("DiffFake_eeRemoved_SymmetricPt.root","READ");
  TH1F *h_diEMPt_candidate=(TH1F*)F2->Get("h_diEMPt_candidate");
  TH1F *h_rho_candidate=(TH1F*)F2->Get("h_rho_candidate");
  TH1F *h_diEMPt_eeSample=(TH1F*)F2->Get("h_diEMPt_eeSample");
  TH1F *h_rho_eeSample_diEMPtReweighted=(TH1F*)F2->Get("h_rho_eeSample_diEMPtReweighted");
  TH1F *h_met_eeSample=(TH1F*)F2->Get("h_met_eeSample");
  TH1F *h_met_ffSample=(TH1F*)F2->Get("h_met_ffSample");
  TH1F *h_rho_ffSample=(TH1F*)F2->Get("h_rho_ffSample");
  TH1F *h_diEMPt_ffSample=(TH1F*)F2->Get("h_diEMPt_ffSample");
  
  h_diEMPt_candidate->Sumw2();
  h_rho_candidate->Sumw2();
  h_diEMPt_eeSample->Sumw2();
  h_rho_eeSample_diEMPtReweighted->Sumw2();
  h_rho_ffSample->Sumw2();
  //h_met_eeSample->Sumw2();
  //h_met_ffSample->Sumw2();
  h_diEMPt_ffSample->Sumw2();
  
  
  float x1=h_diEMPt_candidate->Integral();
  h_diEMPt_candidate->Scale(1./x1);
  h_diEMPt_candidate->SetLineColor(kGreen+3);
  float x2=h_rho_candidate->Integral();
  h_rho_candidate->Scale(1./x2);
  h_rho_candidate->SetLineColor(kGreen+3);
  float x3=h_diEMPt_eeSample->Integral();
  h_diEMPt_eeSample->Scale(1./x3);
  h_diEMPt_eeSample->SetLineColor(kBlue);
  float x4=h_rho_eeSample_diEMPtReweighted->Integral();
  h_rho_eeSample_diEMPtReweighted->Scale(1./x4);
  h_rho_eeSample_diEMPtReweighted->SetLineColor(kBlue);
  float x5=h_rho_ffSample->Integral();
  h_rho_ffSample->Scale(1./x5);
  h_rho_ffSample->SetLineColor(kRed);
  float x6=h_met_eeSample->Integral();
  h_met_eeSample->Scale(1./x6);
  h_met_eeSample->SetLineColor(kBlue);
  float x7=h_met_ffSample->Integral();
  h_met_ffSample->Scale(1./x7);
  h_met_ffSample->SetLineColor(kRed);
  float x8=h_diEMPt_ffSample->Integral();
  h_diEMPt_ffSample->Scale(1./x8);
  h_diEMPt_ffSample->SetLineColor(kRed);


  TCanvas *Q = new TCanvas("Q","Comparison",900,900); //1200,900
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad1->Draw();
   pad1->cd();
   
   gStyle->SetStatTextColor(kBlue);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   //pad1->SetLogy();
   h_diEMPt_eeSample->GetXaxis()->SetRangeUser(0,50);
   h_diEMPt_eeSample->Draw();
   gStyle->SetOptStat(111111);
   
   
   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_diEMPt_candidate->Draw("sames");
   gStyle->SetOptStat(111111);
   
   Q->cd();
   TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad2->Draw();
   pad2->cd();  
   TH1F *h2 = new TH1F("h2","candidate over ee",50,0,50);
   for(int i =0 ; i<50; ++i){
     float y = h_diEMPt_candidate->GetBinContent(i);
     float ery= h_diEMPt_candidate->GetBinError(i);
     float y2= h_diEMPt_eeSample->GetBinContent(i);
     float ery2= h_diEMPt_eeSample->GetBinError(i);
     float erz = 0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     
     h2->SetBinError(i,erz);
     if(y2!=0)h2->SetBinContent(i,y/y2);
   }
   h2->GetYaxis()->SetRangeUser(0,4);
   h2->Draw();
   Q->SaveAs("DiEMPtRatio_eeSymmetric.eps");
   
   
   
   
   TCanvas *R = new TCanvas("R","Comparison",900,900); //1200,900
   TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad3->Draw();
   pad3->cd();
   
   gStyle->SetStatTextColor(kBlue);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   //pad1->SetLogy();
   h_rho_eeSample_diEMPtReweighted->GetXaxis()->SetRangeUser(0,50);
   h_rho_eeSample_diEMPtReweighted->Draw();
   gStyle->SetOptStat(111111);
   
   
   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_rho_candidate->Draw("sames");
   gStyle->SetOptStat(111111);
   
   R->cd();
   TPad *pad4 = new TPad("pad4","pad4",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad4->Draw();
   pad4->cd();  
   TH1F *h3 = new TH1F("h3","candidate over ee",50,0,50);
   for(int i =0 ; i<50; ++i){
     float y = h_rho_candidate->GetBinContent(i);
     float ery= h_rho_candidate->GetBinError(i);
     float y2= h_rho_eeSample_diEMPtReweighted->GetBinContent(i);
     float ery2= h_rho_eeSample_diEMPtReweighted->GetBinError(i);
     float erz = 0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     
     h3->SetBinError(i,erz);
     if(y2!=0)h3->SetBinContent(i,y/y2);
   }
   h3->GetYaxis()->SetRangeUser(0,4);
   h3->Draw();
   R->SaveAs("RhoRatio_eeSymmetric.eps");
   
   
   
   TCanvas *S = new TCanvas("S","Comparison",900,900); //1200,900
   TPad *pad5 = new TPad("pad5","pad5",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad5->Draw();
   pad5->cd();
   
   gStyle->SetStatTextColor(kRed);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   //pad1->SetLogy();
   h_rho_ffSample->GetXaxis()->SetRangeUser(0,50);
   h_rho_ffSample->Draw();
   gStyle->SetOptStat(111111);
   
   
   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_rho_candidate->Draw("sames");
   gStyle->SetOptStat(111111);
   
   S->cd();
   TPad *pad6 = new TPad("pad6","pad6",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad6->Draw();
   pad6->cd();  
   TH1F *h4 = new TH1F("h4","candidate over ff",50,0,50);
   for(int i =0 ; i<50; ++i){
     float y = h_rho_candidate->GetBinContent(i);
     float ery= h_rho_candidate->GetBinError(i);
     float y2= h_rho_ffSample->GetBinContent(i);
     float ery2= h_rho_ffSample->GetBinError(i);
     float erz = 0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     
     h4->SetBinError(i,erz);
     if(y2!=0)h4->SetBinContent(i,y/y2);
     //if(y2!=0)cout << y/y2 << endl;
   }
   h4->GetYaxis()->SetRangeUser(0,4);
   h4->Draw();
   S->SaveAs("RhoRatio_ffSymmetric.eps");
   
   
   
   TCanvas *T = new TCanvas("T","Comparison",900,900); //1200,900
   TPad *pad7 = new TPad("pad7","pad7",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad7->Draw();
   pad7->cd();
   
   gStyle->SetStatTextColor(kRed);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   //pad1->SetLogy();
   h_met_ffSample->GetXaxis()->SetRangeUser(0,50);
   h_met_ffSample->Draw();
   gStyle->SetOptStat(111111);
   
   
   gStyle->SetStatTextColor(kBlue);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_met_eeSample->Draw("sames");
   gStyle->SetOptStat(111111);
   
   T->cd();
   TPad *pad8 = new TPad("pad8","pad8",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad8->Draw();
   pad8->SetLogy();
   pad8->cd(); 
   TH1F *h5 = new TH1F("h5","ee over ff",50, 0,50);
   for(int i =0 ; i<18; ++i){
     float y = h_met_eeSample->GetBinContent(i);
     float ery= h_met_eeSample->GetBinError(i);
     float y2= h_met_ffSample->GetBinContent(i);
     float ery2= h_met_ffSample->GetBinError(i);
     float erz = 0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     
     h5->SetBinError(i,erz);
     if(y2!=0)h5->SetBinContent(i,y/y2);
     //if(y2!=0)cout << y/y2 << endl;
   }
   h5->GetYaxis()->SetRangeUser(0,4);
   h5->Draw();
   T->SaveAs("METRatio_eeffSymmetric.eps");
   TCanvas *U = new TCanvas("U","Comparison",450,450); //1200,900  
   h5->Draw();



   /*TCanvas *V = new TCanvas("V","Comparison",900,900); //1200,900
   TPad *pad9 = new TPad("pad9","pad9",0,0.3,1,1);
   
   //pad1->SetBottomMargin(0);
   pad9->Draw();
   pad9->cd();
   
   gStyle->SetStatTextColor(kRed);
   gStyle->SetStatY(0.9);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   //pad1->SetLogy();
   h_diEMPt_ffSample->GetXaxis()->SetRangeUser(0,50);
   h_diEMPt_ffSample->Draw();
   gStyle->SetOptStat(111111);
   
   
   gStyle->SetStatTextColor(kGreen+3);
   gStyle->SetStatY(0.5);
   gStyle->SetStatX(0.9);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.2);
   h_diEMPt_candidate->Draw("sames");
   gStyle->SetOptStat(111111);
   
   V->cd();
   TPad *pad10 = new TPad("pad10","pad10",0,0,1,0.28);
   //pad2->SetTopMargin(0);
   //pad2->SetBottomMargin(0);
   pad10->Draw();
   pad10->cd();  
   TH1F *h6 = new TH1F("h6","candidate over ff",50,0,50);
   for(int i =0 ; i<50; ++i){
     float y = h_diEMPt_candidate->GetBinContent(i);
     float ery= h_diEMPt_candidate->GetBinError(i);
     float y2= h_diEMPt_ffSample->GetBinContent(i);
     float ery2= h_diEMPt_ffSample->GetBinError(i);
     float erz = 0;
     if(y!=0 && y2!=0)erz= (y/y2)*sqrt((ery/y)*(ery/y)+(ery2/y2)*(ery2/y2));
     
     h6->SetBinError(i,erz);
     if(y2!=0)h6->SetBinContent(i,y/y2);
   }
   h6->GetYaxis()->SetRangeUser(0,4);
   h6->Draw();
   V->SaveAs("DiEMPtRatio_ffSymmetric.eps");*/
   
   
   
   return;
  
}

