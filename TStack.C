#include <TH1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TAttFill.h>
#include <TFile.h>
#include <iostream>
#include <TFractionFitter.h>
#include <TCanvas.h>
#include <TPad.h> 
#include <THStack.h>
#include <TLegend.h>


void TStack(){  

	Float_t bins[] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,150,300};
	Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
	THStack hs("hs","QCD (purity 0.466) and EWK MET with candidate histograms,Symmetric Pt");
	TFile *f = new TFile("SignalSampleSymmetric.root", "READ");
	//Float_t x1 = h_met_eySample_Scaled->Integral(0,300);
	//h_met_eySample_Scaled->Scale(1./x1);
// 	TH1F *h_met_eySample_Scaled_WithErr = (TH1F*)f->Get("rerecoewkcorerrorsym");
// 	TH1F *h_mixture_cor = (TH1F*)f->Get("h_mixture_cor");
// 	TH1F *h_met_candidate = (TH1F*)f->Get("h_met_candidate");
        TH1F *h_met_candidate_error    = (TH1F*)f->Get("h_met_candidate_error");
        TH1F *h_met_QCD_error          = (TH1F*)f->Get("h_met_QCD_error");
        TH1F *h_met_EWK_error          = (TH1F*)f->Get("h_met_EWK_error");
        TH1F *h_met_error              = new TH1F("h_met_error", "QCD+EWk shape error+stat error", binnum, bins);
        TH1F *h_met_QCD_subtract_error = new TH1F("h_met_QCD_subtract_error", "QCD+EWK-error", binnum, bins);
        TH1F *h_met_QCD_plus_EWK     = new TH1F("h_met_QCD_plus_EWK", "QCD+EWK", binnum, bins);

        for(int k=0;k<18;++k){
          float x  = h_met_QCD_error->GetBinContent(k);
          float x1 = h_met_QCD_error->GetBinError(k);
          float d  = h_met_EWK_error->GetBinContent(k);
          float d2 = h_met_EWK_error->GetBinError(k);
          float y  = x + d;
          float y1 = sqrt(x1*x1+d2*d2);
          h_met_error->SetBinContent(k, y1);
          //if((y-0.5*y1)<0)h_met_QCD_subtract_error->SetBinContent(k, 0);
          h_met_QCD_subtract_error->SetBinContent(k, x-0.5*y1);
          h_met_QCD_plus_EWK->SetBinContent(k, y);
        }
  
	h_met_EWK_error->SetFillColor(kBlue);
	hs.Add(h_met_EWK_error);

	h_met_QCD_subtract_error->SetFillColor(kRed);
        h_met_QCD_plus_EWK->SetFillColor(kRed);
        h_met_QCD_plus_EWK->SetFillStyle(3001);
	hs.Add(h_met_QCD_subtract_error);

        h_met_error->SetFillColor(kOrange);
        h_met_error->SetFillStyle(3007);
        hs.Add(h_met_error);

        //hs.Add(h_met_QCD_plus_EWK);

        //h_met_QCD_plus_EWK->SetFillColor(kWhite);
        //h_met_QCD_subtract_error->SetFillColor(kRed);
	
	h_met_candidate_error->SetMarkerStyle(kFullDotLarge);
	h_met_candidate_error->GetXaxis()->SetRangeUser(0,50);
	
	
	TCanvas c1("c1","stacked hists",1200,900);
	
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	
	
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();
	gStyle->SetOptStat(0);
	
        //h_met_QCD_plus_EWK->Draw();
	hs.Draw("hist");
        h_met_QCD_plus_EWK->Draw("hist sames");
	h_met_candidate_error->Draw("sames");
	TLegend *leg = new TLegend(0.9,0.7,0.7,0.9);
	
	leg->AddEntry(h_met_QCD_subtract_error,"QCD MET","f");
	leg->AddEntry(h_met_EWK_error,"EWK MET","f");
        leg->AddEntry(h_met_error, "combined error","f");
	leg->AddEntry(h_met_candidate_error,"candidate MET","lep");
	leg->Draw();
	
	
	
	c1.cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.28);
	
	pad2->Draw();
	pad2->cd();
	TH1F *h2 = new TH1F("h2","candidate over QCD+EWK mixture",binnum,bins);
	
	for(int i=0;i<12;++i){
	  float y = h_met_candidate_error->GetBinContent(i);
	  float erry = h_met_candidate_error->GetBinError(i);
	  float y1  = h_met_QCD_error->GetBinContent(i);
	  float erry1 = h_met_QCD_error->GetBinError(i);
	  float z1  = h_met_EWK_error->GetBinContent(i);
	  float errz1  = h_met_EWK_error->GetBinError(i);
	  float r = y1+z1;
	  float errr = sqrt(erry1*erry1+errz1*errz1);
	
	  if(r!=0)h2->SetBinContent(i,y/r);
	  float erz=0;
	  if(y!=0 && r!=0)erz= (y/r)*sqrt((erry/y)*(erry/y)+(errr/r)*(errr/r));
	
	  h2->SetBinError(i,erz);
	}
	
	h2->GetXaxis()->SetRangeUser(0,300);
	h2->GetYaxis()->SetRangeUser(0.2,2.0);
	
	h2->Draw("ep");
	
	
	c1.SaveAs("CandidateReweighted60Electron0_46Fake0_54MixtureCorPfMetReRecoEWKeeSubtractSymmetricError.eps");
	c1.SaveAs("CandidateReweighted60Electron0_46Fake0_54MixtureCorPfMetReRecoEWKeeSubtractSymmetricError.pdf");
	
	TCanvas *R = new TCanvas("R","Ratio",600,450); //1200,900
	R->cd();
	gStyle->SetOptStat(0);
	h2->Draw();
      

}
