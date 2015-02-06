#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>


void DrawAsymErrorFakeRate(){
  
  TFile *F22                        = new TFile("NotimpSymmetricPt_eeRemoved.root","READ");
  TH1F *h_met_eySample_Scaled       = (TH1F*)F22->Get("h_met_eySample_Scaled");
  TH1F *h_met_eySample_Scaled_highf = (TH1F*)F22->Get("h_met_eySample_Scaled_highf");
  TH1F *h_met_eySample_Scaled_lowf  = (TH1F*)F22->Get("h_met_eySample_Scaled_lowf");
  
  
  Float_t bins2[] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,150,300};
  Int_t  binnum2 = sizeof(bins2)/sizeof(Float_t) - 1; 
 
  TH1F* g1 = new TH1F("rerecoewkcorerrorsym","ewk error propagation, sym pt", binnum2, bins2);

  float metvalue[18]                         = {0};
  float metvalueup[18]                       = {0};
  float metvaluedown[18]                     = {0};
  float metvalueerror[18]                    = {0};
  float errorup[18]                          = {0};
  float errordown[18]                        = {0};
  float totalerror[18]                       = {0};

  for(int j=0 ; j<18 ; ++j){
    float sumdiempt             = 0;
    
    metvalue[j] =  h_met_eySample_Scaled->GetBinContent(j);
    cout << "metvalue[" << j << "]: " << metvalue[j] << endl;
    metvalueup[j] = h_met_eySample_Scaled_highf->GetBinContent(j);
    cout << "metvalueup[" << j << "]: " << metvalueup[j] << endl;
    metvaluedown[j] = h_met_eySample_Scaled_lowf->GetBinContent(j);
    cout << "metvaluedown[" << j << "]: " << metvaluedown[j] << endl;
    metvalueerror[j] =  h_met_eySample_Scaled->GetBinError(j);
   
    errorup[j]   = metvalueup[j]-metvalue[j];
    errordown[j] = metvalue[j]-metvaluedown[j]; 
    
    totalerror[j] = sqrt((errorup[j]+errordown[j])*(errorup[j]+errordown[j])+(metvalueerror[j])*(metvalueerror[j]));

    g1->SetBinContent(j, metvalue[j]);
    g1->SetBinError(j, totalerror[j]);

  }  

  g1->Draw();

  TFile* fout2 = new TFile("ErrorCorEWK_SymPt.root","RECREATE");
  fout2->cd();
  g1->Write();
  

}
