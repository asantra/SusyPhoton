
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFile.h>
#include <TStyle.h>
#include <iostream>


void QCDEWKError(){

  TFile *F22 = new TFile("EWKMixture_eeRemoved.root","READ");
  TH1F *h_met_candidate = (TH1F*)F22->Get("h_met_candidate");
  TH1F *h_met_eySample_Scaled_WithErr = (TH1F*)F22->Get("rerecoewkcorerror");
  TH1F *h_mixture_cor = (TH1F*)F22->Get("h_mixture_cor");
  TH1F *h_mixture_cor_unscaled = (TH1F*)F22->Get("h_mixture_cor_unscaled");
  
  TFile *F23 = new TFile("SignalSample.root", "RECREATE");
  F23->cd();
  Float_t bins2[] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,150,300};
  Int_t  binnum2 = sizeof(bins2)/sizeof(Float_t) - 1;
  TH1F *h_met_signal          = new TH1F("h_met_signal", "signal after MET>=60 GeV", binnum2, bins2); 
  TH1F *h_met_candidate_error = new TH1F("h_met_candidate_error", "candidate sample after proper error considering EWK background", binnum2, bins2);
  TH1F *h_met_QCD_error       = new TH1F("h_met_QCD_error", "QCD background with shape error+stat error", binnum2, bins2);
  TH1F *h_met_EWK_error       = new TH1F("h_met_EWK_error", "EWK background with fake rate error+stat error", binnum2, bins2);

  float metvalue[18]       = {0};
  float metvalueerror[18]  = {0};
  float candidateerror[18] = {0};
  float backgrounderror[18]= {0};
  float qcderror[18]       = {0};
  float qcderroroverall[18]= {0};
  float ewkerror[18]       = {0};
  float candidate60(0), candidateerror60(0), ewk60(0), ewkerror60(0);
  float QCD60(0), QCDError60(0);

  for(int j=0 ; j<=11 ; ++j){
    candidate60      = candidate60 + h_met_candidate->GetBinContent(j);
    candidateerror60 = candidateerror60 + h_met_candidate->GetBinError(j)*h_met_candidate->GetBinError(j);
    ewk60            = ewk60 + h_met_eySample_Scaled_WithErr->GetBinContent(j);
    ewkerror60       = ewkerror60 + h_met_eySample_Scaled_WithErr->GetBinError(j)*h_met_eySample_Scaled_WithErr->GetBinError(j);
    QCD60            = QCD60 + h_mixture_cor_unscaled->GetBinContent(j);
    QCDError60       = QCDError60 + h_mixture_cor_unscaled->GetBinError(j)*h_mixture_cor_unscaled->GetBinError(j);
  }
  //cout << "before sqrt: " << candidateerror60 << endl;
  candidateerror60 = sqrt(candidateerror60); // needed for QCD scale error
  //cout << "after sqrt: " << candidateerror60 << endl;
  ewkerror60       = sqrt(ewkerror60); // needed for QCD scale error
  QCDError60       = sqrt(QCDError60); // needed for QCD scale error
  float fracQCD = sqrt((candidateerror60*candidateerror60+ewkerror60*ewkerror60)/((candidate60-ewk60)*(candidate60-ewk60))+(QCDError60/QCD60)*(QCDError60/QCD60)); // QCD scale error
  //cout << "fracQCD: " << fracQCD << endl;

  for(int k=0; k<18; ++k){
    //cout << "h_mixture_cor_content[" << k << "]: " << h_mixture_cor->GetBinContent(k) << endl;
    //cout << "h_mixture_cor_unscaled_content[" << k << "]: " << h_mixture_cor_unscaled->GetBinContent(k) << endl;
    qcderror[k] = h_mixture_cor->GetBinContent(k)*sqrt(fracQCD*fracQCD + (h_mixture_cor_unscaled->GetBinError(k)*h_mixture_cor_unscaled->GetBinError(k)/(h_mixture_cor_unscaled->GetBinContent(k)*h_mixture_cor_unscaled->GetBinContent(k)))); // total QCD normalization error
    qcderroroverall[k] = sqrt(h_mixture_cor->GetBinError(k)*h_mixture_cor->GetBinError(k)/(h_mixture_cor->GetBinContent(k)*h_mixture_cor->GetBinContent(k))+qcderror[k]*qcderror[k]); // combination of shape and normalization error for QCD
    //cout << "qcderror[" << k <<"]: " << qcderror[k] << endl;
    //cout << "qcderroroverall[" << k <<"]: " << qcderroroverall[k] << endl;
    ewkerror[k] = h_met_eySample_Scaled_WithErr->GetBinError(k);
//     if(k<=11)candidateerror[k] = sqrt(h_met_candidate->GetBinContent(k) - h_met_eySample_Scaled_WithErr->GetBinContent(k));
    candidateerror[k] = h_met_candidate->GetBinError(k);
    //cout << "ewkerror[" << k << "]: "<< ewkerror[k] << endl;
    //cout << "candidate error[" << k << "]: " << candidateerror[k] << endl;
    backgrounderror[k] = sqrt(qcderroroverall[k]*qcderroroverall[k] + ewkerror[k]*ewkerror[k]);
    h_met_candidate_error->SetBinContent(k, h_met_candidate->GetBinContent(k)); // filling candidate
    h_met_candidate_error->SetBinError(k, candidateerror[k]); // filling candidate error
    h_met_QCD_error->SetBinContent(k, h_mixture_cor->GetBinContent(k)); // filling QCD
    h_met_QCD_error->SetBinError(k, qcderroroverall[k]); //filling QCD error
    h_met_EWK_error->SetBinContent(k, h_met_eySample_Scaled_WithErr->GetBinContent(k)); // filling EWK
    h_met_EWK_error->SetBinError(k, ewkerror[k]);  // filling EWK error
  }

  for(int j=12 ; j<18 ; ++j){
    metvalue[j] = h_met_candidate->GetBinContent(j) - ( h_mixture_cor->GetBinContent(j) + h_met_eySample_Scaled_WithErr->GetBinContent(j));
    candidateerror[j] = h_met_candidate->GetBinError(j);
    metvalueerror[j] = sqrt(candidateerror[j]*candidateerror[j]+backgrounderror[j]*backgrounderror[j]);
    h_met_signal->SetBinContent(j,metvalue[j]);
    h_met_signal->SetBinError(j,metvalueerror[j]);
    //cout << "value[" << j <<"] : " << metvalue[j] << endl;
    //cout << "valueerror[" << j <<"] : " << metvalueerror[j] << endl;
  }
  
  TCanvas *R = new TCanvas("R","QCD, EWK candidate",1200,900);
  R->Divide(2,2);
  R->cd(1);
  R->SetLogy();
  h_met_candidate_error->Draw();
  R->cd(2);
  R->SetLogy();
  h_met_QCD_error->Draw();
  R->cd(3);
  R->SetLogy();
  h_met_EWK_error->Draw();
  R->cd(4);
  R->SetLogy();
  h_met_signal->Draw();
  F23->Write();
  //h_met_signal->Write();
  

}
