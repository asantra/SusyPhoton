#include <TH1.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>


void DrawAsymError(){

  TFile *F22 = new TFile("NotimpSymmetricPt_eeRemoved.root","READ");
  TH1F *h_met_candidate = (TH1F*)F22->Get("h_met_candidate");
  TH1F *h_met_eySample_Scaled = (TH1F*)F22->Get("h_met_eySample_Scaled");
  
  Float_t bins2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
  Int_t  binnum2 = sizeof(bins2)/sizeof(Float_t) - 1; 
 
  TH1F* g1 = new TH1F("normalized60rerecoffcorerrorsym","ff error propagation", binnum2, bins2);
  float integral_candidate = h_met_candidate->Integral(0,11);
  //float integral_ey = h_met_eySample_Scaled->Integral(0,11);
  float integral60 = integral_candidate;// - integral_ey;
  printf("integral60: %f\n",integral60);

  //TFile *F2 = new TFile("PromptRecoFakeCorPfMet.root","READ");
  
  TH1F *RhoRatioWeightedff[1000];
  //TH1F *RhoRatioWeightedffPrompt[1000];
  TH1F *h_met_ffSample_RhoReweighted = (TH1F*)F22->Get("h_met_ffSample_RhoReweighted");
  

  for(int g=0; g<1000; ++g){
    char *nameGraph        = new char[50];
    //char *nameGraph2        = new char[50];
    sprintf(nameGraph,"RhoReweightedff%d",g+1);
    RhoRatioWeightedff[g] = (TH1F*)F22->Get(nameGraph);
    //RhoRatioWeightedff[g]->Rebin(4);
    /*Float_t y1 = RhoRatioWeightedff[g]->Integral(0,50);
    RhoRatioWeightedff[g]->Scale(1./y1);*/

//     sprintf(nameGraph2,"RhoReweightedffPrompt%d",g+1);
//     RhoRatioWeightedffPrompt[g] = (TH1F*)F2->Get(nameGraph2);
//     RhoRatioWeightedffPrompt[g]->Rebin(4);
    
  }

  


  float xbin[18][1000]                 = {{0},{0}};
  //float xbinprompt[200][1000]           = {{0},{0}};
  Int_t N                              = 18; 
  float Bin[18]                        = {0};
  float xerror1[18]                    = {0};
  float xerror2[18]                    = {0};
  float metvalue[18]                   = {0};
  //float metvalueprompt[18]             = {0};
  float normalizedmetvalue[18]         = {0};
  
  //float normalizedmetvalueprompt[200]   = {0};
  float metvalueerror[18]              = {0};
  //float metvalueerrorprompt[200]        = {0};
  float normalizedmetvalueerror[18]    = {0};
  //float normalizedmetvalueerrorprompt[200]= {0};
  float totalerror[18]                 = {0};
  //float totalerrorprompt[200]           = {0};
  float normalizedtotalerror[18]       = {0};
  //float normalizedtotalerrorprompt[200] = {0}; 

  float central[18]                    = {0};
  //float centralprompt[200]              = {0};
  float errorup[18]                    = {0};
  //float errorupprompt[200]              = {0};
  float errordown[18]                  = {0};
  //float errordownprompt[200]            = {0};
  float normalizederrorup[18]          = {0};
  //float normalizederrorupprompt[200]    = {0};
  float normalizederrordown[18]        = {0};
  //float normalizederrordownprompt[200]  = {0};
  float integral                        = 0;
  //float integralprompt                  = 0;
  
  

  for(int j=0 ; j<18 ; ++j){
    metvalue[j] = h_met_ffSample_RhoReweighted->GetBinContent(j);
    metvalueerror[j] = h_met_ffSample_RhoReweighted->GetBinError(j);
    
    if(j<=11)integral = integral + metvalue[j] ;// getting value upto 60
    //metvalueprompt[j] = h_met_prompt_diEMPtrhoreweighted->GetBinContent(j);
    //metvalueerrorprompt[j] = h_met_prompt_diEMPtrhoreweighted->GetBinError(j);
    //integralprompt = integralprompt + metvalueprompt[j] ;

    Bin[j] = j; 
    float min(10000), max(0);// minprompt(10000), maxprompt(0);
    float sum = 0;
    //float sumprompt = 0;
    for(int k=0 ; k<1000;++k){
      xbin[j][k]       = RhoRatioWeightedff[k]->GetBinContent(j);
      sum = sum + xbin[j][k];
      if(max < xbin[j][k])max=xbin[j][k];
      if(min > xbin[j][k])min=xbin[j][k];
  
//       xbinprompt[j][k] = RhoRatioWeightedffPrompt[k]->GetBinContent(j);
//       sumprompt = sumprompt+xbinprompt[j][k];
//       if(maxprompt < xbinprompt[j][k])maxprompt=xbinprompt[j][k];
//       if(minprompt > xbinprompt[j][k])minprompt=xbinprompt[j][k];
      
    }
    
    central[j]   = sum/1000;
    //printf("central[%d]: %f\n",j,central[j]);
    errorup[j]   = 0.68*(max-central[j]);
    //printf("errorup[%d]: %f \n",j,errorup[j]);
    errordown[j] = 0.68*(central[j]-min); 
    //printf("errordown[%d]: %f \n",j,errordown[j]);


//     centralprompt[j]   =  sumprompt/1000;
//     errorupprompt[j]   = (maxprompt-centralprompt[j])*0.68;
//     errordownprompt[j] = (centralprompt[j]-minprompt)*0.68; 

    
    
  }

  printf("intergal: %f\n",integral);
  for(int j=0 ; j<18 ; ++j){
    normalizedmetvalue[j] = metvalue[j]*integral60/integral; // normalized to MET 60 of the data
    normalizedmetvalueerror[j] = metvalueerror[j]*integral60/integral;
    //printf("normalizedmetvalue[%d]:%f \n", j,normalizedmetvalue[j]);
    normalizederrorup[j]  = errorup[j]*integral60/integral;
    //printf("normalizederrorup[%d]:%f \n", j,normalizederrorup[j]);
    normalizederrordown[j]= errordown[j]*integral60/integral;
    //printf("normalizederrordown[%d]:%f \n", j,normalizederrordown[j]);

//     normalizedmetvalueprompt[j] = metvalueprompt[j]/integralprompt;
//     normalizedmetvalueerrorprompt[j] = metvalueerrorprompt[j]/integralprompt;
//     normalizederrorupprompt[j]  = errorupprompt[j]/integralprompt;
//     normalizederrordownprompt[j]= errordownprompt[j]/integralprompt;

    totalerror[j] = sqrt((metvalueerror[j]*metvalueerror[j])+(errorup[j]+errordown[j])*(errorup[j]+errordown[j]));

    //totalerrorprompt[j] = sqrt((metvalueerrorprompt[j]*metvalueerrorprompt[j])+(errorupprompt[j]+errordownprompt[j])*(errorupprompt[j]+errordownprompt[j]));


    normalizedtotalerror[j] = sqrt((normalizedmetvalueerror[j]*normalizedmetvalueerror[j])+(normalizederrorup[j]+normalizederrordown[j])*(normalizederrorup[j]+normalizederrordown[j]));

    //normalizedtotalerrorprompt[j] = sqrt((normalizedmetvalueerrorprompt[j]*normalizedmetvalueerrorprompt[j])+(normalizederrorupprompt[j]+normalizederrordownprompt[j])*(normalizederrorupprompt[j]+normalizederrordownprompt[j]));
    g1->SetBinContent(j, normalizedmetvalue[j]);
    g1->SetBinError(j, normalizedtotalerror[j]);
  }
  
  
  //TGraphErrors *g1;// *g2;

  //g1 = new TGraphErrors(N, Bin, normalizedmetvalue, xerror1, normalizedtotalerror);
  //g2 = new TGraphErrors(N, Bin, metvalueprompt, xerror1, totalerrorprompt);

  //float q = g1->Integral(0,50);
  //g1->Scale(1./q);
  //printf("q: %f",q);
  //g1->SetMarkerStyle(kFullDotSmall);
  //g1->SetMarkerColor(kRed);
  //g2->SetMarkerStyle(kFullDotSmall);
  //g2->SetMarkerColor(kBlue);
  //g1->GetXaxis()->SetLimits(0, 80);
  //g2->GetXaxis()->SetLimits(0, 80);

  g1->Draw();
  //g2->Draw();
  TFile* fout = new TFile("NormalizedSymm60ErrorCorFake.root","RECREATE");
  fout->cd();
  g1->Write();
  //g2->Write("unnormalizedpromptffrawerror");

  
  
  

}
