// ----------------------------------------------------------------------------------- //
/*
  ROOT macro for illustrating error propagation using random Gaussian numbers.

  The example is based on first doing the error propagation analytically, and then
  verifying it by running a so-called Monte-Carlo (MC) program.

  For more information on error propagation, see:
    R. J. Barlow: page 48-61 
    P. R. Bevington: page 36-48

  Author: Troels C. Petersen (NBI)
  Email:  petersen@nbi.dk
  Date:   8th of September 2011
*/
// ----------------------------------------------------------------------------------- //

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraphErrors.h>


//---------------------------------------------------------------------------------- 
void ErrorPropagation()
//---------------------------------------------------------------------------------- 
{
  


  // Set parameters:
  TRandom3 r;
  TRandom3 r1;
  TRandom3 r2;
  TRandom3 r3;
  
  
  

  // Make histograms:

  TFile *Fnew = new TFile("ElectronRatioSymmetric.root","READ");
  TH1F *ElectronDiEMPt = (TH1F*)Fnew->Get("ElectronDiEMPt");
  TH1F *ElectrondiEMPtReweightedRho = (TH1F*)Fnew->Get("ElectrondiEMPtReweightedRho");
  TH1F *ElectronDiEMPtPrompt = (TH1F*)Fnew->Get("ElectronDiEMPtPrompt");
  TH1F *ElectrondiEMPtReweightedRhoPrompt=(TH1F*)Fnew->Get("ElectrondiEMPtReweightedRhoPrompt");
  
  TFile* fout = new TFile("RandomRatioEE.root","RECREATE");
  
  
  TH1F *RatioElectronDiEMPt[1000];// for ee diEMPt
  TH1F *RatioElectrondiEMPtReweightedRho[1000]; // for ee rho
  TH1F *RatioElectronDiEMPtPrompt[1000]; //  = new TH1F( "Hist_y",  "Hist_y",  200,  0.0, 200.0);
  TH1F *RatioElectrondiEMPtReweightedRhoPrompt[1000];

  char *histnameElectronDiEMPt                   = new char[50];
  char *histnameElectrondiEMPtReweightedRho      = new char[50];
  char *histnameElectronDiEMPtPrompt                  = new char[50];
  char *histnameElectrondiEMPtReweightedRhoPrompt= new char[50];
  char *histtitleElectronDiEMPt                  = new char[50];
  char *histtitleElectrondiEMPtReweightedRho     = new char[50];
  char *histtitleElectronDiEMPtPrompt            = new char[50];
  char *histtitleElectrondiEMPtReweightedRhoPrompt= new char[50];



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------



  

  double mu1[200]          =  {0};
  double sig1[200]         =  {0};

  double mu2[200]          =  {0};
  double sig2[200]         =  {0};

  double mu3[40]          =  {0};
  double sig3[40]         =  {0};

  double mu4[40]          =  {0};
  double sig4[40]         =  {0};
  

  float rho12 =  0.0;

  for(int p =0 ; p < 200; ++p){
     mu1[p]  = ElectronDiEMPt->GetBinContent(p);
     sig1[p] = ElectronDiEMPt->GetBinError(p);  
     mu2[p]  = ElectronDiEMPtPrompt->GetBinContent(p);
     sig2[p] = ElectronDiEMPtPrompt->GetBinError(p);
  }  

  for(int q=0 ; q < 40;++q){
    mu3[q]  = ElectrondiEMPtReweightedRho->GetBinContent(q);
    sig3[q] = ElectrondiEMPtReweightedRho->GetBinError(q);
    mu4[q]  = ElectrondiEMPtReweightedRhoPrompt->GetBinContent(q); 
    sig4[q] = ElectrondiEMPtReweightedRhoPrompt->GetBinError(q);  
  }


  if ((rho12 < -1.0) || (rho12 > 1.0)) {
    printf("  ERROR: Correlation factor not in interval [-1,1], as it is %6.2f \n", rho12);
    return;
  }



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------

  for(int k=0; k < 1000; ++k){
    if(k%50 == 0)printf("processed diEMPt:%d \n", k);
    sprintf(histnameElectronDiEMPt, "ElectronDiEMPtRandom%d",k+1);
    sprintf(histtitleElectronDiEMPt,"random DiEMPt ratio for ee %d",k+1);
    RatioElectronDiEMPt[k]=new TH1F(histnameElectronDiEMPt,histtitleElectronDiEMPt,200,0,200);

    sprintf(histnameElectronDiEMPtPrompt, "ElectronDiEMPtPromptRandom%d",k+1);
    sprintf(histtitleElectronDiEMPtPrompt,"random prompt DiEMPt ratio for ee %d",k+1);
    RatioElectronDiEMPtPrompt[k]=new TH1F(histnameElectronDiEMPtPrompt,histtitleElectronDiEMPtPrompt,200,0,200);
    
    for(int j=0; j < 200; ++j){
      float u = r.Gaus(mu1[j], sig1[j]);
      RatioElectronDiEMPt[k]->SetBinContent(j,u);
      float v = r1.Gaus(mu2[j], sig2[j]);
      RatioElectronDiEMPtPrompt[k]->SetBinContent(j,v);
    }
    
    RatioElectronDiEMPt[k]->Write();
    RatioElectronDiEMPtPrompt[k]->Write();

  }



  for(int k=0; k < 1000; ++k){
    if(k%50 == 0)printf("processed rho:%d \n", k);//
    sprintf(histnameElectrondiEMPtReweightedRho, "ElectronRhoRandom%d",k+1);
    sprintf(histtitleElectrondiEMPtReweightedRho,"random rho ratio for ee rereco %d",k+1);
    RatioElectrondiEMPtReweightedRho[k]=new TH1F(histnameElectrondiEMPtReweightedRho, histtitleElectrondiEMPtReweightedRho,40,0.0,40.0);

    sprintf(histnameElectrondiEMPtReweightedRhoPrompt, "ElectronRhoRandomPrompt%d",k+1);
    sprintf(histtitleElectrondiEMPtReweightedRhoPrompt,"random rho ratio for ee prompt %d",k+1);
    RatioElectrondiEMPtReweightedRhoPrompt[k]=new TH1F(histnameElectrondiEMPtReweightedRhoPrompt,histtitleElectrondiEMPtReweightedRhoPrompt,40,0.0,40.0);
    



    for(int j=0; j < 40; ++j){  
      float w = r2.Gaus(mu3[j], sig3[j]);
      RatioElectrondiEMPtReweightedRho[k]->SetBinContent(j,w);
      float f = r3.Gaus(mu4[j], sig4[j]);
      RatioElectrondiEMPtReweightedRhoPrompt[k]->SetBinContent(j,f);
    } 
    RatioElectrondiEMPtReweightedRho[k]->Write();
    RatioElectrondiEMPtReweightedRhoPrompt[k]->Write();
  }

  
  
}


