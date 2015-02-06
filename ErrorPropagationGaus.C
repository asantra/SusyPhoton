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
void ErrorPropagationGaus()
//---------------------------------------------------------------------------------- 
{
  


  // Set parameters:
  TRandom3 r;
  TRandom3 r1;
  TRandom3 r2;
  
  
  

  // Make histograms:

  TFile *Fnew = new TFile("ElectronRatioSymmetric.root","READ");
  TH1F *candidateeediEMPt = (TH1F*)Fnew->Get("candidateeediEMPt");
  TH1F *candidateeediEMPtweightedrho = (TH1F*)Fnew->Get("candidateeediEMPtweightedrho");
  TH1F *candidateffrho = (TH1F*)Fnew->Get("candidateffrho");
  
  TFile* fout = new TFile("RandomRatioSymmetric.root","RECREATE");
  
  
  TH1F *RatioGaus[1000];// for ee diEMPt
  TH1F *RatioeeDiEMPtRatioWeightedRho[1000]; // for ee rho
  TH1F *RatioffRho[1000]; // ff rho  = new TH1F( "Hist_y",  "Hist_y",  200,  0.0, 200.0);


  char *histnameY        = new char[25];
  char *histnameRatio    = new char[25];
  char *histnameX        = new char[25];
  char *histtitleY       = new char[50];
  char *histtitleRatio   = new char[50];
  char *histtitleX       = new char[50];



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------



  

  double mu1[200]          =  {0};
  double sig1[200]         =  {0};

  double mu2[40]          =  {0};
  double sig2[40]         =  {0};

  double mu3[40]          =  {0};
  double sig3[40]         =  {0};
  

  float rho12 =  0.0;

  for(int p =0 ; p < 200; ++p){
     mu1[p]  = candidateeediEMPt->GetBinContent(p);
     sig1[p] = candidateeediEMPt->GetBinError(p);  
  }  

  for(int q=0 ; q < 40;++q){
    mu2[q]  = candidateeediEMPtweightedrho->GetBinContent(q);
    sig2[q] = candidateeediEMPtweightedrho->GetBinError(q);
    mu3[q]  = candidateffrho->GetBinContent(q); 
    sig3[q] = candidateffrho->GetBinError(q);  
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
    sprintf(histnameRatio, "DiEMPtee%d",k+1);
    sprintf(histtitleRatio,"random DiEMPt ratio for ee %d",k+1);
    RatioGaus[k]=new TH1F(histnameRatio,histtitleRatio,200,0,200);
    
    for(int j=0; j < 200; ++j){
      float u = r.Gaus(mu1[j], sig1[j]);
      RatioGaus[k]->SetBinContent(j,u);
    }
    
    RatioGaus[k]->Write();

  }



  for(int k=0; k < 1000; ++k){
    if(k%50 == 0)printf("processed rho:%d \n", k);//
    sprintf(histnameY, "DiEMPtweightedeeRho%d",k+1);
    sprintf(histtitleY,"random rho ratio for ee %d",k+1);
    sprintf(histnameX, "ffRho%d",k+1);
    sprintf(histtitleX,"random rho ratio for ff %d",k+1);
    RatioeeDiEMPtRatioWeightedRho[k]=new TH1F(histnameY,histtitleY,40,0.0,40.0);
    RatioffRho[k]=new TH1F(histnameX,histtitleX,40,0.0,40.0);    



    for(int j=0; j < 40; ++j){  
      float v = r1.Gaus(mu2[j], sig2[j]);
      RatioeeDiEMPtRatioWeightedRho[k]->SetBinContent(j,v);
      float w = r2.Gaus(mu3[j], sig3[j]);
      RatioffRho[k]->SetBinContent(j,w);
    } 
    RatioeeDiEMPtRatioWeightedRho[k]->Write();
    RatioffRho[k]->Write();
  }

  
  
}


