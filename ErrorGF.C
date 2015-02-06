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
void ErrorGF()
//---------------------------------------------------------------------------------- 
{
  


  // Set parameters:
  TRandom3 r;
  TRandom3 r1;
  
  
  
  

  // Make histograms:

  TFile *Fnew = new TFile("ElectronRatio.root","READ");
  TH1F *candidategfrho=(TH1F*)Fnew->Get("candidategfrho");
  TH1F *candidatefgrho=(TH1F*)Fnew->Get("candidatefgrho");
  
  TFile* fout = new TFile("RandomRatioGF.root","RECREATE");

  TH1F *RatioGFRho[1000];
  TH1F *RatioFGRho[1000];

  char *histnameGF                   = new char[50];
  char *histnameFG                   = new char[50];
  
  char *histtitleGF                  = new char[50];
  char *histtitleFG                  = new char[50];
  



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------

  double mu3[40]          =  {0};
  double sig3[40]         =  {0};

  double mu4[40]          =  {0};
  double sig4[40]         =  {0};
  

  float rho12 =  0.0; 

  for(int q=0 ; q < 40;++q){
    mu3[q]  = candidategfrho->GetBinContent(q);
    sig3[q] = candidategfrho->GetBinError(q);
    mu4[q]  = candidatefgrho->GetBinContent(q); 
    sig4[q] = candidatefgrho->GetBinError(q);  
  }


  if ((rho12 < -1.0) || (rho12 > 1.0)) {
    printf("  ERROR: Correlation factor not in interval [-1,1], as it is %6.2f \n", rho12);
    return;
  }



//----------------------------------------------------------------------------------
// Loop over process:
//----------------------------------------------------------------------------------



  for(int k=0; k < 1000; ++k){
    if(k%50 == 0)printf("processed rho:%d \n", k);//
    sprintf(histnameGF, "GFRhoRandom%d",k+1);
    sprintf(histtitleGF,"random rho ratio for gf rereco %d",k+1);
    RatioGFRho[k]=new TH1F(histnameGF, histtitleGF,40,0.0,40.0);

    sprintf(histnameFG, "FGRhoRandom%d",k+1);
    sprintf(histtitleFG,"random rho ratio for fg rereco %d",k+1);
    RatioFGRho[k]=new TH1F(histnameFG, histtitleFG,40,0.0,40.0);

    for(int j=0; j < 40; ++j){  
      float w = r.Gaus(mu3[j], sig3[j]);
      RatioGFRho[k]->SetBinContent(j,w);
      float f = r1.Gaus(mu4[j], sig4[j]);
      RatioFGRho[k]->SetBinContent(j,f);
    } 
    RatioGFRho[k]->Write();
    RatioFGRho[k]->Write();
  }

  
  
}


