#include <TH1.h>
#include <TFile.h>




void RatioProducer(){  

  //////// electron reweighting  ////////////
   TFile f8("DiffFake_eeRemoved_SymmetricPt.root");
   /*h_diEMPt_candidate->Sumw2();
   Float_t y1=h_diEMPt_candidate->Integral(0,200);
   h_diEMPt_candidate->Scale(1.0/y1);*/


   h_rho_candidate->Sumw2();
   Float_t d2=h_rho_candidate->Integral(0,50);
   h_rho_candidate->Scale(1./d2);
   TH1F *h_rho_source = (TH1F*)h_rho_candidate->Clone("h_rho_source");


   TFile f9("Signal_Gluino1350_LSP375_SymmetricPt.root");
   h_rho_candidate->Sumw2();
   Float_t d3=h_rho_candidate->Integral(0,50);
   h_rho_candidate->Scale(1./d3);
   

   /*h_rho_eeSample->Sumw2();
   Float_t d3=h_rho_eeSample->Integral(0,50);
   h_rho_eeSample->Scale(1./d3);

   h_rho_ffSample->Sumw2();
   Float_t d4=h_rho_ffSample->Integral(0,50);
   h_rho_ffSample->Scale(1./d4);

   h_rho_gfSample->Sumw2();
   Float_t d5=h_rho_gfSample->Integral(0,50);
   h_rho_gfSample->Scale(1./d5);

   h_rho_fgSample->Sumw2();
   Float_t d6=h_rho_fgSample->Integral(0,50);
   h_rho_fgSample->Scale(1./d6);*/

   



   /*h_rho_eeSample_diEMPtReweighted->Sumw2();
   Float_t d1=h_rho_eeSample_diEMPtReweighted->Integral(0,50);
   h_rho_eeSample_diEMPtReweighted->Scale(1.0/d1);
   

   h_diEMPt_ffSample->Sumw2();
   Float_t y12=h_diEMPt_ffSample->Integral(0,200);
   h_diEMPt_ffSample->Scale(1.0/y12);

   h_diEMPt_eeSample->Sumw2();
   Float_t y13=h_diEMPt_eeSample->Integral(0,200);
   h_diEMPt_eeSample->Scale(1./y13);

   TH1F *candidateeediEMPt = new TH1F("candidateeediEMPt","candidate over diEMPt",200,0,200);
   TH1F *candidateeediEMPtweightedrho = new TH1F("candidateeediEMPtweightedrho","candidate over diEMPtweighted ee rho",50,0,50);
   TH1F *candidateffrho = new TH1F("candidateffrho","candidate over ff rho",50,0,50);
   TH1F *candidategfrho = new TH1F("candidategfrho","candidate over gf rho",50,0,50);*/
   TH1F *candidatesignalrho = new TH1F("candidatesignalrho","candidate over signal rho",50,0,50);



   /*candidateeediEMPt->Sumw2();
   candidateeediEMPt->Divide(h_diEMPt_candidate, h_diEMPt_eeSample);
   

   candidateeediEMPtweightedrho->Sumw2();
   candidateeediEMPtweightedrho->Divide(h_rho_candidate, h_rho_eeSample_diEMPtReweighted);


   candidateffrho->Sumw2();
   candidateffrho->Divide(h_rho_candidate, h_rho_ffSample);

   candidategfrho->Sumw2();
   candidategfrho->Divide(h_rho_candidate, h_rho_gfSample);*/
  
   candidatesignalrho->Sumw2();
   candidatesignalrho->Divide(h_rho_source, h_rho_candidate);

   TFile* fout = new TFile("SignalRhoRatioSymmetric.root","RECREATE");
   fout->cd();
   /*candidateeediEMPt->Write();
   candidateeediEMPtweightedrho->Write();
   candidateffrho->Write();
   candidategfrho->Write();*/
   candidatesignalrho->Write();
   
   

   
   
   
   

}
