#define SusyLimit_cxx
#include "SusyLimit.h"
#include <TGraphErrors.h> 



void
photonEffectiveAreas(double _eta, double* _effA)
{
  double& effACH(_effA[0]);
  double& effANH(_effA[1]);
  double& effAPh(_effA[2]);
  double& effAW(_effA[3]);

  // Source: CutBasedPhotonID2012 twiki
  if(_eta < 1.){
    effACH = 0.012;
    effANH = 0.03;
    effAPh = 0.148;
    effAW  = 0.075;
  }
  else if(_eta < 1.479){
    effACH = 0.010;
    effANH = 0.057;
    effAPh = 0.13;
    effAW  = 0.0617;
  }
}

float dRCalc(float etaLead, float phiLead, float etaTrail, float phiTrail){
    
  float dphi = fabs(phiLead - phiTrail);
  if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  float deta = fabs(etaLead - etaTrail);
  float dR = sqrt(deta*deta + dphi*dphi);
  return dR;
    
}

float dPhiCalc(float phiLead, float phiTrail){
  float dphi = fabs(phiLead - phiTrail);
  if(dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  return dphi;
}

float findDiEMPt(float ELead, float EtaLead, float PhiLead, float ETrail, float EtaTrail, float PhiTrail){
  float theta1 = 2*atan(exp(-EtaLead));
  float theta2 = 2*atan(exp(-EtaTrail));
  float PX1 = ELead*sin(theta1)*cos(PhiLead);
  float PY1 = ELead*sin(theta1)*sin(PhiLead);
  float PX2 = ETrail*sin(theta2)*cos(PhiTrail);
  float PY2 = ETrail*sin(theta2)*sin(PhiTrail);
  float DiEMPt = sqrt((PX1+PX2)*(PX1+PX2)+(PY1+PY2)*(PY1+PY2));
  return DiEMPt;

}

float findMass(float ELead, float EtaLead, float PhiLead, float ETrail, float EtaTrail, float PhiTrail){
  float theta1 = 2*atan(exp(-EtaLead));
  float theta2 = 2*atan(exp(-EtaTrail));
  float PX1 = ELead*sin(theta1)*cos(PhiLead);
  float PY1 = ELead*sin(theta1)*sin(PhiLead);
  float PX2 = ETrail*sin(theta2)*cos(PhiTrail);
  float PY2 = ETrail*sin(theta2)*sin(PhiTrail);
  float PZ1 = ELead*cos(theta1);
  float PZ2 = ETrail*cos(theta2);
  float Mass = sqrt((ELead+ETrail)*(ELead+ETrail)-((PX1+PX2)*(PX1+PX2)+(PY1+PY2)*(PY1+PY2)+(PZ1+PZ2)*(PZ1+PZ2)));
  return Mass;

}

float findPt(float Energy, float Eta, float Phi){
  float theta = 2*atan(exp(-Eta));
  float PX1 = Energy*sin(theta)*cos(Phi);
  float PY1 = Energy*sin(theta)*sin(Phi);
  float PT  = sqrt(PX1*PX1+PY1*PY1);
  return PT;
}

void SusyLimit::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SusyLimit.C
//      Root > SusyLimit t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   TFile *f1 = new TFile("Signal_GluinoXXX_LSPYYY_AsymmetricPt.root","RECREATE");
   /* Defining histograms */
  
   
   //TH1F *h2 = new TH1F("h2","photon over electron",binnum,bins);
   Float_t bins[]={0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,150,300};
   Int_t binnum = sizeof(bins)/sizeof(Float_t) - 1;
   
   TH1F* h_rho_Signal(new TH1F("h_rho_Signal","candidate rho, Asymmetric Pt ;Rho;Events", 50,0,50));
   TH1F* h_met_candidate_rhoreweighted_MET60_100(new TH1F("h_met_candidate_rhoreweighted_MET60_100","candidate rho reweighted #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));
   TH1F* h_met_candidate_rhoreweighted_MET100_Up(new TH1F("h_met_candidate_rhoreweighted_MET100_Up","candidate rho reweighted #slash{E}_{T}, Asymmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));


   h_rho_Signal->Sumw2();
   h_met_candidate_rhoreweighted_MET60_100->Sumw2();
   h_met_candidate_rhoreweighted_MET100_Up->Sumw2();
 
   if (fChain == 0) return;
   FILE *SignalLimit = fopen("file_SusyLimit_Asymmetric.txt","a");

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   printf("Total no of entry: %lld\n",nentries);
   for (Long64_t jentry=0; jentry<nentries;jentry++){
      Long64_t ientry = LoadTree(jentry);
      if(jentry%10000==0)printf("Entry processed: %lld\n", jentry);//
      if (ientry < 0) break;
      //if(jentry>600000)break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int   realphoton = 0;
      vector<int> realphoton_tracker;
      for(size_t i=0;i < photon_neutraliso->size();++i){
        float absEta = fabs(photon_eta->at(i));
        double effA[4];
        photonEffectiveAreas(absEta, effA);
        float theta = 2*atan(exp(-photon_eta->at(i)));
        float pt = fabs(photon_e->at(i)*sin(theta)); 

        bool neutral     = (photon_neutraliso->at(i) - rho * effA[1] - 0.04 * pt) < 3.5;
        bool charged     = (photon_chargeiso->at(i) - rho * effA[0]) < 2.6;
        bool photoniso   = (photon_photoniso->at(i) - rho * effA[2] - 0.005 * pt) < 1.3;
        bool showercut   = photon_showershape->at(i) < 0.011 ;
        bool pixelcut    = photon_pixelseed->at(i) == 0;
	//bool SymPtCut    = pt > 40; //for symmetric pt cut
	bool SymPtCut    = pt > 0; // for asymmetric p_t cut
  
        if(neutral && charged && photoniso && showercut && pixelcut && SymPtCut ){
          realphoton++;
          realphoton_tracker.push_back(i);
        } //worstiso not used
      }
      if(realphoton >= 2){
        float DR=dRCalc(photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
        float MASS=findMass(photon_e->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_e->at(realphoton_tracker.at(1)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
        if(DR > 0.6 && (MASS<=80 || MASS>=100)){
          h_rho_Signal->Fill(rho); 
        }  
      }   
   }
   
   Float_t x1 = h_rho_Signal->Integral();
   h_rho_Signal->Scale(1./x1);
   //TFile *Fnew = new TFile("DiffFake_eeRemoved_SymmetricPt.root","READ"); // for Symmetric P_T CutBasedPhotonID2012
   TFile *Fnew = new TFile("Notimp_eeRemoved.root", "READ"); // for Asymmetric P_T 
   TH1F *h_rho_candidate = (TH1F*)Fnew->Get("h_rho_candidate");
   Float_t x2 = h_rho_candidate->Integral();
   h_rho_candidate->Scale(1./x2);
   TH1F* candidatesignalrho = new TH1F("candidatesignalrho","candidate over signal rho",50,0,50);
   candidatesignalrho->Sumw2();
   candidatesignalrho->Divide(h_rho_Signal, h_rho_candidate); 
   
   for (Long64_t jentry=0; jentry<nentries;jentry++){
      Long64_t ientry = LoadTree(jentry);
      if(jentry%10000==0)printf("Again Entry processed: %lld\n", jentry);//
      if (ientry < 0) break;
      //if(jentry>600000)break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      int   realphoton = 0;
      vector<int> realphoton_tracker;
      for(size_t i=0;i < photon_neutraliso->size();++i){
        float absEta = fabs(photon_eta->at(i));
        double effA[4];
        photonEffectiveAreas(absEta, effA);
        float theta = 2*atan(exp(-photon_eta->at(i)));
        float pt = fabs(photon_e->at(i)*sin(theta)); 

        bool neutral     = (photon_neutraliso->at(i) - rho * effA[1] - 0.04 * pt) < 3.5;
        bool charged     = (photon_chargeiso->at(i) - rho * effA[0]) < 2.6;
        bool photoniso   = (photon_photoniso->at(i) - rho * effA[2] - 0.005 * pt) < 1.3;
        bool showercut   = photon_showershape->at(i) < 0.011 ;
        bool pixelcut    = photon_pixelseed->at(i) == 0;
	//bool SymPtCut    = pt > 40; // for symmetric P_T cut
        bool SymPtCut    = pt > 0; // for asymmetric P_T cut
  
        if(neutral && charged && photoniso && showercut && pixelcut && SymPtCut ){
          realphoton++;
          realphoton_tracker.push_back(i);
        } //worstiso not used
      }
      if(realphoton >= 2){
        float DR=dRCalc(photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
        float MASS=findMass(photon_e->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_e->at(realphoton_tracker.at(1)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
        if(DR > 0.6 && (MASS<=80 || MASS>=100)){
          int BinC=candidatesignalrho->FindBin(rho);
          float WeightC=candidatesignalrho->GetBinContent(BinC);
          if(met_et > 60 && met_et <= 100)h_met_candidate_rhoreweighted_MET60_100->Fill(met_et, WeightC);
          if(met_et > 100)h_met_candidate_rhoreweighted_MET100_Up->Fill(met_et, WeightC); 
        }  
      }   
   }
   float NumberEvent60  = h_met_candidate_rhoreweighted_MET60_100->Integral()/nentries;
   float NumberEvent100 = h_met_candidate_rhoreweighted_MET100_Up->Integral()/nentries;
   float Lumi(19.6), Xsec(0), XsecError(0);
   int MGlu(GGG), GluMass(0);
   int MLSP(LLL);
   FILE *Cross=fopen("CrossSec.txt","r");
   char oneline[100];
   while(fgets(oneline, sizeof oneline, Cross) != NULL){
     fscanf(Cross, "%d\t%f\t\t%f", &GluMass, &Xsec, &XsecError);
     //printf("%d\t%f\t\t%f\n",GluMass,Xsec,XsecError);
     if(GluMass == MGlu){
       //printf("%d\t%f\t%f\n",GluMass,Xsec,XsecError);
       break;
     }
   }
   fclose(Cross);
   double Expect60  = NumberEvent60*Lumi*Xsec;
   double Expect100 = NumberEvent100*Lumi*Xsec;
   //printf("NumberEvent60:%f Lumi:%f, Xsec:%f\n",NumberEvent60, Lumi, Xsec);
   //printf("MGlu:%d MLSP:%d Expect60:%f expect100:%f\n",MGlu, MLSP, Expect60, Expect100);
   
   fprintf(SignalLimit,"%d\t%d\t%f\t%f\t%f\n",MGlu,MLSP,Expect60,Expect100,XsecError);
    
   fclose(SignalLimit);
   f1->cd();
   f1->Write();
   f1->Close();
   
}
