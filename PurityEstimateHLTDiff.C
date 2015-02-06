#define PurityEstimateHLTDiff_cxx
#include "PurityEstimateHLTDiff.h"



void
photonEffectiveAreas(double _eta, double* _effA)
{
  double& effACH(_effA[0]);
  double& effANH(_effA[1]);
  double& effAPh(_effA[2]);

  // Source: CutBasedPhotonID2012 twiki
  if(_eta < 1.){
    effACH = 0.012;
    effANH = 0.03;
    effAPh = 0.148;
  }
  else if(_eta < 1.479){
    effACH = 0.010;
    effANH = 0.057;
    effAPh = 0.13;
  }
}

float dRCalc(float etaLead, float phiLead, float etaTrail, float phiTrail){
    
  float dphi = fabs(phiLead - phiTrail);
  if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
    /*if (debug){
        cout << "LeadPhi: " << phiLead << endl;
        cout << "TrailPhi: " << phiTrail << endl;
        cout << "dphi: " << dphi << endl;
        cout << "deta: " << fabs(etaLead-etaTrail) << endl;
    }*/
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

void PurityEstimateHLTDiff::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L PurityEstimate.C
//      Root > PurityEstimate t
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

   TFile *f1 = new TFile("BhulBhal.root","RECREATE");
   /* Defining histograms */

   // met histograms
   TH1F* h_met_ffSample_HLT(new TH1F("h_met_ffSample_HLT","ff #slash{E}_{T},Symmetric Pt Cut, different HLT ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eeSample_HLT(new TH1F("h_met_eeSample_HLT","ee #slash{E}_{T},Symmetric Pt Cut, different HLT ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_HLT(new TH1F("h_met_candidate_HLT","candidate #slash{E}_{T},Symmetric Pt Cut, different HLT ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   /*TH1F* h_met_eeSample_RhoReweighted_HLT(new TH1F("h_met_eeSample_RhoReweighted_HLT","ee #slash{E}_{T}, rho reweighted to candidate, different HLT ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_RhoReweighted_HLT(new TH1F("h_met_ffSample_RhoReweighted_HLT","ff #slash{E}_{T}, rho reweighted to candidate, different HLT ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eeSample_diEMPtRhoReweighted_HLT(new TH1F("h_met_eeSample_diEMPtRhoReweighted_HLT","ee #slash{E}_{T}, rho and diEMPt reweighted to candidate, different HLT ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_dR_candidate_HLT(new TH1F("h_dR_candidate_HLT","candidate #Delta R, different HLT ;#Delta R ;Events ", 100, 0., 5.));
   TH1F* h_dPhi_candidate_dr06_HLT(new TH1F("h_dPhi_candidate_dr06_HLT","candidate #Delta#phi,#Delta R > 0.6 different HLT ;#Delta#phi ;Events ", 80, 0., 4.));*/
   // histograms for purity
   // candidate
   /*TH1F* h_et_all_passing_shower_cut_candidate_HLT(new TH1F("h_et_all_passing_shower_cut_candidate_HLT","all passing shower cut,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_HLT(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_HLT","lead passing trail failing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_HLT(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_HLT","lead failing trail passing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_HLT(new TH1F("h_et_all_failing_shower_cut_candidate_HLT","all failing shower cut,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));*/

   TH1F* h_et_all_passing_chargeiso_candidate_HLT(new TH1F("h_et_all_passing_chargeiso_candidate_HLT","all passing charge iso,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_chargeiso_candidate_HLT(new TH1F("h_et_lead_passing_trail_failing_chargeiso_candidate_HLT","lead passing trail failing charge iso, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_chargeiso_candidate_HLT(new TH1F("h_et_lead_failing_trail_passing_chargeiso_candidate_HLT","lead failing trail passing charge iso, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_chargeiso_candidate_HLT(new TH1F("h_et_all_failing_chargeiso_candidate_HLT","all failing charge iso,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));


   // fake
   /*TH1F* h_et_all_passing_shower_cut_fake_HLT(new TH1F("h_et_all_passing_shower_cut_fake_HLT","all passing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_fake_HLT(new TH1F("h_et_lead_passing_trail_failing_shower_cut_fake_HLT","lead passing trail failing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_fake_HLT(new TH1F("h_et_lead_failing_trail_passing_shower_cut_fake_HLT","lead failing trail passing shower cut,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_fake_HLT(new TH1F("h_et_all_failing_shower_cut_fake_HLT","all failing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));*/

   TH1F* h_et_all_passing_chargeiso_fake_HLT(new TH1F("h_et_all_passing_chargeiso_fake_HLT","all passing charge iso,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_chargeiso_fake_HLT(new TH1F("h_et_lead_passing_trail_failing_chargeiso_fake_HLT","lead passing trail failing charge iso, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_chargeiso_fake_HLT(new TH1F("h_et_lead_failing_trail_passing_chargeiso_fake_HLT","lead failing trail passing charge iso, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_chargeiso_fake_HLT(new TH1F("h_et_all_failing_chargeiso_fake_HLT","all failing charge iso,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   // electron
   TH1F* h_et_all_passing_shower_cut_electron_HLT(new TH1F("h_et_all_passing_shower_cut_electron_HLT","all passing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_electron_HLT(new TH1F("h_et_lead_passing_trail_failing_shower_cut_electron_HLT","lead passing trail failing shower cut, Symmetric Pt Cut,different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_electron_HLT(new TH1F("h_et_lead_failing_trail_passing_shower_cut_electron_HLT","lead failing trail passing shower cut,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_electron_HLT(new TH1F("h_et_all_failing_shower_cut_electron_HLT","all failing shower cut,Symmetric Pt Cut, different HLT ;Pt (GeV);Events / GeV", 200, 0., 200.));

   // 2D histograms
  /* TH2F* h_leadeta_trailphi_candidate_HLT(new TH2F("h_leadeta_trailphi_candidate_HLT","lead #eta vs trail #phi, different HLT ;lead #eta ;trail #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_leadeta_leadphi_candidate_HLT(new TH2F("h_leadeta_leadphi_candidate_HLT","lead #eta vs lead #phi, different HLT ;lead #eta ;lead #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_traileta_leadphi_candidate_HLT(new TH2F("h_traileta_leadphi_candidate_HLT","trail #eta vs lead #phi, different HLT ;trail #eta ;lead #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_traileta_trailphi_candidate_HLT(new TH2F("h_traileta_trailphi_candidate_HLT","trail #eta vs trail #phi, different HLT ;trail #eta ;trail #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_leadeta_traileta_candidate_HLT(new TH2F("h_leadeta_traileta_candidate_HLT","lead #eta vs trail #eta, different HLT ;lead #eta ;trail #eta ", 80, -2., 2.,80,-2.,2.));
   TH2F* h_leadphi_trailphi_candidate_HLT(new TH2F("h_leadphi_trailphi_candidate_HLT","lead #phi vs trail #phi, different HLT ;lead #phi ;trail #phi ", 160, -4., 4.,160,-4.,4.));
   TH2F* h_leadpt_trailpt_candidate_HLT(new TH2F("h_leadpt_trailpt_candidate_HLT","lead Pt vs trail Pt ;lead Pt, different HLT ;trail Pt ", 300, 0., 300.,300,0.,300.));*/
   // rho histograms
   TH1F* h_rho_ffSample_HLT(new TH1F("h_rho_ffSample_HLT","ff Rho25, Symmetric Pt Cut, different HLT ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_eeSample_HLT(new TH1F("h_rho_eeSample_HLT","ee Rho25, Symmetric Pt Cut, different HLT ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_candidate_HLT(new TH1F("h_rho_candidate_HLT","candidate Rho25, Symmetric Pt Cut, different HLT ;Rho25;Events ", 50, 0., 50.));
   // nVtx histograms
   //TH1F* h_nVtx_ffSample(new TH1F("h_nVtx_ffSample","ff nVtx ;nVtx;Events ", 50, 0., 50.));
   //TH1F* h_nVtx_eeSample(new TH1F("h_nVtx_eeSample","ee nVtx ;nVtx;Events ", 50, 0., 50.));
   //TH1F* h_nVtx_candidate(new TH1F("h_nVtx_candidate","candidate nVtx ;nVtx;Events ", 50, 0., 50.));

   // diEMPt histograms
   TH1F* h_diEMPt_eeSample_HLT(new TH1F("h_diEMPt_eeSample_HLT","ee diEMPt, Symmetric Pt Cut, different HLT ;diEMPt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_diEMPt_candidate_HLT(new TH1F("h_diEMPt_candidate_HLT","candidate diEMPt, Symmetric Pt Cut, different HLT;diEMPt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_diEMPt_ffSample_HLT(new TH1F("h_diEMPt_ffSample_HLT","ff diEMPt, Symmetric Pt Cut, different HLT ;diEMPt (GeV);Events / GeV", 200, 0., 200.));

   // Saving the error weights //
   /*h_met_eeSample_HLT->Sumw2();
   h_met_ffSample_HLT->Sumw2();
   h_met_candidate_HLT->Sumw2();
   h_rho_eeSample_HLT->Sumw2();
   h_rho_ffSample_HLT->Sumw2();
   h_rho_candidate_HLT->Sumw2();
   h_diEMPt_eeSample_HLT->Sumw2();
   h_diEMPt_candidate_HLT->Sumw2();
   h_diEMPt_ffSample_HLT->Sumw2();*/
   

   // getting the reweighting //

   /*TFile *Fnew = new TFile("SymmetricCheckHLTDiff.root","READ");
   TH1F *h_rho_eeSample_HLT = (TH1F*)Fnew->Get("h_rho_eeSample_HLT");
   TH1F *h_rho_ffSample_HLT = (TH1F*)Fnew->Get("h_rho_ffSample_HLT");
   TH1F *h_rho_candidate_HLT = (TH1F*)Fnew->Get("h_rho_candidate_HLT");
   h_rho_eeSample_HLT->Scale(1./h_rho_eeSample_HLT->Integral(0,50));
   h_rho_ffSample_HLT->Scale(1./h_rho_ffSample_HLT->Integral(0,50));
   h_rho_candidate_HLT->Scale(1./h_rho_candidate_HLT->Integral(0,50));
   TH1F *candidateOvereeRho = new TH1F("candidateOvereeRho", "candidate  to ee Rho ratio", 50, 0., 50.);
   candidateOvereeRho->Sumw2(); 
   candidateOvereeRho->Divide(h_rho_candidate_HLT,h_rho_eeSample_HLT);
   TH1F *candidateOverffRho = new TH1F("candidateOverffRho", "candidate  to ff Rho ratio", 50, 0., 50.);
   candidateOverffRho->Sumw2();
   candidateOverffRho->Divide(h_rho_candidate_HLT,h_rho_ffSample_HLT);

   //TFile *Fnew2 = new TFile("Weighted.root","READ");
   TH1F *h_diEMPt_eeSample_HLT = (TH1F*)Fnew->Get("h_diEMPt_eeSample_HLT");
   TH1F *h_diEMPt_candidate_HLT = (TH1F*)Fnew->Get("h_diEMPt_candidate_HLT");
   h_diEMPt_eeSample_HLT->Scale(1./h_diEMPt_eeSample_HLT->Integral(0,200));
   h_diEMPt_candidate_HLT->Scale(1./h_diEMPt_candidate_HLT->Integral(0,200));
   TH1F *candidateOvereediEMPt = new TH1F("candidateOvereediEMPt", "candidate  to ee diEMPt ratio", 200, 0., 200.);
   candidateOvereediEMPt->Sumw2();
   candidateOvereediEMPt->Divide(h_diEMPt_candidate_HLT,h_diEMPt_eeSample_HLT);*/
   
    
  
 

   if (fChain == 0) return;
   
   

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if(jentry%100000==0)printf("Entry processed: %lld\n", jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      /* /////// Photon, electron and fake Selection /////// */

      int   realphoton = 0;
      int   candidate  = 0;
      int   fakephoton = 0;
      int   fake       = 0;
      int   elephoton  = 0;
      int   electron   = 0;
      for(unsigned int i=0;i < photon_neutraliso->size();++i){
        float absEta = fabs(photon_eta->at(i));
        double effA[3];
        photonEffectiveAreas(absEta, effA);
        float theta = 2*atan(exp(-photon_eta->at(i)));
        float pt = fabs(photon_e->at(i)*sin(theta)); 

        bool neutral   = photon_neutraliso->at(i) - rho * effA[1] - 0.04 * pt < 3.5;
        bool charged   = photon_chargeiso->at(i) - rho * effA[0] < 2.6;
        bool photoniso = photon_photoniso->at(i) - rho * effA[2] - 0.005 * pt < 1.3;
        bool showercut = photon_showershape->at(i) < 0.011 ; 
        bool pixelcut  = photon_pixelseed->at(i) == 0;
        bool SymPtCut  = pt > 40;
  
        if(neutral && charged && photoniso && showercut && pixelcut && SymPtCut  ) realphoton++;
        if(neutral && !charged && photoniso && showercut && pixelcut && SymPtCut) fakephoton++;
        if(neutral && charged && photoniso && showercut && !pixelcut && SymPtCut) elephoton++;

        if(neutral && showercut && photoniso && pixelcut  && SymPtCut) candidate++;
        if(neutral && !showercut && photoniso && pixelcut && SymPtCut)  fake++;
        if(neutral && charged && photoniso && !pixelcut && SymPtCut)  electron++;
       
        
      }
      if(realphoton >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        //h_dR_candidate_HLT->Fill(DR);
        float DPHI=dPhiCalc(photon_phi->at(0),photon_phi->at(1));
        //h_dPhi_candidate_HLT->Fill(DPHI);
        if(DR > 0.6){
          //printf("Photon phi: %f\n",photon_phi->at(0));
          /*h_leadeta_trailphi_candidate_HLT->Fill(photon_eta->at(0),photon_phi->at(1));
          h_leadeta_leadphi_candidate_HLT->Fill(photon_eta->at(0),photon_phi->at(0));
          h_traileta_leadphi_candidate_HLT->Fill(photon_eta->at(1),photon_phi->at(0));
          h_traileta_trailphi_candidate_HLT->Fill(photon_eta->at(1),photon_phi->at(1));
          h_leadeta_traileta_candidate_HLT->Fill(photon_eta->at(0),photon_eta->at(1));
          h_leadphi_trailphi_candidate_HLT->Fill(photon_phi->at(0),photon_phi->at(1));
          float leadPt  = findPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0));
          float trailPt = findPt(photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
          h_leadpt_trailpt_candidate_HLT->Fill(leadPt,trailPt);
          

          h_dPhi_candidate_dr06_HLT->Fill(DPHI);
          //printf("realphoton Entry: %lld\n", jentry);*/
          h_met_candidate_HLT->Fill(met_et);
          float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
          h_diEMPt_candidate_HLT->Fill(DIEMPT);
          h_rho_candidate_HLT->Fill(rho);
        }
      }
      if(fakephoton >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          //printf("fakephoton Entry: %lld\n", jentry);
          h_met_ffSample_HLT->Fill(met_et);
          float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
          h_diEMPt_ffSample_HLT->Fill(DIEMPT);
          h_rho_ffSample_HLT->Fill(rho);
          //int Bin=candidateOverffRho->FindBin(rho);
          //float Weight=candidateOverffRho->GetBinContent(Bin);
          //h_met_ffSample_RhoReweighted_HLT->Fill(met_et, Weight);
          
        }
      }
      if(elephoton >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          //printf("elephoton Entry: %lld\n", jentry);
          h_met_eeSample_HLT->Fill(met_et);
          float DIEMPT=findDiEMPt(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
          //int Bin=candidateOvereeRho->FindBin(rho);
          //float Weight=candidateOvereeRho->GetBinContent(Bin);
          //h_met_eeSample_RhoReweighted_HLT->Fill(met_et, Weight);
          //int Bin2=candidateOvereediEMPt->FindBin(DIEMPT);
          //float Weight2=candidateOvereediEMPt->GetBinContent(Bin2);
          //h_met_eeSample_diEMPtRhoReweighted_HLT->Fill(met_et, Weight*Weight2);
          h_diEMPt_eeSample_HLT->Fill(DIEMPT);
          h_rho_eeSample_HLT->Fill(rho);
        }
      }
      // calculating purity from matrix method
      if(candidate >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        float absEta = fabs(photon_eta->at(0));
        double effA[3];
        photonEffectiveAreas(absEta, effA);
        float absEta2 = fabs(photon_eta->at(1));
        double effA2[3];
        photonEffectiveAreas(absEta2, effA2);
        
        if(DR>0.6){
          
          /*if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)<0.011)h_et_all_passing_shower_cut_candidate_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)>0.011)h_et_lead_passing_trail_failing_shower_cut_candidate_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)<0.011)h_et_lead_failing_trail_passing_shower_cut_candidate_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)>0.011)h_et_all_failing_shower_cut_candidate_HLT->Fill(photon_e->at(0));*/
          bool charged1   = photon_chargeiso->at(0) - rho * effA[0] < 2.6;
          bool charged2   = photon_chargeiso->at(1) - rho * effA2[0] < 2.6;
          if(charged1 && charged2)h_et_all_passing_chargeiso_candidate_HLT->Fill(photon_e->at(0));
          if(charged1 && !charged2)h_et_lead_passing_trail_failing_chargeiso_candidate_HLT->Fill(photon_e->at(0));
          if(!charged1 && charged2)h_et_lead_failing_trail_passing_chargeiso_candidate_HLT->Fill(photon_e->at(0));
          if(!charged1 && !charged2)h_et_all_failing_chargeiso_candidate_HLT->Fill(photon_e->at(0));
        }
      }

      if(fake >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        float absEta = fabs(photon_eta->at(0));
        double effA[3];
        photonEffectiveAreas(absEta, effA);
        float absEta2 = fabs(photon_eta->at(1));
        double effA2[3];
        photonEffectiveAreas(absEta2, effA2);
        if(DR>0.6){
          bool charged1   = photon_chargeiso->at(0) - rho * effA[0] < 2.6;
          bool charged2   = photon_chargeiso->at(1) - rho * effA2[0] < 2.6;
          if(charged1 && charged2)h_et_all_passing_chargeiso_fake_HLT->Fill(photon_e->at(0));
          if(charged1 && !charged2)h_et_lead_passing_trail_failing_chargeiso_fake_HLT->Fill(photon_e->at(0));
          if(!charged1 && charged2)h_et_lead_failing_trail_passing_chargeiso_fake_HLT->Fill(photon_e->at(0));
          if(!charged1 && !charged2)h_et_all_failing_chargeiso_fake_HLT->Fill(photon_e->at(0));
          /*if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)<0.011)h_et_all_passing_shower_cut_fake_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)>0.011)h_et_lead_passing_trail_failing_shower_cut_fake_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)<0.011)h_et_lead_failing_trail_passing_shower_cut_fake_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)>0.011)h_et_all_failing_shower_cut_fake_HLT->Fill(photon_e->at(0));*/
        }
      }

      if(electron >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        float MASS=findMass(photon_e->at(0),photon_eta->at(0),photon_phi->at(0),photon_e->at(1),photon_eta->at(1),photon_phi->at(1));
        //printf("MASS:%f for entry:%lld\n",MASS,jentry);
        if(DR>0.6 && MASS>=75 && MASS<=105){
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)<0.011)h_et_all_passing_shower_cut_electron_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)<0.011 && photon_showershape->at(1)>0.011)h_et_lead_passing_trail_failing_shower_cut_electron_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)<0.011)h_et_lead_failing_trail_passing_shower_cut_electron_HLT->Fill(photon_e->at(0));
          if(photon_showershape->at(0)>0.011 && photon_showershape->at(1)>0.011)h_et_all_failing_shower_cut_electron_HLT->Fill(photon_e->at(0));
        }
      }
      // if (Cut(ientry) < 0) continue;*/
      
   }
   f1->cd();
   f1->Write();
   
}
