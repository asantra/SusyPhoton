#define PurityEstimate_cxx
#include "PurityEstimate.h"
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

void PurityEstimate::Loop()
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

   TFile *f1 = new TFile("DiffFake_eeRemoved_SymmetricPt.root","RECREATE");
   /* Defining histograms */
   ////// event no hisotgram ////////
   TH1F* h_met_rereco_all(new TH1F("h_met_rereco_all","rereco #slash{E}_{T}, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_rereco(new TH1F("h_met_rereco","rereco #slash{E}_{T}, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_PhoFake_rereco(new TH1F("h_met_PhoFake_rereco","rereco #slash{E}_{T} gamma-fake, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_FakePho_rereco(new TH1F("h_met_FakePho_rereco","rereco #slash{E}_{T} fake-gamma, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));

   TH1F* h_met_PhoFake_rereco_RhoReweighted(new TH1F("h_met_PhoFake_rereco_RhoReweighted","rereco #slash{E}_{T} gamma-fake RhoReweighted, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_FakePho_rereco_RhoReweighted(new TH1F("h_met_FakePho_rereco_RhoReweighted","rereco #slash{E}_{T} fake-gamma RhoReweighted, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));

   TH1F* h_metX_rereco(new TH1F("h_metX_rereco","rereco #slash{E}_{T} X component, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{TX} (GeV);Events / GeV", 400, -200., 200.));
   TH1F* h_metY_rereco(new TH1F("h_metY_rereco","rereco #slash{E}_{T} Y component, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{TY} (GeV);Events / GeV", 400, -200., 200.));
   TH1F* h_lead_photon_pt_rereco(new TH1F("h_lead_photon_pt_rereco","rereco lead photon pt, Mass 80-100 GeV excluded, Symmetric Pt ;#P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_trail_photon_pt_rereco(new TH1F("h_trail_photon_pt_rereco","rereco trail photon pt, Mass 80-100 GeV excluded, Symmetric Pt ;#P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_lead_photon_phi_rereco(new TH1F("h_lead_photon_phi_rereco","rereco lead photon #Phi, Mass 80-100 GeV excluded, Symmetric Pt ;#Phi ;Events", 80, -4., 4.));
   TH1F* h_trail_photon_phi_rereco(new TH1F("h_trail_photon_phi_rereco","rereco trail photon #Phi, Mass 80-100 GeV excluded, Symmetric Pt ;#Phi ;Events ", 80, -4., 4.));
   TH1F* h_lead_photon_eta_rereco(new TH1F("h_lead_photon_eta_rereco","rereco lead photon #eta, Mass 80-100 GeV excluded, Symmetric Pt ;#eta ;Events", 40, -2., 2.));
   TH1F* h_trail_photon_eta_rereco(new TH1F("h_trail_photon_eta_rereco","rereco trail photon #eta, Symmetric Pt ;#eta ;Events ", 40, -2., 2.));
   TH1F* h_diEMPt_rereco(new TH1F("h_diEMPt_rereco","rereco diEMPt, Mass 80-100 GeV excluded, Symmetric Pt ;DiEM #P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_reduced_met_rereco(new TH1F("h_reduced_met_rereco","rereco MET, Mass 80-100 GeV excluded, Symmetric Pt ;reduced MET (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_jetpt_rereco(new TH1F("h_jetpt_rereco","rereco jet #P_T, Symmetric Pt ; #P_T (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_jetphi_rereco(new TH1F("h_jetphi_rereco","rereco jet #Phi, Symmetric Pt ; #Phi;Events ", 80, -4., 4.));
   TH1F* h_jeteta_rereco(new TH1F("h_jeteta_rereco","rereco jet #eta, Symmetric Pt ; #eta;Events ", 40, -2., 2.));



   TH1F *h_met_eeSample_DiEMPtRhoReweighted_ErrorPropagatedFromDiEMPt[1000];
   TH1F *h_met_eeSample_DiEMPtRhoReweighted_ErrorPropagatedFromRho[1000];

   TH1F *h_met_ffSample_RhoReweighted_ErrorPropagated[1000];
   TH1F *h_met_ffSample_RhoReweightedGF_ErrorPropagated[1000];

   Float_t bins[] = {0,5,10,15,20,25,30,35,40,45,50,60,70,80,100,150,300};
   Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1; 
  
   ///// met histograms  //////
   TH1F* h_met_ffSample(new TH1F("h_met_ffSample","ff #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));
   TH1F* h_met_eeSample(new TH1F("h_met_eeSample","ee #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));
   
   //////  checking between ee and ey Sample ///////
   TH1F* h_eySample_OpeAng(new TH1F("h_eySample_OpeAng","ey opening angle, Mass 80-100 GeV excluded, Symmetric Pt ;#Delta#phi;Events", 60,0,4));
   TH1F* h_eeSample_OpeAng(new TH1F("h_eeSample_OpeAng","ee opening angle, Symmetric Pt ;#Delta#phi;Events", 60,0,4));
   TH1F* h_eySample_OpeAng_Eta(new TH1F("h_eySample_OpeAng_Eta","ey opening angle Eta, Mass 80-100 GeV excluded, Symmetric Pt ;#Delta#eta;Events", 30,0,2));
   TH1F* h_eeSample_OpeAng_Eta(new TH1F("h_eeSample_OpeAng_Eta","ee opening angle Eta, Symmetric Pt ;#Delta#eta;Events", 30,0,2));
   TH1F* h_eySample_DeltaR(new TH1F("h_eySample_DeltaR","ey DeltaR, Mass 80-100 GeV excluded, Symmetric Pt ;#Delta R;Events", 60,0,4));
   TH1F* h_eeSample_DeltaR(new TH1F("h_eeSample_DeltaR","ee DeltaR, Symmetric Pt ;#Delta R;Events", 60,0,4));


   TH1F* h_eySample_LeadPhi(new TH1F("h_eySample_LeadPhi","ey lead #phi distribution, Mass 80-100 GeV excluded, Symmetric Pt ;#phi;Events", 90,-3.2,3.2));
   TH1F* h_eySample_TrailPhi(new TH1F("h_eySample_TrailPhi","ey trail #phi distribution, Mass 80-100 GeV excluded, Symmetric Pt ;#phi;Events", 90,-3.2,3.2));
   TH1F* h_eySample_LeadEta(new TH1F("h_eySample_LeadEta","ey lead #eta distribution, Mass 80-100 GeV excluded, Symmetric Pt ;#eta;Events", 45,-1.5,1.5));
   TH1F* h_eySample_TrailEta(new TH1F("h_eySample_TrailEta","ey trail #eta distribution, Mass 80-100 GeV excluded, Symmetric Pt ;#eta;Events", 45,-1.5,1.5));
   TH1F* h_eeSample_LeadPhi(new TH1F("h_eeSample_LeadPhi","ee lead #phi distribution, Symmetric Pt ;#phi;Events", 90,-3.2,3.2));
   TH1F* h_eeSample_TrailPhi(new TH1F("h_eeSample_TrailPhi","ee trail #phi distribution, Symmetric Pt ;#phi;Events", 90,-3.2,3.2));
   TH1F* h_eeSample_LeadEta(new TH1F("h_eeSample_LeadEta","ee lead #eta distribution, Symmetric Pt ;#eta;Events", 45,-1.5,1.5));
   TH1F* h_eeSample_TrailEta(new TH1F("h_eeSample_TrailEta","ee trail #eta distribution, Symmetric Pt ;#eta;Events", 45,-1.5,1.5));
   

   //// electroweak background ////

   
   //TH1F *h2 = new TH1F("h2","photon over electron",binnum,bins);
   TH1F* h_met_candidate(new TH1F("h_met_candidate","candidate #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));
   TH1F* h_met_candidate_MET60_100(new TH1F("h_met_candidate_MET60_100","candidate #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));
   TH1F* h_met_candidate_MET100_Up(new TH1F("h_met_candidate_MET100_Up","candidate #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum,bins));

   TH1F* h_met_eySample(new TH1F("h_met_eySample","ey #slash{E}_{T}, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));
   TH1F* h_met_eySample_MET60_100(new TH1F("h_met_eySample_MET60_100","ey #slash{E}_{T}, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));
   TH1F* h_met_eySample_MET100_Up(new TH1F("h_met_eySample_MET100_Up","ey #slash{E}_{T}, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));
   TH1F* h_mass_eySample(new TH1F("h_mass_eySample","Invariant Mass for e+#gamma, Mass 80-100 GeV excluded, Symmetric Pt ;Invariant Mass (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_mass_eeSample(new TH1F("h_mass_eeSample","Invariant Mass for e+e, Symmetric Pt ;Invariant Mass (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_mass_candidate(new TH1F("h_mass_candidate","Invariant Mass for #gamma+#gamma, Symmetric Pt ;Invariant Mass (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eySample_Scaled(new TH1F("h_met_eySample_Scaled","ey #slash{E}_{T} scaled with fake rate, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));
   TH1F* h_met_eySample_Scaled_highf(new TH1F("h_met_eySample_Scaled_highf","ey #slash{E}_{T} scaled with fake rate, high f, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));
   TH1F* h_met_eySample_Scaled_lowf(new TH1F("h_met_eySample_Scaled_lowf","ey #slash{E}_{T} scaled with fake rate, low f, Mass 80-100 GeV excluded, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));
   TH1F* h_mass_eeSample_Insideey(new TH1F("h_mass_eeSample_Insideey","Invariant Mass for e+e in e-y, Symmetric Pt ;Invariant Mass (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_mass_yySample_Insideey(new TH1F("h_mass_yySample_Insideey","Invariant Mass for y+y in e-y, Symmetric Pt ;Invariant Mass (GeV);Events / GeV", 200, 0., 200.));




   TH1F* h_met_candidate_NoCut(new TH1F("h_met_candidate_NoCut","candidate No Cut #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoIsoCut(new TH1F("h_met_candidate_PhoIsoCut","candidate with PhoIsoCut #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoNeutralIso(new TH1F("h_met_candidate_PhoNeutralIso","candidate with Photon and Neutral iso #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoNeutralChargedIso(new TH1F("h_met_candidate_PhoNeutralChargedIso","candidate with Photon , Charged and Neutral iso #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_candidate_PhoNeutralChargedWorstIso(new TH1F("h_met_candidate_PhoNeutralChargedWorstIso","candidate with Photon, Charged , Worst and Neutral iso #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));

   /////// mixing of two mets  ///////
   TH1F* h_met_eeSample_PurityReweighted(new TH1F("h_met_eeSample_PurityReweighted","ee Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eeSample_MET60_100(new TH1F("h_met_eeSample_MET60_100","ee Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_eeSample_MET100_Up(new TH1F("h_met_eeSample_MET100_Up","ee Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_PurityReweighted(new TH1F("h_met_ffSample_PurityReweighted","ff Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_MET60_100(new TH1F("h_met_ffSample_MET60_100","ff Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_MET100_Up(new TH1F("h_met_ffSample_MET100_Up","ff Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_PhoFake_PurityReweighted(new TH1F("h_met_PhoFake_PurityReweighted","yf Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_FakePho_PurityReweighted(new TH1F("h_met_FakePho_PurityReweighted","fy Sample Purity Reweighted #slash{E}_{T}, Symmetric Pt ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));


   TH1F *h_met_PhoFake_rereco_ErrorFromRho[1000];
   TH1F *h_met_FakePho_rereco_ErrorFromRho[1000];
   

   ///// reduced met histograms //////
   
   TH1F* h_reduced_met_ffSample(new TH1F("h_reduced_met_ffSample","ff reduced #slash{E}_{T}, Symmetric Pt ;reduced #slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_reduced_met_eeSample(new TH1F("h_reduced_met_eeSample","ee reduced #slash{E}_{T}, Symmetric Pt ;reduced #slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_reduced_met_eySample(new TH1F("h_reduced_met_eySample","ey reduced #slash{E}_{T}, Symmetric Pt ;reduced #slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));


   TH1F* h_met_eeSample_RhoReweighted(new TH1F("h_met_eeSample_RhoReweighted","ee #slash{E}_{T}, rho reweighted to candidate ;#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_met_ffSample_RhoReweighted(new TH1F("h_met_ffSample_RhoReweighted","ff #slash{E}_{T}, rho reweighted to candidate;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));// changing to variable binning
   TH1F* h_met_eeSample_diEMPtRhoReweighted(new TH1F("h_met_eeSample_diEMPtRhoReweighted","ee #slash{E}_{T}, rho and diEMPt reweighted to candidate;#slash{E}_{T} (GeV);Events / GeV", binnum, bins));// changing to variable binning
   /*TH1F* h_dR_candidate(new TH1F("h_dR_candidate","candidate #Delta R, Symmetric Pt;#Delta R ;Events ", 100, 0., 5.));
   TH1F* h_dPhi_candidate(new TH1F("h_dPhi_candidate","candidate #Delta#phi, Symmetric Pt ;#Delta#phi ;Events ", 80, 0., 4.));
   TH1F* h_dPhi_candidate_dr06(new TH1F("h_dPhi_candidate_dr06","candidate #Delta#phi,#Delta R > 0.6 , Symmetric Pt;#Delta#phi ;Events ", 80, 0., 4.));*/
   // histograms for purity
   ///// candidate //////
   TH1F* h_et_all_passing_shower_cut_candidate(new TH1F("h_et_all_passing_shower_cut_candidate","all passing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate","lead passing trail failing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate","lead failing trail passing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate(new TH1F("h_et_all_failing_shower_cut_candidate","all failing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// varying purity with MET bin MET <= 20  //////
   TH1F* h_et_all_passing_shower_cut_candidate_20(new TH1F("h_et_all_passing_shower_cut_candidate_20","all passing shower cut,MET <= 20, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_20(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_20","lead passing trail failing shower cut,MET <= 20, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_20(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_20","lead failing trail passing shower cut,MET <= 20, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_20(new TH1F("h_et_all_failing_shower_cut_candidate_20","all failing shower cut,MET <= 20, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// varying purity with MET bin 20 < MET <= 40  //////
   TH1F* h_et_all_passing_shower_cut_candidate_40(new TH1F("h_et_all_passing_shower_cut_candidate_40","all passing shower cut,20 < MET <= 40, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_40(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_40","lead passing trail failing shower cut,20 < MET <= 40, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_40(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_40","lead failing trail passing shower cut,20 < MET <= 40, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_40(new TH1F("h_et_all_failing_shower_cut_candidate_40","all failing shower cut,20 < MET <= 40, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// varying purity with MET bin 40 < MET <= 60  //////
   TH1F* h_et_all_passing_shower_cut_candidate_60(new TH1F("h_et_all_passing_shower_cut_candidate_60","all passing shower cut,40 < MET <= 60, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_60(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_60","lead passing trail failing shower cut,40 < MET <= 60, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_60(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_60","lead failing trail passing shower cut,40 < MET <= 60, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_60(new TH1F("h_et_all_failing_shower_cut_candidate_60","all failing shower cut,40 < MET <= 60, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// varying purity with MET bin 60 < MET <= 80  //////
   TH1F* h_et_all_passing_shower_cut_candidate_80(new TH1F("h_et_all_passing_shower_cut_candidate_80","all passing shower cut,60 < MET <= 80, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_80(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_80","lead passing trail failing shower cut,60 < MET <= 80, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_80(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_80","lead failing trail passing shower cut,60 < MET <= 80, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_80(new TH1F("h_et_all_failing_shower_cut_candidate_80","all failing shower cut,60 < MET <= 80, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// varying purity with MET bin 80 < MET <= 100  //////
   TH1F* h_et_all_passing_shower_cut_candidate_100(new TH1F("h_et_all_passing_shower_cut_candidate_100","all passing shower cut,80 < MET <= 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_100(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_100","lead passing trail failing shower cut,80 < MET <= 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_100(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_100","lead failing trail passing shower cut,80 < MET <= 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_100(new TH1F("h_et_all_failing_shower_cut_candidate_100","all failing shower cut,80 < MET <= 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// varying purity with MET bin  MET > 100  //////
   TH1F* h_et_all_passing_shower_cut_candidate_100Up(new TH1F("h_et_all_passing_shower_cut_candidate_100Up","all passing shower cut,MET > 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_candidate_100Up(new TH1F("h_et_lead_passing_trail_failing_shower_cut_candidate_100Up","lead passing trail failing shower cut,MET > 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_candidate_100Up(new TH1F("h_et_lead_failing_trail_passing_shower_cut_candidate_100Up","lead failing trail passing shower cut,MET > 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_candidate_100Up(new TH1F("h_et_all_failing_shower_cut_candidate_100Up","all failing shower cut,MET > 100, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   
   ////// fake ///////
   TH1F* h_et_all_passing_shower_cut_fake(new TH1F("h_et_all_passing_shower_cut_fake","all passing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_fake(new TH1F("h_et_lead_passing_trail_failing_shower_cut_fake","lead passing trail failing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_fake(new TH1F("h_et_lead_failing_trail_passing_shower_cut_fake","lead failing trail passing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_fake(new TH1F("h_et_all_failing_shower_cut_fake","all failing shower cut, Symmetric Pt;Pt (GeV);Events / GeV", 200, 0., 200.));
   ////// electron //////
   TH1F* h_et_all_passing_shower_cut_electron(new TH1F("h_et_all_passing_shower_cut_electron","all passing shower cut, Symmetric Pt;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_passing_trail_failing_shower_cut_electron(new TH1F("h_et_lead_passing_trail_failing_shower_cut_electron","lead passing trail failing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_lead_failing_trail_passing_shower_cut_electron(new TH1F("h_et_lead_failing_trail_passing_shower_cut_electron","lead failing trail passing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_et_all_failing_shower_cut_electron(new TH1F("h_et_all_failing_shower_cut_electron","all failing shower cut, Symmetric Pt ;Pt (GeV);Events / GeV", 200, 0., 200.));

   ////// MHT histograms ////////
   TH1F* h_mht_candidate(new TH1F("h_mht_candidate","candidate MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_NoCut(new TH1F("h_mht_candidate_NoCut","candidate No Cut MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoIsoCut(new TH1F("h_mht_candidate_PhoIsoCut","candidate with PhoIsoCut MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoNeutralIso(new TH1F("h_mht_candidate_PhoNeutralIso","candidate with Photon and Neutral iso MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoNeutralChargedIso(new TH1F("h_mht_candidate_PhoNeutralChargedIso","candidate with Photon , Charged and Neutral iso MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_candidate_PhoNeutralChargedWorstIso(new TH1F("h_mht_candidate_PhoNeutralChargedWorstIso","candidate with Photon, Charged , Worst and Neutral iso MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_fake(new TH1F("h_mht_fake","fake MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));
   TH1F* h_mht_electron(new TH1F("h_mht_electron","electron MHT, Symmetric Pt ;MHT (GeV);Events / GeV", 500, 0., 500.));

   /////// MHT Delta Phi histograms ///////
   TH1F* h_mht_deltaphi_candidate(new TH1F("h_mht_deltaphi_candidate","#Delta#Phi between MHT and nearest photon, Symmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_NoCut(new TH1F("h_mht_deltaphi_candidate_NoCut","#Delta#Phi MHT and nearest canddiate with no cut, Symmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoIsoCut(new TH1F("h_mht_deltaphi_candidate_PhoIsoCut","#Delta#Phi MHT and candidate with PhoIsoCut, Symmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoNeutralIso(new TH1F("h_mht_deltaphi_candidate_PhoNeutralIso","#Delta#Phi MHT and candidate with Photon and Neutral iso, Symmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoNeutralChargedIso(new TH1F("h_mht_deltaphi_candidate_PhoNeutralChargedIso","#Delta#Phi MHT and candidate with Photon , Charged and Neutral iso, Symmetric Pt ; #Delta#Phi;Events", 80, 0., 4.));
   TH1F* h_mht_deltaphi_candidate_PhoNeutralChargedWorstIso(new TH1F("h_mht_deltaphi_candidate_PhoNeutralChargedWorstIso","#Delta#Phi MHT and candidate with Photon, Charged , Worst and Neutral iso, Symmetric Pt ;#Delta#Phi;Events", 80, 0., 4.));


   // 2D histograms
   /*TH2F* h_leadeta_trailphi_candidate(new TH2F("h_leadeta_trailphi_candidate_HLT","lead #eta vs trail #phi, Symmetric Pt ;lead #eta ;trail #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_leadeta_leadphi_candidate(new TH2F("h_leadeta_leadphi_candidate","lead #eta vs lead #phi, Symmetric Pt ;lead #eta ;lead #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_traileta_leadphi_candidate(new TH2F("h_traileta_leadphi_candidate","trail #eta vs lead #phi, Symmetric Pt;trail #eta ;lead #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_traileta_trailphi_candidate(new TH2F("h_traileta_trailphi_candidate","trail #eta vs trail #phi, Symmetric Pt;trail #eta ;trail #phi ", 80, -2., 2.,160,-4.,4.));
   TH2F* h_leadeta_traileta_candidate(new TH2F("h_leadeta_traileta_candidate","lead #eta vs trail #eta, Symmetric Pt;lead #eta ;trail #eta ", 80, -2., 2.,80,-2.,2.));
   TH2F* h_leadphi_trailphi_candidate(new TH2F("h_leadphi_trailphi_candidate","lead #phi vs trail #phi, Symmetric Pt;lead #phi ;trail #phi ", 160, -4., 4.,160,-4.,4.));
   TH2F* h_leadpt_trailpt_candidate(new TH2F("h_leadpt_trailpt_candidate","lead Pt vs trail Pt , Symmetric Pt;lead Pt;trail Pt ", 300, 0., 300.,300,0.,300.));*/
   // rho histograms
   TH1F* h_rho_ffSample(new TH1F("h_rho_ffSample","ff Rho25, Symmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_gfSample(new TH1F("h_rho_gfSample","gf Rho25, Symmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_fgSample(new TH1F("h_rho_fgSample","fg Rho25, Symmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_eeSample(new TH1F("h_rho_eeSample","ee Rho25, Symmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_eeSample_diEMPtReweighted(new TH1F("h_rho_eeSample_diEMPtReweighted","ee Rho25 diEMPtReweighted to candidate, Symmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   TH1F* h_rho_candidate(new TH1F("h_rho_candidate","candidate Rho25, Symmetric Pt Cut ;Rho25;Events ", 50, 0., 50.));
   // nVtx histograms
   //TH1F* h_nVtx_ffSample(new TH1F("h_nVtx_ffSample","ff nVtx ;nVtx;Events ", 50, 0., 50.));
   //TH1F* h_nVtx_eeSample(new TH1F("h_nVtx_eeSample","ee nVtx ;nVtx;Events ", 50, 0., 50.));
   //TH1F* h_nVtx_candidate(new TH1F("h_nVtx_candidate","candidate nVtx ;nVtx;Events ", 50, 0., 50.));

   // diEMPt histograms
   TH1F* h_diEMPt_eeSample(new TH1F("h_diEMPt_eeSample","ee diEMPt, Symmetric Pt Cut ;diEMPt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_diEMPt_candidate(new TH1F("h_diEMPt_candidate","candidate diEMPt, Symmetric Pt Cut;diEMPt (GeV);Events / GeV", 200, 0., 200.));
   TH1F* h_diEMPt_ffSample(new TH1F("h_diEMPt_ffSample","ff diEMPt, Symmetric Pt Cut ;diEMPt (GeV);Events / GeV", 200, 0., 200.));



   char *histNameee        = new char[50];
   char *histtitleee       = new char[50];
   char *histNameee2       = new char[50];
   char *histtitleee2      = new char[50];
   char *histNameff        = new char[50];
   char *histtitleff       = new char[50];
   char *histNamegf        = new char[50];
   char *histtitlegf       = new char[50];

   for (int d=0;d<1000; ++d) {
     sprintf(histNameee, "DiEMPtReweightedee%d",d+1);
     sprintf(histtitleee,"Reweighted ee MET with DIEMPT %d",d+1);
     h_met_eeSample_DiEMPtRhoReweighted_ErrorPropagatedFromDiEMPt[d]=new TH1F(histNameee,histtitleee,binnum, bins);
     sprintf(histNameee2, "RhoReweightedee%d",d+1);
     sprintf(histtitleee2,"Reweighted ee MET with RHO %d",d+1);
     h_met_eeSample_DiEMPtRhoReweighted_ErrorPropagatedFromRho[d]=new TH1F(histNameee2,histtitleee2,binnum, bins);
     
     sprintf(histNameff, "RhoReweightedff%d",d+1);
     sprintf(histtitleff,"Reweighted ff MET %d",d+1);
     h_met_ffSample_RhoReweighted_ErrorPropagated[d]=new TH1F(histNameff,histtitleff,binnum, bins);
     sprintf(histNamegf, "ReweightedforGF%d",d+1);
     sprintf(histtitlegf,"Reweighted ff MET with purity 0.18 %d",d+1);
     h_met_ffSample_RhoReweightedGF_ErrorPropagated[d]=new TH1F(histNamegf,histtitlegf,binnum, bins);
   }


   char *histNameGF       = new char[50];
   char *histtitleGF      = new char[50];
   char *histNameFG       = new char[50];
   char *histtitleFG      = new char[50];

   for (int d=0;d<1000; ++d) {
     sprintf(histNameGF, "GFRhoReweighted%d",d+1);
     sprintf(histtitleGF,"Reweighted gf MET with RHO %d",d+1);
     h_met_PhoFake_rereco_ErrorFromRho[d]=new TH1F(histNameGF,histtitleGF,200,0.0,200.0);
     sprintf(histNameFG, "FGRhoReweighted%d",d+1);
     sprintf(histtitleFG,"Reweighted fg MET with RHO %d",d+1);
     h_met_FakePho_rereco_ErrorFromRho[d]=new TH1F(histNameFG,histtitleFG,200,0.0,200.0);
   }
   
   // Saving the error weights //
   h_met_eeSample->Sumw2();
   h_met_ffSample->Sumw2();
   h_met_candidate->Sumw2();
   /*h_rho_eeSample->Sumw2();
   h_rho_ffSample->Sumw2();
   h_rho_candidate->Sumw2();
   h_diEMPt_eeSample->Sumw2();
   h_diEMPt_candidate->Sumw2();
   h_diEMPt_ffSample->Sumw2();*/
   //h_met_ffSample_RhoReweighted->Sumw2();
   //h_met_eeSample_diEMPtRhoReweighted->Sumw2();
   h_met_eeSample_PurityReweighted->Sumw2();
   h_met_ffSample_PurityReweighted->Sumw2();
   

   // getting the reweighting //

   TFile *Fnew = new TFile("ElectronRatioSymmetric.root","READ");
   TH1F *candidateeediEMPt = (TH1F*)Fnew->Get("candidateeediEMPt");
   TH1F *candidateeediEMPtweightedrho = (TH1F*)Fnew->Get("candidateeediEMPtweightedrho");
   TH1F *candidateffrho = (TH1F*)Fnew->Get("candidateffrho");
   TH1F *candidategfrho = (TH1F*)Fnew->Get("candidategfrho");
   TH1F *candidatefgrho = (TH1F*)Fnew->Get("candidatefgrho");
   //Fnew->Close();

   TFile *T2 = new TFile("SignalRhoRatioSymmetric.root");
   TH1F *candidatesignalrho = (TH1F*)T2->Get("candidatesignalrho");

   TFile *F2 = new TFile("RandomRatioSymmetric.root","READ");
   TH1F *DiEMPtRatioGaus[1000]; //= (TH1F*)Fnew->Get("candidateeediEMPt");
   TH1F *EERhoRatioGaus[1000];
   TH1F *FFRhoRatioGaus[1000];
   
   for(int g=0; g<1000; ++g){
     //printf("g:%d\n",g);
     char *nameGraph        = new char[40];
     char *nameGraph2       = new char[40];
     char *nameGraph3       = new char[40];

     sprintf(nameGraph,"DiEMPtee%d",g+1);
     DiEMPtRatioGaus[g] = (TH1F*)F2->Get(nameGraph);
     
     sprintf(nameGraph2,"DiEMPtweightedeeRho%d",g+1);
     EERhoRatioGaus[g] = (TH1F*)F2->Get(nameGraph2);
     
     sprintf(nameGraph3,"ffRho%d",g+1);
     FFRhoRatioGaus[g] = (TH1F*)F2->Get(nameGraph3);
   }
   
   
  // F2->Close();


   /*TH1F *h_rho_ffSample = (TH1F*)Fnew->Get("h_rho_ffSample");
   TH1F *h_rho_candidate = (TH1F*)Fnew->Get("h_rho_candidate");
   h_rho_eeSample->Scale(1./h_rho_eeSample->Integral(0,50));
   h_rho_ffSample->Scale(1./h_rho_ffSample->Integral(0,50));
   h_rho_candidate->Scale(1./h_rho_candidate->Integral(0,50));
   TH1F *candidateOvereeRho = new TH1F("candidateOvereeRho", "candidate  to ee Rho ratio", 50, 0., 50.);
   candidateOvereeRho->Sumw2(); 
   candidateOvereeRho->Divide(h_rho_candidate,h_rho_eeSample);
   TH1F *candidateOverffRho = new TH1F("candidateOverffRho", "candidate  to ff Rho ratio", 50, 0., 50.);
   candidateOverffRho->Sumw2();
   candidateOverffRho->Divide(h_rho_candidate,h_rho_ffSample);

   //TFile *Fnew2 = new TFile("Weighted.root","READ");
   TH1F *h_diEMPt_eeSample = (TH1F*)Fnew->Get("h_diEMPt_eeSample");
   TH1F *h_diEMPt_candidate = (TH1F*)Fnew->Get("h_diEMPt_candidate");
   h_diEMPt_eeSample->Scale(1./h_diEMPt_eeSample->Integral(0,200));
   h_diEMPt_candidate->Scale(1./h_diEMPt_candidate->Integral(0,200));
   TH1F *candidateOvereediEMPt = new TH1F("candidateOvereediEMPt", "candidate  to ee diEMPt ratio", 200, 0., 200.);
   candidateOvereediEMPt->Sumw2();
   candidateOvereediEMPt->Divide(h_diEMPt_candidate,h_diEMPt_eeSample);*/
   
    


   TFile *F4 = new TFile("RandomRatioGF.root","READ");
   TH1F *GFRhoRatio[1000];
  
   TH1F *FGRhoRatio[1000];
  
   for(int g=0; g<1000; ++g){
     //printf("g:%d\n",g);
     char *nameGraph4        = new char[50];
     char *nameGraph5        = new char[50];

     sprintf(nameGraph4,"GFRhoRandom%d",g+1);
     GFRhoRatio[g] = (TH1F*)F4->Get(nameGraph4);
     
     sprintf(nameGraph5,"FGRhoRandom%d",g+1);
     FGRhoRatio[g] = (TH1F*)F4->Get(nameGraph5);
   }

   
 

   if (fChain == 0) return;
   FILE *RecoPhoton   = fopen("file_photon_ReReco.txt", "w");
   FILE *RecoFake     = fopen("file_fake_ReReco.txt", "w");
   FILE *RecoAll1     = fopen("file_All_ReReco1.txt", "w");
   FILE *RecoAll2     = fopen("file_All_ReReco2.txt", "w");
   FILE *RecoAll3     = fopen("file_All_ReReco3.txt", "w");
   FILE *RecoAll4     = fopen("file_All_ReReco4.txt", "w");
   FILE *RecoAll5     = fopen("file_All_ReReco5.txt", "w");
   FILE *RecoAll6     = fopen("file_All_ReReco6.txt", "w");
   FILE *RecoAll7     = fopen("file_All_ReReco7.txt", "w");
   FILE *RecoAll8     = fopen("file_All_ReReco8.txt", "w");
   FILE *RecoAll9     = fopen("file_All_ReReco9.txt", "w");
   FILE *RecoAll10    = fopen("file_All_ReReco10.txt", "w");
   FILE *RecoAll11    = fopen("file_All_ReReco11.txt", "w");

   FILE *RecoEle1     = fopen("file_ele_ReReco1.txt", "w");
   FILE *RecoEle2     = fopen("file_ele_ReReco2.txt", "w");
   FILE *RecoEle3     = fopen("file_ele_ReReco3.txt", "w");
   FILE *RecoEle4     = fopen("file_ele_ReReco4.txt", "w");
   FILE *RecoEle5     = fopen("file_ele_ReReco5.txt", "w");
   FILE *RecoEle6     = fopen("file_ele_ReReco6.txt", "w");
   FILE *RecoEle7     = fopen("file_ele_ReReco7.txt", "w");
   FILE *RecoEle8     = fopen("file_ele_ReReco8.txt", "w");
   FILE *RecoEle9     = fopen("file_ele_ReReco9.txt", "w");
   FILE *RecoEle10    = fopen("file_ele_ReReco10.txt", "w");
   FILE *RecoEle11    = fopen("file_ele_ReReco11.txt", "w");
   FILE *RecoEle12    = fopen("file_ele_ReReco12.txt", "w");


   FILE *RecoPhoFake  = fopen("file_phofake_ReReco.txt", "w");
   FILE *RecoFakePho  = fopen("file_fakepho_ReReco.txt", "w");
   
   
   //FILE *RecoRun = fopen("file_ReReco.txt", "w");
   ////// Reading from a file /////////
   ////// Reading from a file /////////
   /*FILE *myfile;
   myfile=fopen("file_Prompt.txt", "r");
   char line[1200];
   int z=0;
   Long64_t ReRecoEventNo[103489]={0};
   Long64_t ReRecoEntry[103489]={0};
   Int_t ReRecoRunNo[103489]={0};
   while(fgets(line, sizeof line, myfile) != NULL){
     fscanf(myfile,"%d\t%lld\t%lld", &ReRecoRunNo[z], &ReRecoEventNo[z], &ReRecoEntry[z]);
     //printf("%d\t%lld\t%lld\n", ReRecoRunNo[z], ReRecoEventNo[z], ReRecoEntry[z]);
     z++;
   }*/

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   printf("Total no of entry: %lld\n",nentries);
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if(jentry%100000==0)printf("Entry processed: %lld\n", jentry);//
      if (ientry < 0) break;
      //if(jentry>600000)break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      /////////  Writing to a file /////////////
      //fprintf(RecoRun,"%lld\t%lld\t%lld\n",runNo,eventNo,jentry);
      
      //fprintf(RecoAll,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);

      if(runNo <= 191235)fprintf(RecoAll1,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 191235 && runNo <= 192014)fprintf(RecoAll2,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 192014 && runNo <= 193572)fprintf(RecoAll3,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 193572 && runNo <= 193768)fprintf(RecoAll4,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 193768 && runNo <= 193963)fprintf(RecoAll5,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 193963 && runNo <= 194158)fprintf(RecoAll6,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 194158 && runNo <= 194353)fprintf(RecoAll7,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 194353 && runNo <= 194548)fprintf(RecoAll8,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 194548 && runNo <= 194743)fprintf(RecoAll9,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 194743 && runNo <= 194938)fprintf(RecoAll10,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      if(runNo > 194938 && runNo <= 195130)fprintf(RecoAll11,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
      

      
      /* /////// Photon, electron and fake Selection /////// */

      int   realphoton = 0;
      int   candidate  = 0;
      int   fakephoton = 0;
      int   fake       = 0;
      int   elephoton  = 0;
      int   electron   = 0;
      int   noCut      = 0;
      int   OnlyPho    = 0;
      int   PhoNeutral = 0;
      int   PhoNeuCharg= 0;
      int   PhoNChargeW= 0;
      vector<int> realphoton_tracker, fakephoton_tracker, elephoton_tracker;
      vector<int> candidate_tracker, fake_tracker, electron_tracker; 
      for(size_t i=0;i < photon_neutraliso->size();++i){
        float absEta = fabs(photon_eta->at(i));
        double effA[4];
        photonEffectiveAreas(absEta, effA);
        float theta = 2*atan(exp(-photon_eta->at(i)));
        float pt = fabs(photon_e->at(i)*sin(theta)); 

        bool neutral     = (photon_neutraliso->at(i) - rho * effA[1] - 0.04 * pt) < 3.5;
        bool charged     = (photon_chargeiso->at(i) - rho * effA[0]) < 2.6;
        bool chargedfake = (photon_chargeiso->at(i) - rho * effA[0]) < 2.6;
        bool chargedlim  = (photon_chargeiso->at(i) - rho * effA[0]) < 15.;
        bool photoniso   = (photon_photoniso->at(i) - rho * effA[2] - 0.005 * pt) < 1.3;
        //printf("absEta: %f and effACh: %f effPh:%f and effN:%f\n", absEta,effA[0],effA[2],effA[1]);
        bool worstiso    = (photon_worstiso->at(i) -rho*effA[3]) < 2.6;
        bool showercut   = photon_showershape->at(i) < 0.011 ; 
        bool shower      = photon_showershape->at(i) < 0.014 ;
        bool pixelcut    = photon_pixelseed->at(i) == 0;
        bool SymPtCut    = pt > 40;
        //if(jentry==535519)printf("yy pixelseed head:%d\n",photon_pixelseed->at(i));
        if(pixelcut && SymPtCut )noCut++;
        if(photoniso && pixelcut && SymPtCut )OnlyPho++;
        if(neutral && photoniso && pixelcut && SymPtCut )PhoNeutral++;
        if(neutral && photoniso && charged && pixelcut && SymPtCut )PhoNeuCharg++;
        if(neutral && photoniso && charged && worstiso && pixelcut && SymPtCut )PhoNChargeW++;
  
        if(neutral && charged && photoniso && showercut && pixelcut && SymPtCut ){
          realphoton++;
          realphoton_tracker.push_back(i);
        } //worstiso not used
        if(neutral && !charged && photoniso && showercut && pixelcut && chargedlim && SymPtCut ){
          fakephoton++;
          fakephoton_tracker.push_back(i);
        }
        if(neutral && charged && photoniso && showercut && !pixelcut && SymPtCut ){
          elephoton++;
          elephoton_tracker.push_back(i);
        }
        if(neutral && charged && photoniso && shower && pixelcut && SymPtCut ){
          candidate++; // worstiso not used.
          candidate_tracker.push_back(i);
        }
        if(neutral && !chargedfake && photoniso && shower && pixelcut && chargedlim && SymPtCut ){
          fake++;
          fake_tracker.push_back(i);
        }
        if(neutral && charged && photoniso && shower && !pixelcut && SymPtCut ){
          electron++;
          electron_tracker.push_back(i);
        }
      }
      if(realphoton >= 2){
        //printf("Event number: %lld\n",eventNo);
        float DR=dRCalc(photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
        //h_dR_candidate->Fill(DR);
        //float DPHI=dPhiCalc(photon_phi->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(1)));
        //h_dPhi_candidate->Fill(DPHI);
        float MASS=findMass(photon_e->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_e->at(realphoton_tracker.at(1)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
        if(DR > 0.6)h_mass_candidate->Fill(MASS);
        if(DR > 0.6 && (MASS<=80 || MASS>=100)){
          /////////  Writing to a file /////////////
          fprintf(RecoPhoton,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          //h_met_rereco_all->Fill(met_et);
          
          float DIEMPT=findDiEMPt(photon_e->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_e->at(realphoton_tracker.at(1)),photon_eta->at(realphoton_tracker.at(1)),photon_phi->at(realphoton_tracker.at(1)));
          h_diEMPt_candidate->Fill(DIEMPT);
          h_rho_candidate->Fill(rho);
          float PxTotal = 0;
          float PyTotal = 0;
          for(int k = 0; k < realphoton; ++k){
            float theta = 2*atan(exp(-photon_eta->at(realphoton_tracker.at(k))));
            float PX    = photon_e->at(realphoton_tracker.at(k))*sin(theta)*cos(photon_phi->at(realphoton_tracker.at(k))); 
            float PY    = photon_e->at(realphoton_tracker.at(k))*sin(theta)*sin(photon_phi->at(realphoton_tracker.at(k)));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          float redX    = met_X+PxTotal;
          float redY    = met_Y+PyTotal;
          float redMET  = sqrt(redX*redX+redY*redY);
          h_lead_photon_pt_rereco->Fill(photon_e->at(realphoton_tracker.at(0))*2*atan(exp(-photon_eta->at(realphoton_tracker.at(0)))));
          h_trail_photon_pt_rereco->Fill(photon_e->at(realphoton_tracker.at(1))*2*atan(exp(-photon_eta->at(realphoton_tracker.at(1)))));
          h_lead_photon_phi_rereco->Fill(photon_phi->at(realphoton_tracker.at(0)));
          h_trail_photon_phi_rereco->Fill(photon_phi->at(realphoton_tracker.at(1)));
          h_lead_photon_eta_rereco->Fill(photon_eta->at(realphoton_tracker.at(0)));
          h_trail_photon_eta_rereco->Fill(photon_eta->at(realphoton_tracker.at(1)));
          h_diEMPt_rereco->Fill(DIEMPT);
          h_met_rereco->Fill(met_et);
          h_met_candidate->Fill(met_et);
          int BinC=candidatesignalrho->FindBin(rho);
          float WeightC=candidatesignalrho->GetBinContent(BinC);
          if(met_et > 60 && met_et <= 100)h_met_candidate_MET60_100->Fill(met_et, WeightC);
          if(met_et > 100)h_met_candidate_MET100_Up->Fill(met_et, WeightC); 
          h_metX_rereco->Fill(met_X);
          h_metY_rereco->Fill(met_Y);
          h_reduced_met_rereco->Fill(redMET);
          for(size_t s=0; s<jet_pt->size(); ++s){
            float dr=0;
            int count = 0;
            for(int j=0; j<realphoton; ++j){
              dr=dRCalc(jet_eta->at(s),jet_phi->at(s),photon_eta->at(realphoton_tracker.at(j)),photon_phi->at(realphoton_tracker.at(j)));
              if(dr<0.6)count++;
            }
            if(count == 0){
              h_jetpt_rereco->Fill(jet_pt->at(s));
              h_jeteta_rereco->Fill(jet_eta->at(s));
              h_jetphi_rereco->Fill(jet_phi->at(s));
           
            }
          }
        }
        
      }
      
      if(realphoton>=1 && fakephoton>=1){
        float DR=dRCalc(photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_eta->at(fakephoton_tracker.at(0)),photon_phi->at(fakephoton_tracker.at(0)));
        if(DR > 0.6){
          
          if(realphoton_tracker.at(0)<fakephoton_tracker.at(0)){
            
            fprintf(RecoPhoFake,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
            h_met_PhoFake_rereco->Fill(met_et);
            h_rho_gfSample->Fill(rho);
            //// reweighting ///
            int Bin1=candidategfrho->FindBin(rho);
            float Weight1=candidategfrho->GetBinContent(Bin1);
            h_met_PhoFake_rereco_RhoReweighted->Fill(met_et, Weight1);


            for( int s=0 ; s<1000; ++s){
              float Weightprop = GFRhoRatio[s]->GetBinContent(Bin1);
              h_met_PhoFake_rereco_ErrorFromRho[s]->Fill(met_et, Weightprop);
            }
         
            //h_met_PhoFake_PurityReweighted->Fill(met_et, 0.260527);
          }
          if(realphoton_tracker.at(0)>fakephoton_tracker.at(0)){
            //printf("phofake event second:%lld\n",jentry);
            fprintf(RecoFakePho,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
            h_met_FakePho_rereco->Fill(met_et);
            h_rho_fgSample->Fill(rho);
            //// reweighting ///
            int Bin2=candidatefgrho->FindBin(rho);
            float Weight2=candidatefgrho->GetBinContent(Bin2);
            h_met_FakePho_rereco_RhoReweighted->Fill(met_et, Weight2);

            for( int t=0 ; t<1000; ++t){
              float WeightFromRho = FGRhoRatio[t]->GetBinContent(Bin2);
              h_met_FakePho_rereco_ErrorFromRho[t]->Fill(met_et, WeightFromRho);
            }
  
            //h_met_FakePho_PurityReweighted->Fill(met_et, 0.09079);
          }
        } 
          
      }

      // getting EWK background
      
      if(realphoton>=1 && elephoton>=1){
        float DR=dRCalc(photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_eta->at(elephoton_tracker.at(0)),photon_phi->at(elephoton_tracker.at(0)));
        float MASS=findMass(photon_e->at(realphoton_tracker.at(0)),photon_eta->at(realphoton_tracker.at(0)),photon_phi->at(realphoton_tracker.at(0)),photon_e->at(elephoton_tracker.at(0)),photon_eta->at(elephoton_tracker.at(0)),photon_phi->at(elephoton_tracker.at(0)));
        
        if(photon_pixelseed->at(realphoton_tracker.at(0))!=0 && photon_pixelseed->at(elephoton_tracker.at(0))!=0){
          h_mass_eeSample_Insideey->Fill(MASS);
          printf("ee event: %lld\n",jentry);
        }
        if(photon_pixelseed->at(realphoton_tracker.at(0))==0 && photon_pixelseed->at(elephoton_tracker.at(0))==0){
          h_mass_yySample_Insideey->Fill(MASS);
          printf("yy event: %lld\n",jentry);
        }
        if(DR > 0.6 && (MASS<=80 || MASS>=100)){
          h_mass_eySample->Fill(MASS);
          float OpAng = fabs(photon_phi->at(realphoton_tracker.at(0))-photon_phi->at(elephoton_tracker.at(0)));
          float OpAngEta = fabs(photon_eta->at(realphoton_tracker.at(0))-photon_eta->at(elephoton_tracker.at(0)));
          float OpAngDr = sqrt(OpAng*OpAng+OpAngEta*OpAngEta);
          h_eySample_OpeAng_Eta->Fill(OpAngEta);
          h_eySample_DeltaR->Fill(OpAngDr);
          if(realphoton_tracker.at(0)<elephoton_tracker.at(0)){
            h_eySample_LeadPhi->Fill(photon_phi->at(realphoton_tracker.at(0)));
            h_eySample_LeadEta->Fill(photon_eta->at(realphoton_tracker.at(0)));
            h_eySample_TrailPhi->Fill(photon_phi->at(elephoton_tracker.at(0)));
            h_eySample_TrailEta->Fill(photon_eta->at(elephoton_tracker.at(0)));
          }
          else{
            h_eySample_LeadPhi->Fill(photon_phi->at(elephoton_tracker.at(0)));
            h_eySample_LeadEta->Fill(photon_eta->at(elephoton_tracker.at(0)));
            h_eySample_TrailPhi->Fill(photon_phi->at(realphoton_tracker.at(0)));
            h_eySample_TrailEta->Fill(photon_eta->at(realphoton_tracker.at(0)));
          }
          if(OpAng > TMath::Pi())h_eySample_OpeAng->Fill(2*TMath::Pi() - OpAng);
          else h_eySample_OpeAng->Fill(OpAng);
 
          float PxTotal = 0;
          float PyTotal = 0;
          for(int j=0; j<elephoton; ++j){
            float theta = 2*atan(exp(-photon_eta->at(elephoton_tracker.at(j))));
            float PX    = photon_e->at(elephoton_tracker.at(j))*sin(theta)*cos(photon_phi->at(elephoton_tracker.at(j))); 
            float PY    = photon_e->at(elephoton_tracker.at(j))*sin(theta)*sin(photon_phi->at(elephoton_tracker.at(j)));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(int j=0; j<realphoton; ++j){
            float theta = 2*atan(exp(-photon_eta->at(realphoton_tracker.at(j))));
            float PX    = photon_e->at(realphoton_tracker.at(j))*sin(theta)*cos(photon_phi->at(realphoton_tracker.at(j))); 
            float PY    = photon_e->at(realphoton_tracker.at(j))*sin(theta)*sin(photon_phi->at(realphoton_tracker.at(j)));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          float redX    = met_X+PxTotal;
          float redY    = met_Y+PyTotal;
          float redMET  = sqrt(redX*redX+redY*redY);
          h_reduced_met_eySample->Fill(redMET);


          float pte(0);
          // selecting pt_e
          for(size_t h=0; h < photon_eta->size();++h){
            if(photon_pixelseed->at(h)!=0){
              float thetae = 2*atan(exp(-photon_eta->at(h)));
              pte = fabs(photon_e->at(h)*sin(thetae));
              break;  
            }
          }
          if(pte!=0){
            float fetogamma = 1-(1 - 0.00700)*(1-TMath::Power((pte/2.9+1),-2.4))*(1-0.23*exp(-0.277*tracks_n))*(1-0.000566*vertices_n); // from AN-13-240
            float scalefactor = fetogamma/(1 - fetogamma);
            float highscale   = (fetogamma+0.11*fetogamma)/(1-(fetogamma+0.11*fetogamma));
            float lowscale    = (fetogamma-0.11*fetogamma)/(1-(fetogamma-0.11*fetogamma));
            //printf("scale factor: %f\n", scalefactor); 
            h_met_eySample->Fill(met_et);
            if(met_et > 60 && met_et <= 100)h_met_eySample_MET60_100->Fill(met_et);
            if(met_et > 100)h_met_eySample_MET100_Up->Fill(met_et);
            h_met_eySample_Scaled->Fill(met_et,scalefactor);
            h_met_eySample_Scaled_highf->Fill(met_et,highscale);
            h_met_eySample_Scaled_lowf->Fill(met_et,lowscale);
          }
        }        
      }




      if(fakephoton >= 2){
        float DR=dRCalc(photon_eta->at(fakephoton_tracker.at(0)),photon_phi->at(fakephoton_tracker.at(0)),photon_eta->at(fakephoton_tracker.at(1)),photon_phi->at(fakephoton_tracker.at(1)));
        if(DR > 0.6){
          /////////  Writing to a file /////////////
          fprintf(RecoFake,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          //
          h_met_ffSample->Fill(met_et);
          if(met_et > 60 && met_et <= 100)h_met_ffSample_MET60_100->Fill(met_et);
          if(met_et > 100)h_met_ffSample_MET100_Up->Fill(met_et);
          
          float PxTotal = 0;
          float PyTotal = 0;
          for(int j = 0; j < fakephoton; ++j){
            float theta = 2*atan(exp(-photon_eta->at(fakephoton_tracker.at(j))));
            float PX    = photon_e->at(fakephoton_tracker.at(j))*sin(theta)*cos(photon_phi->at(fakephoton_tracker.at(j))); 
            float PY    = photon_e->at(fakephoton_tracker.at(j))*sin(theta)*sin(photon_phi->at(fakephoton_tracker.at(j)));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          float redX    = met_X+PxTotal;
          float redY    = met_Y+PyTotal;
          float redMET  = sqrt(redX*redX+redY*redY);
          h_reduced_met_ffSample->Fill(redMET);
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(int j=0; j< fakephoton; ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(fakephoton_tracker.at(j)),photon_phi->at(fakephoton_tracker.at(j)));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_fake->Fill(MHT);
          }
          float DIEMPT=findDiEMPt(photon_e->at(fakephoton_tracker.at(0)),photon_eta->at(fakephoton_tracker.at(0)),photon_phi->at(fakephoton_tracker.at(0)),photon_e->at(fakephoton_tracker.at(1)),photon_eta->at(fakephoton_tracker.at(1)),photon_phi->at(fakephoton_tracker.at(1)));
          h_diEMPt_ffSample->Fill(DIEMPT);
          h_rho_ffSample->Fill(rho);
          int Bin11=candidateffrho->FindBin(rho);
          float Weight11=candidateffrho->GetBinContent(Bin11);
          float Binerror=candidateffrho->GetBinError(Bin11);
          
         
          h_met_ffSample_RhoReweighted->Fill(met_et, Weight11);
          
          //h_met_ffSample_PurityReweighted->Fill(met_et, Weight*0.48459);
          //h_met_ffSample_PurityReweighted->Fill(met_et, Weight11*0.182506);     
          for( int t=0 ; t<1000; ++t){
            
            float WeightFromRho = FFRhoRatioGaus[t]->GetBinContent(Bin11);
            //printf("%f\n",Weightprop);
            h_met_ffSample_RhoReweighted_ErrorPropagated[t]->Fill(met_et, WeightFromRho);
            //h_met_ffSample_RhoReweightedGF_ErrorPropagated[t]->Fill(met_et, 0.182506*WeightFromRho);
            
          }
          
        }
      }
      if(elephoton >= 2){
        float DR=dRCalc(photon_eta->at(elephoton_tracker.at(0)),photon_phi->at(elephoton_tracker.at(0)),photon_eta->at(elephoton_tracker.at(1)),photon_phi->at(elephoton_tracker.at(1)));
        float MASS=findMass(photon_e->at(elephoton_tracker.at(0)),photon_eta->at(elephoton_tracker.at(0)),photon_phi->at(elephoton_tracker.at(0)),photon_e->at(elephoton_tracker.at(1)),photon_eta->at(elephoton_tracker.at(1)),photon_phi->at(elephoton_tracker.at(1)));
        h_mass_eeSample->Fill(MASS);
        if(DR > 0.6 && MASS>=75 && MASS<=105 ){
          float OpAng = fabs(photon_phi->at(elephoton_tracker.at(0))-photon_phi->at(elephoton_tracker.at(1)));
          float OpAngEta = fabs(photon_eta->at(elephoton_tracker.at(0))-photon_eta->at(elephoton_tracker.at(1)));
          float OpAngDr = sqrt(OpAng*OpAng+OpAngEta*OpAngEta);
          h_eeSample_OpeAng_Eta->Fill(OpAngEta);
          h_eeSample_DeltaR->Fill(OpAngDr);
          if(OpAng > TMath::Pi())h_eeSample_OpeAng->Fill(2*TMath::Pi() - OpAng);
          else h_eeSample_OpeAng->Fill(OpAng);


          h_eeSample_LeadPhi->Fill(photon_phi->at(elephoton_tracker.at(0)));
          h_eeSample_LeadEta->Fill(photon_eta->at(elephoton_tracker.at(0)));
          h_eeSample_TrailPhi->Fill(photon_phi->at(elephoton_tracker.at(1)));
          h_eeSample_TrailEta->Fill(photon_eta->at(elephoton_tracker.at(1)));

          /////////  Writing to a file /////////////
          if(runNo <= 193572)fprintf(RecoEle1,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 193572 && runNo <= 195130)fprintf(RecoEle2,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 195130 && runNo <= 196688)fprintf(RecoEle3,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 196688 && runNo <= 198246)fprintf(RecoEle4,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 198246 && runNo <= 199804)fprintf(RecoEle5,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 199804 && runNo <= 201362)fprintf(RecoEle6,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 201362 && runNo <= 202920)fprintf(RecoEle7,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 202920 && runNo <= 204478)fprintf(RecoEle8,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 204478 && runNo <= 206036)fprintf(RecoEle9,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 206036 && runNo <= 206815)fprintf(RecoEle10,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 206815 && runNo <= 207594)fprintf(RecoEle11,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          if(runNo > 207594 && runNo <= 209152)fprintf(RecoEle12,"%d\t%lld\t%lld\n",runNo,eventNo,jentry);
          //printf("elephoton Entry: %lld\n", jentry);
          h_met_eeSample->Fill(met_et);
          
          float DIEMPT=findDiEMPt(photon_e->at(elephoton_tracker.at(0)),photon_eta->at(elephoton_tracker.at(0)),photon_phi->at(elephoton_tracker.at(0)),photon_e->at(elephoton_tracker.at(1)),photon_eta->at(elephoton_tracker.at(1)),photon_phi->at(elephoton_tracker.at(1)));
          float PxTotal = 0;
          float PyTotal = 0;
          for(int j=0; j<elephoton; ++j){
            float theta = 2*atan(exp(-photon_eta->at(elephoton_tracker.at(j))));
            float PX    = photon_e->at(elephoton_tracker.at(j))*sin(theta)*cos(photon_phi->at(elephoton_tracker.at(j))); 
            float PY    = photon_e->at(elephoton_tracker.at(j))*sin(theta)*sin(photon_phi->at(elephoton_tracker.at(j)));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          float redX    = met_X+PxTotal;
          float redY    = met_Y+PyTotal;
          float redMET  = sqrt(redX*redX+redY*redY);
          h_reduced_met_eeSample->Fill(redMET);
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(int j=0; j<elephoton; ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(elephoton_tracker.at(j)),photon_phi->at(elephoton_tracker.at(j)));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            } 
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_electron->Fill(MHT);
          }
          int Bin=candidateeediEMPt->FindBin(DIEMPT);
          float Weight=candidateeediEMPt->GetBinContent(Bin);
          int Bin2=candidateeediEMPtweightedrho->FindBin(rho);
          float Weight2=candidateeediEMPtweightedrho->GetBinContent(Bin2);
          h_met_eeSample_diEMPtRhoReweighted->Fill(met_et, Weight*Weight2);
          if(met_et > 60 && met_et <= 100)h_met_eeSample_MET60_100->Fill(met_et);
          if(met_et > 100)h_met_eeSample_MET100_Up->Fill(met_et);
          //h_met_eeSample_PurityReweighted->Fill(met_et, Weight*Weight2*0.466175);

          h_diEMPt_eeSample->Fill(DIEMPT);
          h_rho_eeSample->Fill(rho);
          h_rho_eeSample_diEMPtReweighted->Fill(rho, Weight);
          for( int s=0 ; s<1000; ++s){
            float Weightprop = DiEMPtRatioGaus[s]->GetBinContent(Bin);
            h_met_eeSample_DiEMPtRhoReweighted_ErrorPropagatedFromDiEMPt[s]->Fill(met_et, Weightprop);
          }
          for( int t=0 ; t<1000; ++t){
            float WeightFromRho = EERhoRatioGaus[t]->GetBinContent(Bin2);
            h_met_eeSample_DiEMPtRhoReweighted_ErrorPropagatedFromRho[t]->Fill(met_et, WeightFromRho);
          }
        }
      }
      // calculating purity from matrix method
      if(candidate >= 2){
        float DR=dRCalc(photon_eta->at(candidate_tracker.at(0)),photon_phi->at(candidate_tracker.at(0)),photon_eta->at(candidate_tracker.at(1)),photon_phi->at(candidate_tracker.at(1)));
        float MASS=findMass(photon_e->at(candidate_tracker.at(0)),photon_eta->at(candidate_tracker.at(0)),photon_phi->at(candidate_tracker.at(0)),photon_e->at(candidate_tracker.at(1)),photon_eta->at(candidate_tracker.at(1)),photon_phi->at(candidate_tracker.at(1)));
        if(DR>0.6){
          if(photon_showershape->at(candidate_tracker.at(0))<0.011 && photon_showershape->at(candidate_tracker.at(1))<0.011 && (MASS < 80 || MASS > 100)){
            h_et_all_passing_shower_cut_candidate->Fill(photon_e->at(candidate_tracker.at(0))); 
            if(met_et <=20)h_et_all_passing_shower_cut_candidate_20->Fill(photon_e->at(candidate_tracker.at(0)));
            if(20 < met_et && met_et <=40)h_et_all_passing_shower_cut_candidate_40->Fill(photon_e->at(candidate_tracker.at(0)));
            if(40 < met_et && met_et <=60)h_et_all_passing_shower_cut_candidate_60->Fill(photon_e->at(candidate_tracker.at(0)));
            if(60 < met_et && met_et <=80)h_et_all_passing_shower_cut_candidate_80->Fill(photon_e->at(candidate_tracker.at(0)));
            if(80 < met_et && met_et <=100)h_et_all_passing_shower_cut_candidate_100->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 100)h_et_all_passing_shower_cut_candidate_100Up->Fill(photon_e->at(candidate_tracker.at(0)));
            
            float PxTotal = 0;
            float PyTotal = 0;
            for(int j=0; j<candidate; ++j){
              float theta = 2*atan(exp(-photon_eta->at(candidate_tracker.at(j))));
              float PX    = photon_e->at(candidate_tracker.at(j))*sin(theta)*cos(photon_phi->at(candidate_tracker.at(j))); 
              float PY    = photon_e->at(candidate_tracker.at(j))*sin(theta)*sin(photon_phi->at(candidate_tracker.at(j)));
              PxTotal    += PX;
              //printf("PxTotal: %f and PX: %f for entry: %lld\n", PxTotal, PX, jentry);
              PyTotal    += PY;
              //printf("PyTotal: %f and PY: %f for entry: %lld\n", PyTotal, PY, jentry);
            }
            for(size_t k=0; k<jet_pt->size(); ++k){
              float dr=0;
              int counter = 0;
              for(int j=0; j<candidate; ++j){
                dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(candidate_tracker.at(j)),photon_phi->at(candidate_tracker.at(j)));
                if(dr<0.6)counter++;
              }
              if(counter == 0){
                float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
                float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
                PxTotal    += PXj;
                //printf("PxTotal: %f and PXj: %f for entry: %lld\n", PxTotal, PXj, jentry);
                PyTotal    += PYj;
                //printf("PyTotal: %f and PYj: %f for entry: %lld\n", PyTotal, PYj, jentry);
              }
            }
            if(PxTotal!=0 && PyTotal!=0){
              float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
              //printf("final PxTotal: %f PyTotal: %f for entry: %lld and MHT: %f\n",PxTotal, PyTotal , jentry, MHT);
              h_mht_candidate->Fill(MHT); 
              float MHTPhi=atan2(-PyTotal,-PxTotal);
              float deltaphi[10]={0};
              for(int i=0; i< candidate;++i){
                deltaphi[i]=fabs(MHTPhi-photon_phi->at(candidate_tracker.at(i)));
                if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
              }
              float mindphi=deltaphi[0];
              for(int i=0;i< candidate;++i){
                if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
              }
              h_mht_deltaphi_candidate->Fill(mindphi);
               
            }
            float DIEMPT=findDiEMPt(photon_e->at(candidate_tracker.at(0)),photon_eta->at(candidate_tracker.at(0)),photon_phi->at(candidate_tracker.at(0)),photon_e->at(candidate_tracker.at(1)),photon_eta->at(candidate_tracker.at(1)),photon_phi->at(candidate_tracker.at(1)));
            // 
          }
          if(photon_showershape->at(candidate_tracker.at(0))<0.011 && photon_showershape->at(candidate_tracker.at(1))>0.011 && (MASS < 80 || MASS > 100)){
            h_et_lead_passing_trail_failing_shower_cut_candidate->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et <=20)h_et_lead_passing_trail_failing_shower_cut_candidate_20->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 20 && met_et <=40)h_et_lead_passing_trail_failing_shower_cut_candidate_40->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 40 && met_et <=60)h_et_lead_passing_trail_failing_shower_cut_candidate_60->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 60 && met_et <=80)h_et_lead_passing_trail_failing_shower_cut_candidate_80->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 80 && met_et <=100)h_et_lead_passing_trail_failing_shower_cut_candidate_100->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 100)h_et_lead_passing_trail_failing_shower_cut_candidate_100Up->Fill(photon_e->at(candidate_tracker.at(0)));
          }
          if(photon_showershape->at(candidate_tracker.at(0))>0.011 && photon_showershape->at(candidate_tracker.at(1))<0.011 && (MASS < 80 || MASS > 100)){
            h_et_lead_failing_trail_passing_shower_cut_candidate->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et <=20)h_et_lead_failing_trail_passing_shower_cut_candidate_20->Fill(photon_e->at(candidate_tracker.at(0)));
            if(20 < met_et && met_et <=40)h_et_lead_failing_trail_passing_shower_cut_candidate_40->Fill(photon_e->at(candidate_tracker.at(0)));
            if(40 < met_et && met_et <=60)h_et_lead_failing_trail_passing_shower_cut_candidate_60->Fill(photon_e->at(candidate_tracker.at(0)));
            if(60 < met_et && met_et <=80)h_et_lead_failing_trail_passing_shower_cut_candidate_80->Fill(photon_e->at(candidate_tracker.at(0)));
            if(80 < met_et && met_et <=100)h_et_lead_failing_trail_passing_shower_cut_candidate_100->Fill(photon_e->at(candidate_tracker.at(0)));
            if( met_et > 100)h_et_lead_failing_trail_passing_shower_cut_candidate_100Up->Fill(photon_e->at(candidate_tracker.at(0)));
          }
          if(photon_showershape->at(candidate_tracker.at(0))>0.011 && photon_showershape->at(candidate_tracker.at(1))>0.011 && (MASS < 80 || MASS > 100)){
            h_et_all_failing_shower_cut_candidate->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et <= 20)h_et_all_failing_shower_cut_candidate_20->Fill(photon_e->at(candidate_tracker.at(0)));
            if(20 < met_et && met_et <= 40)h_et_all_failing_shower_cut_candidate_40->Fill(photon_e->at(candidate_tracker.at(0)));
            if(40 < met_et && met_et <= 60)h_et_all_failing_shower_cut_candidate_60->Fill(photon_e->at(candidate_tracker.at(0)));
            if(60 < met_et && met_et <= 80)h_et_all_failing_shower_cut_candidate_80->Fill(photon_e->at(candidate_tracker.at(0)));
            if(80 < met_et && met_et <= 100)h_et_all_failing_shower_cut_candidate_100->Fill(photon_e->at(candidate_tracker.at(0)));
            if(met_et > 100)h_et_all_failing_shower_cut_candidate_100Up->Fill(photon_e->at(candidate_tracker.at(0)));
          }
        }
      }
      if(fake >= 2){
        float DR=dRCalc(photon_eta->at(fake_tracker.at(0)),photon_phi->at(fake_tracker.at(0)),photon_eta->at(fake_tracker.at(1)),photon_phi->at(fake_tracker.at(1)));
        if(DR>0.6){
          if(photon_showershape->at(fake_tracker.at(0))<0.011 && photon_showershape->at(fake_tracker.at(1))<0.011)h_et_all_passing_shower_cut_fake->Fill(photon_e->at(fake_tracker.at(0)));
          if(photon_showershape->at(fake_tracker.at(0))<0.011 && photon_showershape->at(fake_tracker.at(1))>0.011)h_et_lead_passing_trail_failing_shower_cut_fake->Fill(photon_e->at(fake_tracker.at(0)));
          if(photon_showershape->at(fake_tracker.at(0))>0.011 && photon_showershape->at(fake_tracker.at(1))<0.011)h_et_lead_failing_trail_passing_shower_cut_fake->Fill(photon_e->at(fake_tracker.at(0)));
          if(photon_showershape->at(fake_tracker.at(0))>0.011 && photon_showershape->at(fake_tracker.at(1))>0.011)h_et_all_failing_shower_cut_fake->Fill(photon_e->at(fake_tracker.at(0)));
        }
      }

      if(electron >= 2){
        float DR=dRCalc(photon_eta->at(electron_tracker.at(0)),photon_phi->at(electron_tracker.at(0)),photon_eta->at(electron_tracker.at(1)),photon_phi->at(electron_tracker.at(1)));
        float MASS=findMass(photon_e->at(electron_tracker.at(0)),photon_eta->at(electron_tracker.at(0)),photon_phi->at(electron_tracker.at(0)),photon_e->at(electron_tracker.at(1)),photon_eta->at(electron_tracker.at(1)),photon_phi->at(electron_tracker.at(1)));
        
        if(DR>0.6 && MASS>=75 && MASS<=105){
          if(photon_showershape->at(electron_tracker.at(0))<0.011 && photon_showershape->at(electron_tracker.at(1))<0.011)h_et_all_passing_shower_cut_electron->Fill(photon_e->at(electron_tracker.at(0)));
          if(photon_showershape->at(electron_tracker.at(0))<0.011 && photon_showershape->at(electron_tracker.at(1))>0.011)h_et_lead_passing_trail_failing_shower_cut_electron->Fill(photon_e->at(electron_tracker.at(0)));
          if(photon_showershape->at(electron_tracker.at(0))>0.011 && photon_showershape->at(electron_tracker.at(1))<0.011)h_et_lead_failing_trail_passing_shower_cut_electron->Fill(photon_e->at(electron_tracker.at(0)));
          if(photon_showershape->at(electron_tracker.at(0))>0.011 && photon_showershape->at(electron_tracker.at(1))>0.011)h_et_all_failing_shower_cut_electron->Fill(photon_e->at(electron_tracker.at(0)));
        }
      }

      if(noCut >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            //if(counter > 0)printf("counter: %d for entry:%lld\n ",counter, jentry);
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_NoCut->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(size_t i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(size_t i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_NoCut->Fill(mindphi);  
          }
          h_met_candidate_NoCut->Fill(met_et);
          
        }
      }
      if(OnlyPho >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoIsoCut->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(size_t i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(size_t i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoIsoCut->Fill(mindphi);  
          }
          h_met_candidate_PhoIsoCut->Fill(met_et);
        }
      }
      if(PhoNeutral >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoNeutralIso->Fill(MHT); 
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(size_t i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(size_t i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoNeutralIso->Fill(mindphi); 
          }
          h_met_candidate_PhoNeutralIso->Fill(met_et);
        }
      }
      if(PhoNeuCharg >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj; 
            }
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoNeutralChargedIso->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(size_t i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(size_t i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoNeutralChargedIso->Fill(mindphi);  
          }
          h_met_candidate_PhoNeutralChargedIso->Fill(met_et);
        }
      }
      if(PhoNChargeW >= 2){
        float DR=dRCalc(photon_eta->at(0),photon_phi->at(0),photon_eta->at(1),photon_phi->at(1));
        if(DR > 0.6){
          float PxTotal = 0;
          float PyTotal = 0;
          for(size_t j=0; j<photon_e->size(); ++j){
            float theta = 2*atan(exp(-photon_eta->at(j)));
            float PX    = photon_e->at(j)*sin(theta)*cos(photon_phi->at(j)); 
            float PY    = photon_e->at(j)*sin(theta)*sin(photon_phi->at(j));
            PxTotal    += PX;
            PyTotal    += PY;
          }
          for(size_t k=0; k<jet_pt->size(); ++k){
            float dr=0;
            int counter = 0;
            for(size_t j=0; j<photon_e->size(); ++j){
              dr=dRCalc(jet_eta->at(k),jet_phi->at(k),photon_eta->at(j),photon_phi->at(j));
              if(dr<0.6)counter++;
            }
            if(counter == 0){
              float PXj   = jet_pt->at(k)*cos(jet_phi->at(k));
              float PYj   = jet_pt->at(k)*sin(jet_phi->at(k));
              PxTotal    += PXj;
              PyTotal    += PYj;
            } 
          }
          if(PxTotal!=0 && PyTotal!=0){
            float MHT=sqrt(PxTotal*PxTotal+PyTotal*PyTotal);
            h_mht_candidate_PhoNeutralChargedWorstIso->Fill(MHT);
            float MHTPhi=atan2(-PyTotal,-PxTotal);
            float deltaphi[10]={0};
            for(size_t i=0; i< photon_e->size();++i){
              deltaphi[i]=fabs(MHTPhi-photon_phi->at(i));
              if(deltaphi[i] > TMath::Pi() ) deltaphi[i] = 2.*TMath::Pi()-deltaphi[i]; 
            }
            float mindphi=deltaphi[0];
            for(size_t i=0;i<photon_e->size();++i){
              if(deltaphi[i]<mindphi)mindphi=deltaphi[i];
            }
            h_mht_deltaphi_candidate_PhoNeutralChargedWorstIso->Fill(mindphi);  
          }
          h_met_candidate_PhoNeutralChargedWorstIso->Fill(met_et);
        }
      }

      // if (Cut(ientry) < 0) continue;*/
     
      
   }
   fclose(RecoAll1);
   fclose(RecoAll2);
   fclose(RecoAll3);
   fclose(RecoAll4);
   fclose(RecoAll5);
   fclose(RecoAll6);
   fclose(RecoAll7);
   fclose(RecoAll8);
   fclose(RecoAll9);
   fclose(RecoAll10);
   fclose(RecoAll11);
   fclose(RecoPhoton);
   fclose(RecoFake);
   fclose(RecoEle1);
   fclose(RecoEle2);
   fclose(RecoEle3);
   fclose(RecoEle4);
   fclose(RecoEle5);
   fclose(RecoEle6);
   fclose(RecoEle7);
   fclose(RecoEle8);
   fclose(RecoEle9);
   fclose(RecoEle10);
   fclose(RecoEle11);
   fclose(RecoEle12);
   
   f1->cd();
   f1->Write();
   //f1->Close();
   
}
