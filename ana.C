// Original Author:  Dongwook Jang
// $Id: ana.C,v 1.12 2013/05/07 20:29:39 yiiyama Exp $
#include "TString.h"
#include "TChain.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>


void ana(TString outputName="analysis_SUSYSignal_1350_375"){
   TChain chain("susyTree");//susyTree

  // REMOVE THE LINE BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
  gSystem->Load("libCondFormatsJetMETObjects.so");

  gSystem->Load("libSusyEvent.so");

  gSystem->AddIncludePath("-I" + TString(gSystem->Getenv("CMSSW_RELEASE_BASE")) + "/src");

  // Analysis macro
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");
  //gSystem->Load("SusyEventAnalyzer_cc.so");

  // chain of inputs
 
  //chain.Add("uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/susyEvents_100_1_Zsk.root");
  /*chain.Add("dcap::///pnfs/cms/WAX/resilient/bfrancis/SusyNtuples/cms538v1/DiPhotonJets_8TeV-madgraph-tarball-v2/susyEvents_182_1_A1U.root"); // MonteCarlo DiPhoton Jet cms538v1 dcache
  chain.Add("dcap::///pnfs/cms/WAX/resilient/bfrancis/SusyNtuples/cms538v1/DiPhotonJets_8TeV-madgraph-tarball-v2/susyEvents_100_1_x7J.root");*/ 

  /*chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v0p1/Summer12_DR53X/GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/susyEvents_118_2_8q7.root");// MonteCarlo GJet cms538v0p1
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v0p1/Summer12_DR53X/GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6/susyEvents_333_2_D8u.root");*/

  /*chain.Add("/eos/uscms/store/user/bfrancis/SusyNtuples/cms533v1/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/susyEvents_320_1_cmG.root"); // MonteCarlo QCD DiEMEnriched cms533v1
  chain.Add("/eos/uscms/store/user/bfrancis/SusyNtuples/cms533v1/QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/susyEvents_1188_1_gbb.root");*/


   //chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms533v1/Run2012D-PromptReco-v1/DoublePhoton/Runs207920-209151/susyEvents_222_2_smE.root"); // data Tag cms533v1

  /*chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1000_1_Xkh.root");// data Tag cms538v1
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1001_1_kKI.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1002_1_Kye.root");*/

  
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_156.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_161.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_187.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_198.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_214.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_227.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_23.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_242.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_264.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_273.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_276.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_292.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_296.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_300.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_301.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_314.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_345.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_381.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_396.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_405.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_414.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_424.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_436.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_445.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_459.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_473.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_486.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_522.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_546.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_552.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_584.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_586.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_587.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_595.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_599.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_613.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_620.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_622.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_634.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_637.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_65.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_665.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_668.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_673.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_686.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_689.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_703.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_709.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_71.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_733.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_752.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_766.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_784.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_788.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_803.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_810.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_819.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_831.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_88.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-1350to1500_mLSP-25to1475_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_1350_75_89.root");


  
  
  /*chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1003_3_KvH.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1004_1_C9b.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1005_1_7SL.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1006_1_MxE.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1007_2_VhO.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1008_2_Dth.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_1009_2_ONU.root");*/

  //chain.Add("dcap://cmsdca1.fnal.gov:24140/pnfs/fnal.gov/usr/cms/WAX/11/store/user/askew/EMTriggers/doublePho_run2012D_prompt_nvtx_pfmet_seltrig/pfPhoIDHists_65_1_meP.root"); // Andrew's nTuple


  // Disabling unused branches will speed up the processing significantly, but risks inconsistencies if wrong trees are turned off.
  // Make sure you know what you are doing.
  // Especially be careful when producing a skim - disabled branches will not be copied.
  // You can either enable all branches, or have two Event objects (one for event selection with reduced branch configuration, and the
  // other for event copying with branches fully enabled).
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("runNumber", 1);
  chain.SetBranchStatus("luminosityBlockNumber", 1);
  chain.SetBranchStatus("eventNumber", 1);
  chain.SetBranchStatus("isRealData", 1);
  chain.SetBranchStatus("metFilterBit", 1);
  chain.SetBranchStatus("rho", 1);
  chain.SetBranchStatus("rho25", 1);
  chain.SetBranchStatus("hlt*", 1);
  chain.SetBranchStatus("vertices*", 1);
  chain.SetBranchStatus("tracks*", 1);
  chain.SetBranchStatus("photons_photons*", 1);
  chain.SetBranchStatus("muons_muons*", 1);
  chain.SetBranchStatus("electrons_gsfElectrons*", 1);
  chain.SetBranchStatus("pfJets_ak5*", 1);
  //chain.SetBranchStatus("met_pfType01SysShiftCorrectedMet*", 1);
  //chain.SetBranchStatus("met_pfMet*", 1);// yutaro's idea on April 18, 2014
  chain.SetBranchStatus("met_pfType01CorrectedMet*", 1);// corrected from Yutaro's new code, commented out to check yutaro's idea on march 3, 2014, comment3d in on June 25, 2014 to get the fake rate
  //chain.SetBranchStatus("met_pfType1CorrectedMet*", 1);
  chain.SetBranchStatus("gridParams*", 1);

  if(chain.LoadTree(0) != 0){
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  SusyEventAnalyzer sea(chain);

  sea.SetOutput(outputName);
  sea.SetLogFile("cout"); // set to a full path to a file to output log to a file
  sea.SetPrintInterval(10000);
  sea.SetPrintLevel(0);
  //sea.AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50"); // uncomment if working with data
  //sea.AddHltName("HLT_Photon36_R9Id85_Photon22_R9Id85"); // adding new trigger to see if we are getting consistent result
  //sea.IncludeAJson(""); // put your favourite JSON
  sea.CopyEvents(false);
  sea.SetProcessNEvents(-1);

  TStopwatch ts;

  ts.Start();

  sea.Run();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
