// Original Author:  Dongwook Jang
// $Id: ana.C,v 1.12 2013/05/07 20:29:39 yiiyama Exp $
#include "TString.h"
#include "TChain.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>


void ana_2(TString outputName="analysis_SUSYSignal_2"){
   TChain chain("susyTree");//susyTree

  // REMOVE THE LINE BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
  gSystem->Load("libCondFormatsJetMETObjects.so");

  gSystem->Load("libSusyEvent.so");

  gSystem->AddIncludePath("-I" + TString(gSystem->Getenv("CMSSW_RELEASE_BASE")) + "/src");

  // Analysis macro
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");
  //gSystem->Load("SusyEventAnalyzer_cc.so");

  // chain of inputs
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v1/SMS-T5gg_2J_mGl-800to950_mLSP-25to925_TuneZ2star_8TeV-madgraph-tauola_Summer12-START53_V19_FSIM_PU_S12/susyEvents_T5gg_800_325_48.root");
  


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
