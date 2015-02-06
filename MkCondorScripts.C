#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
void MkCondorScripts(void){
  ifstream fnames("RootFilesB.txt");
  
  vector<string> fullpath;

  while (!fnames.eof()){
    string temp;
    if (!(fnames >> temp)) break;
    fullpath.push_back(temp);
  }
  cout << "read " << fullpath.size() << " files." << endl;

  ofstream BigJob("RunJobB.sh");
  for (int i=0;i<int(fullpath.size());++i){
    char ScriptName[240];
    sprintf(ScriptName, "condor_jobsB/Script_%i.sh",i);
    char JDLName[240];
    sprintf(JDLName, "condor_jobsB/AnaCondor_%i.jdl",i);
    ofstream scri(ScriptName);
    ofstream jdl(JDLName);
    BigJob << "condor_submit " << JDLName << endl;

    //Make the JDL file:
    jdl << "Universe                = vanilla" << endl;
    jdl << "Environment             = CONDORJOBID=$(Process)" << endl;
    jdl << "notification            = Error" << endl;
    jdl << "requirements            = Memory >= 199 &&OpSys == \"LINUX\"&& (Arch != \"DUMMY\" )" << endl;
    jdl << "Transfer_Input_Files    =  /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/ana.C,/uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/SusyEvent.h,/uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/JetCorrectorParameters.h, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/FactorizedJetCorrector.h, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/SusyEventAnalyzer.cc, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/SusyEventAnalyzer.h,/uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/SusyEventPrinter.h,/uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/SusyEventPrinter.cc, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/libSusyEvent.so, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/libJetMETObjects.so, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/reweights.root,/uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/fakeRateReweights.root, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/fakeRateReweights_2D.root,  /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/PUweightsForBino.root, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/MetAndMassShapes.root,/uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/Ratio.root, /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/Cert_190456-208686_8TeV_PromptReco_Collisions12.txt" << endl;
    jdl << "should_transfer_files   = YES" << endl;
    jdl << "when_to_transfer_output = ON_EXIT" << endl;

    jdl << "Executable              = /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/" << ScriptName << endl;
    jdl << "output                  = /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/logsB/an_" << i << ".out" << endl;
    jdl << "error                   = /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/logsB/an_" << i <<".err" << endl;
    jdl << "log                     = /uscms_data/d3/asantra4/CMSSW_5_3_3/src/SusyAnalysis/SusyNtuplizer/macroNEW/logsB/an_"<<i<<".log" << endl;
    jdl << "Queue 1" << endl;
    //Done with condor JDL

    //Write analysis script
    scri << "#! /bin/bash" << endl;
    scri << "source /uscmst1/prod/sw/cms/shrc prod" << endl;
    scri << "export SCRAM_ARCH=slc5_amd64_gcc462" << endl;
    scri << "cd /uscms_data/d3/asantra4/CMSSW_5_3_3/src" << endl;
    scri << "eval `scramv1 runtime -sh`" << endl;
    scri << "cd -" << endl;
    scri << "root -l -b << EOF " << endl;
    scri << ".L ana.C "<< endl;
    scri << "ana(\"" <<  fullpath[i] << "\", \"" << i << "\" );" << endl;//
    scri << ".q" << endl;
    scri << "EOF" << endl;
    /*scri << "   TString makeshared(gSystem->GetMakeSharedLib());" << endl;
    scri << "   TString dummy = makeshared.ReplaceAll(\"-W \", \"\");" << endl;
    scri << "   gSystem->SetMakeSharedLib(makeshared);" << endl;
    scri << "   cout<<\"here \"<<endl;" << endl;
    scri << "   .L UCandidateSkimmer.cpp++" << endl;
    scri << "   UCandidateSkimmer(\""<<fullpath[i] <<"\");" << endl;
    scri << "   .q" << endl;*/
    //  scri << "EOF" <<endl;

  }
    




}
