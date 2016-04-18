#ifndef EvtWgtAnalysis_h
#define EvtWgtAnalysis_h

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TH1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"

#include "AnaTree/AnaBNB.h"

class EvtWgtAnalysis {
  public :
  
  AnaBNB * anaBNBtree;
  
  int Nfiles = 0;
  int Nfunctions = 0;
  vector<string> *functionsName;
  vector<TH1D*> histoPmuVec;            // histograms for Pmu (1 for each function used)
  vector<TH1D*> histoPmuVec_p1;         // histograms for Pmu (+ 1sigma)
  vector<TH1D*> histoPmuVec_m1;         // histograms for Pmu (- 1sigma)
  vector<TH1D*> histoCosThetaMuVec;     // histograms for CosThetaMu (1 for each function used)
  vector<TH1D*> histoCosThetaMuVec_p1;  // histograms for CosThetaMu (+ 1sigma)
  vector<TH1D*> histoCosThetaMuVec_m1;  // histograms for CosThetaMu (- 1sigma)
  int treeNumber = -1;
  
  
  // Efficiency
  
  double all_evts_nominal = 0;  // all numu cc events (w/o selection) -> nominal values of the parameters
  vector<double> all_evts_p1;   // all numu cc events (w/o selection) -> + 1 sigma parameters
  vector<double> all_evts_m1;   // all numu cc events (w/o selection) -> - 1 sigma parameters
  
  double sel_evts_nominal = 0;  // selected numu cc events (w/o selection) -> nominal values of the parameters
  vector<double> sel_evts_p1;   // selected numu cc events (w/o selection) -> + 1 sigma parameters
  vector<double> sel_evts_m1;   // selected numu cc events (w/o selection) -> - 1 sigma parameters
  
  TH1F *xsec_mom_truth = new TH1F("xsec_mom_truth", "", 10, 0, 2);
  TH1F *xsec_mom_data = new TH1F("xsec_mom_data", "", 10, 0, 2);
  TH1F *xsec_mom_bg = new TH1F("xsec_mom_bg", "", 10, 0, 2);
  TH1F *xsec_mom_eff = new TH1F("xsec_mom_eff", "", 10, 0, 2);
  
  TH1F *xsec_mom_reco_truth = new TH1F("xsec_mom_reco_truth", "", 10, 0, 2);
  TH1F *xsec_mom_reco_data = new TH1F("xsec_mom_reco_data", "", 10, 0, 2);
  TH1F *xsec_mom_reco_bg = new TH1F("xsec_mom_reco_bg", "", 10, 0, 2);
  TH1F *xsec_mom_reco_eff = new TH1F("xsec_mom_reco_eff", "", 10, 0, 2);
  
  TH1F *xsec_theta_truth = new TH1F("xsec_theta_truth", "", 10, -1, 1);
  TH1F *xsec_theta_data = new TH1F("xsec_theta_data", "", 10, -1, 1);
  TH1F *xsec_theta_bg = new TH1F("xsec_theta_bg", "", 10, -1, 1);
  TH1F *xsec_theta_eff = new TH1F("xsec_theta_eff", "", 10, -1, 1);
  
  TH1F *xsec_theta_reco_truth = new TH1F("xsec_theta_reco_truth", "", 10, -1, 1);
  TH1F *xsec_theta_reco_data = new TH1F("xsec_theta_reco_data", "", 10, -1, 1);
  TH1F *xsec_theta_reco_bg = new TH1F("xsec_theta_reco_bg", "", 10, -1, 1);
  TH1F *xsec_theta_reco_eff = new TH1F("xsec_theta_reco_eff", "", 10, -1, 1);
  
  vector<TH1D*> xsec_mom_truth_p1;
  vector<TH1D*> xsec_mom_truth_m1;
  vector<TH1D*> xsec_theta_truth_p1;
  vector<TH1D*> xsec_theta_truth_m1;
  vector<TH1D*> xsec_mom_data_p1;
  vector<TH1D*> xsec_mom_data_m1;
  vector<TH1D*> xsec_mom_reco_data_p1;
  vector<TH1D*> xsec_mom_reco_data_m1;
  vector<TH1D*> xsec_theta_data_p1;
  vector<TH1D*> xsec_theta_data_m1;
  vector<TH1D*> xsec_theta_reco_data_p1;
  vector<TH1D*> xsec_theta_reco_data_m1;
  
  vector<TH1D*> xsec_mom_eff_p1;
  vector<TH1D*> xsec_mom_eff_m1;
  vector<TH1D*> xsec_mom_reco_eff_p1;
  vector<TH1D*> xsec_mom_reco_eff_m1;
  vector<TH1D*> xsec_theta_eff_p1;
  vector<TH1D*> xsec_theta_eff_m1;
  vector<TH1D*> xsec_theta_reco_eff_p1;
  vector<TH1D*> xsec_theta_reco_eff_m1;
  
  
  TH1F *pmu_numu_cc_reco_histo = new TH1F("pmu_numu_cc_reco_histo", "", 10, 0, 2);
  TH1F *costhetamu_numu_cc_reco_histo = new TH1F("costhetamu_numu_cc_reco_histo", "", 10, -1, 1);
  TH1F *pmu_anumu_cc_reco_histo = new TH1F("pmu_anumu_cc_reco_histo", "", 10, 0, 2);
  TH1F *costhetamu_anumu_cc_reco_histo = new TH1F("costhetamu_anumu_cc_reco_histo", "", 10, -1, 1);
  TH1F *pmu_nue_cc_reco_histo = new TH1F("pmu_nue_cc_reco_histo", "", 10, 0, 2);
  TH1F *costhetamu_nue_cc_reco_histo = new TH1F("costhetamu_nue_cc_reco_histo", "", 10, -1, 1);
  TH1F *pmu_nc_reco_histo = new TH1F("pmu_nc_reco_histo", "", 10, 0, 2);
  TH1F *costhetamu_nc_reco_histo = new TH1F("costhetamu_nc_reco_histo", "", 10, -1, 1);

  
  
  
  vector<TH1D*> pmu_numu_cc_reco_histo_p1;
  vector<TH1D*> pmu_numu_cc_reco_histo_m1;
  vector<TH1D*> costhetamu_numu_cc_reco_histo_p1;
  vector<TH1D*> costhetamu_numu_cc_reco_histo_m1;
  vector<TH1D*> pmu_anumu_cc_reco_histo_p1;
  vector<TH1D*> pmu_anumu_cc_reco_histo_m1;
  vector<TH1D*> costhetamu_anumu_cc_reco_histo_p1;
  vector<TH1D*> costhetamu_anumu_cc_reco_histo_m1;
  vector<TH1D*> pmu_nue_cc_reco_histo_p1;
  vector<TH1D*> pmu_nue_cc_reco_histo_m1;
  vector<TH1D*> costhetamu_nue_cc_reco_histo_p1;
  vector<TH1D*> costhetamu_nue_cc_reco_histo_m1;
  vector<TH1D*> pmu_nc_reco_histo_p1;
  vector<TH1D*> pmu_nc_reco_histo_m1;
  vector<TH1D*> costhetamu_nc_reco_histo_p1;
  vector<TH1D*> costhetamu_nc_reco_histo_m1;
  
  
  // Cross Section
  TH1D *        XSec_pmu_nominal;
  vector<TH1D*> XSec_pmu_p1;
  vector<TH1D*> XSec_pmu_m1;
  
  TH1D * background_pmu_nominal;
  vector<TH1D*> background_pmu_p1;
  vector<TH1D*> background_pmu_m1;
  
  vector<TH1D*> XSec_pmu_percDiff_p1;
  vector<TH1D*> XSec_pmu_percDiff_m1;
  // ***

  
  
  
  TH1D *temp = new TH1D("temp", "", 100, 0,100);
  int Ncoh = 0; // number of Coherent production interactions
  
  TChain *cflux;
  
  EvtWgtAnalysis(string pattern="/pnfs/uboone/scratch/users/mdeltutt/v04_30_03/anatree_bnb_eventWeight_all/reunion/standard_reco_hist*.root"/*"/uboone/app/users/mdeltutt/eventWeight/standard_reco_hist.root"*/);
  virtual ~EvtWgtAnalysis();
  
  void CalcEfficiency();
  bool inFV(double x, double y, double z);
  void MakeHistograms();
  //void MakePlotsPmu(bool normalised = false);
  void MakePlots(bool normalised = false, int variable = 0);
  void MakeBackgroundPlots(int variable = 0);
  void CalculateXSecPercDifference();
  void MakeXsecDiffPlots();
  //void MakePlotsCosThetaMu(bool normalised = false);
  TString GetLegendName(string fName);
  void InstantiateHistograms(Int_t nFunc, vector<string> *funcName);
  void InstantiateHistogramsEfficiency(Int_t nFunc, vector<string> *funcName);
};

#endif

#ifdef EvtWgtAnalysis_cxx

EvtWgtAnalysis::EvtWgtAnalysis(string pattern) {
  
  const char* path = "/uboone/app/users/mdeltutt/eventWeight/AnaTree/";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    //libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());
    gSystem->Load("AnaBNB_C.so");
  }
  
  cflux = new TChain("analysistree/anatree");
  cflux->Add(pattern.c_str());
  
  Nfiles = cflux->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;
  
}

EvtWgtAnalysis::~EvtWgtAnalysis() {
  
}

#endif // #ifdef EvtWgtAnalysis_cxx
