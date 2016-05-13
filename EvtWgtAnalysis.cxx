#define EvtWgtAnalysis_cxx
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

using namespace std;

#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TRotation.h"
#include "TMath.h"
#include "THStack.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"

#include "EvtWgtAnalysis.h"
#include "AnaTree/AnaBNB.h"


void EvtWgtAnalysis::CalcEfficiency() {
  
  std::cout << "EvtWgtAnalysis::CalcEfficiency starts" << std::endl;
  
  
  const int maxentries = 35000;
  const int maxtracks = 2000;
  const int maxvtx = 200;
  const int maxnu = 10;
  const int kMaxFlashes = 2000;
  
  
  anaBNBtree = new AnaBNB(cflux);
  
  // It's important to active only the branches we are interested in
  // Some branches cause segmentation violation...
  anaBNBtree->fChain->SetBranchStatus("*",0);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_funcname",1);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_nfunc",1);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_weight",1);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_nweight",1);
  anaBNBtree->fChain->SetBranchStatus("genie_primaries_pdg",1);
  anaBNBtree->fChain->SetBranchStatus("genie_P",1);
  anaBNBtree->fChain->SetBranchStatus("genie_Pz",1);
  anaBNBtree->fChain->SetBranchStatus("genie_status_code",1);
  anaBNBtree->fChain->SetBranchStatus("genie_no_primaries",1);
  anaBNBtree->fChain->SetBranchStatus("ccnc_truth",1);
  anaBNBtree->fChain->SetBranchStatus("nuPDG_truth",1);
  anaBNBtree->fChain->SetBranchStatus("mode_truth",1);
  anaBNBtree->fChain->SetBranchStatus("nuvtxx_truth",1);
  anaBNBtree->fChain->SetBranchStatus("nuvtxy_truth",1);
  anaBNBtree->fChain->SetBranchStatus("nuvtxz_truth",1);
  anaBNBtree->fChain->SetBranchStatus("lep_mom_truth",1);
  anaBNBtree->fChain->SetBranchStatus("lep_dcosz_truth",1);
  anaBNBtree->fChain->SetBranchStatus("enu_truth",1);
  anaBNBtree->fChain->SetBranchStatus("ntracks_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkstartx_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkstarty_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkstartz_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkendx_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkendy_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkendz_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trklen_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trktheta_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkmomrange_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkcosmicscore_tagger_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkcosmicscore_flashmatch_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkorigin_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("trkpidbestplane_pandoraNuKHit",1);
  anaBNBtree->fChain->SetBranchStatus("nvtx_pandoraNu",1);
  anaBNBtree->fChain->SetBranchStatus("vtxx_pandoraNu",1);
  anaBNBtree->fChain->SetBranchStatus("vtxy_pandoraNu",1);
  anaBNBtree->fChain->SetBranchStatus("vtxz_pandoraNu",1);
  anaBNBtree->fChain->SetBranchStatus("no_flashes",1);
  anaBNBtree->fChain->SetBranchStatus("flash_time",1);
  anaBNBtree->fChain->SetBranchStatus("flash_pe",1);
  
  
  int Size = cflux -> GetEntries();
  cout << "number of events used is: " << Size << endl;
  
  bool contained = false;
  int index_longest_track = -1;
  double longest_track = 0;
  
  int vbest = 0;
  double distmin = 100000;
  double dist = 0;
  bool            vertexatstart[maxtracks];
  bool            vertexatend[maxtracks];
  
  int n_neutrinotag = 0;
  int n_nvtxtag = 0;
  int n_trackstag = 0;
  int n_eventstag = 0;
  int n_longesttag = 0;
  int n_longestnutag = 0;
  int n_longestcosmictag = 0;
  int n_longestothertag = 0;
  int n_longestSNtag = 0;
  int n_nc = 0;
  int n_numubar = 0;
  int n_nue = 0;
  bool neutrinotag = false;
  bool eventtag = false;
  bool tracktag = false;
  int n_trackclose = 0;
  int n_flashmatchtag = 0;
  int n_flashtag = 0;
  int n_geometrytag = 0;
  bool geometrytag = false;
  bool flashtag = false;
  bool flashmatchtag = false;
  
  
  // Looping over the entries now
  for(int i = 0; i < Size; i++) {
    if(i%1000 == 0) cout << "\t... " << i << endl;
    cflux->GetEntry(i);
    
    if(i == 0) {
      InstantiateHistogramsEfficiency(anaBNBtree->evtwgt_nfunc, anaBNBtree->evtwgt_funcname);
    }
    
    // Coping here the relevant variables
    int k = 0;
    Int_t           event;
    Int_t           ccnc_truth = anaBNBtree->ccnc_truth[k];
    Int_t           mode_truth = anaBNBtree->mode_truth[k];
    Int_t           nuPDG_truth = anaBNBtree->nuPDG_truth[k];
    Float_t         nuvtxx_truth = anaBNBtree->nuvtxx_truth[k];
    Float_t         nuvtxy_truth = anaBNBtree->nuvtxy_truth[k];
    Float_t         nuvtxz_truth = anaBNBtree->nuvtxz_truth[k];
    Float_t         lep_mom_truth = anaBNBtree->lep_mom_truth[k];
    Float_t         lep_dcosz_truth = anaBNBtree->lep_dcosz_truth[k];
    Float_t         enu_truth = anaBNBtree->enu_truth[k];
    
    Int_t           ntracks_reco = anaBNBtree->ntracks_pandoraNuKHit;
    Short_t         trkorigin[maxtracks][3];
    for (int itrk=0; itrk<ntracks_reco; itrk++) {
      trkorigin[itrk][0] = anaBNBtree->trkorigin_pandoraNuKHit[itrk][0];
      trkorigin[itrk][1] = anaBNBtree->trkorigin_pandoraNuKHit[itrk][1];
      trkorigin[itrk][2] = anaBNBtree->trkorigin_pandoraNuKHit[itrk][2];
    }
    
    //    Short_t         trkbestplane[ntracks_reco]; trkbestplane= anaBNBtree->trkpidbestplane_pandoraNuKHit;
    //  Float_t         *trklen = anaBNBtree->trklen_pandoraNuKHit;   //[ntracks_trackkalmanhit]
    //  cout<<"i = "<<i<<"  trklen: "<<trklen[0]<<" "<<trklen[1]<<" "<<trklen[2]<<endl;
    
    /*    Float_t         trkstartx[ntracks_reco];  trkstartx = anaBNBtree->trkstartx_pandoraNuKHit;
     Float_t         *trkstarty = anaBNBtree->trkstarty_pandoraNuKHit;
     Float_t         *trkstartz = anaBNBtree->trkstartz_pandoraNuKHit;
     Float_t         *trkendx = anaBNBtree->trkendx_pandoraNuKHit;
     Float_t         *trkendy = anaBNBtree->trkendy_pandoraNuKHit;
     Float_t         *trkendz = anaBNBtree->trkendz_pandoraNuKHit;
     Float_t         *trktheta = anaBNBtree->trktheta_pandoraNuKHit;
     Float_t         *trkmomrange = anaBNBtree->trkmomrange_pandoraNuKHit;
     Float_t         *trkcosmicscore_tagger = anaBNBtree->trkcosmicscore_tagger_pandoraNuKHit;
     Float_t         *trkcosmicscore_flashmatch = anaBNBtree->trkcosmicscore_flashmatch_pandoraNuKHit;
     */
    Short_t         nvtx = anaBNBtree->nvtx_pandoraNu;
    Float_t         vtxx[maxvtx];
    Float_t         vtxy[maxvtx];
    Float_t         vtxz[maxvtx];
    
    for (int ivertex=0; ivertex<nvtx; ivertex++) {
      vtxx[ivertex] = anaBNBtree->vtxx_pandoraNu[ivertex];
      vtxy[ivertex] = anaBNBtree->vtxy_pandoraNu[ivertex];
      vtxz[ivertex] = anaBNBtree->vtxz_pandoraNu[ivertex];
    }
    
    //Int_t           no_flashes = anaBNBtree->no_flashes;                            //number of flashes in the event
    //Float_t         *flash_time = anaBNBtree->flash_time;       //flash time (in microseconds)
    //Float_t         *flash_pe = anaBNBtree->flash_pe;          //total number of photoelectrons corresponding to the flash
    
    
    
    if(ccnc_truth == 0 && nuPDG_truth == 14 && enu_truth > 0.4) {
      if(inFV(nuvtxx_truth, nuvtxy_truth, nuvtxz_truth)) {
        xsec_mom_truth -> Fill(lep_mom_truth);
        xsec_theta_truth -> Fill(lep_dcosz_truth);
        all_evts_nominal++;
        for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
          all_evts_p1[function] += (anaBNBtree->evtwgt_weight->at(function)).at(0);
          all_evts_m1[function] += (anaBNBtree->evtwgt_weight->at(function)).at(1);
          
          xsec_mom_truth_p1[function] -> Fill(lep_mom_truth, (anaBNBtree->evtwgt_weight->at(function)).at(0));
          xsec_mom_truth_m1[function] -> Fill(lep_mom_truth, (anaBNBtree->evtwgt_weight->at(function)).at(1));
          xsec_theta_truth_p1[function] -> Fill(lep_dcosz_truth, (anaBNBtree->evtwgt_weight->at(function)).at(0));
          xsec_theta_truth_m1[function] -> Fill(lep_dcosz_truth, (anaBNBtree->evtwgt_weight->at(function)).at(1));
        } // end loop re-weighting functions
      }
    }
    
    distmin = 100000;
    contained = true;
    longest_track = 0;
    index_longest_track = -1;
    neutrinotag = false;
    eventtag = false;
    geometrytag = false;
    flashtag = false;
    flashmatchtag = false;
    tracktag = false;
    
    
    for(int j = 0; j < anaBNBtree->ntracks_pandoraNuKHit; j++) {
      if(anaBNBtree->trkcosmicscore_tagger_pandoraNuKHit[j] < 0.4) geometrytag = true;
      if(anaBNBtree->trkcosmicscore_flashmatch_pandoraNuKHit[j] < 0.5) flashmatchtag = true;
    }
    
    for(int f = 0; f < anaBNBtree->no_flashes; f++) {
      if((anaBNBtree->flash_time[f] > 0 && anaBNBtree->flash_time[f] < 1.6) && anaBNBtree->flash_pe[f] > 50) flashtag = true;
      //  if(i==1) cout << f << " " << flash_time[f] << " " << flash_pe[f] << endl;
    }
    
    
    if(flashtag == true && geometrytag == true && flashmatchtag == true) neutrinotag = true;
    
    if(neutrinotag == true) {
      n_neutrinotag++;
      if(nvtx > 0) n_nvtxtag++;
    }
    
    if(neutrinotag == true) {
      for(int j = 0; j < anaBNBtree->ntracks_pandoraNuKHit; j++) {
        if((anaBNBtree->trkcosmicscore_flashmatch_pandoraNuKHit[j] < 0.5) && (anaBNBtree->trkcosmicscore_tagger_pandoraNuKHit[j] < 0.4)) {
          if(nvtx > 0) {
            distmin = 100000;
            //loop over all vertices and calculate the closest distance of the longest contained track to each vertex. find the minimum
            for(int v = 0; v < nvtx; v++) {
              dist = sqrt((vtxx[v] - anaBNBtree->trkstartx_pandoraNuKHit[j])*(vtxx[v] - anaBNBtree->trkstartx_pandoraNuKHit[j]) + (vtxy[v] - anaBNBtree->trkstarty_pandoraNuKHit[j])*(vtxy[v] - anaBNBtree->trkstarty_pandoraNuKHit[j]) + (vtxz[v] - anaBNBtree->trkstartz_pandoraNuKHit[j])*(vtxz[v] - anaBNBtree->trkstartz_pandoraNuKHit[j]));
              if(dist < distmin) {
                distmin = dist;
                vbest = v;
                vertexatstart[j] = true;
                vertexatend[j] = false;
              }
              dist = sqrt((vtxx[v] - anaBNBtree->trkendx_pandoraNuKHit[j])*(vtxx[v] - anaBNBtree->trkendx_pandoraNuKHit[j]) + (vtxy[v] - anaBNBtree->trkendy_pandoraNuKHit[j])*(vtxy[v] - anaBNBtree->trkendy_pandoraNuKHit[j]) + (vtxz[v] - anaBNBtree->trkendz_pandoraNuKHit[j])*(vtxz[v] - anaBNBtree->trkendz_pandoraNuKHit[j]));
              if(dist < distmin) {
                distmin = dist;
                vbest = v;
                vertexatstart[j] = false;
                vertexatend[j] = true;
              }
            }//loop over all vertices
            if(distmin < 5 && inFV(vtxx[vbest], vtxy[vbest], vtxz[vbest])) { //cut on 5cm distance longest track (start or end) to closest reco vertex
              if(eventtag == false) eventtag = true;
              if(anaBNBtree->trklen_pandoraNuKHit[j] > longest_track) {
                longest_track = anaBNBtree->trklen_pandoraNuKHit[j];
                index_longest_track = j;
              }
            }
          }//nvtx > 0
        }//neutrino tag
      }//loop over all tracks
      //cout<<"longest_track: "<<longest_track<<endl;
      if(longest_track > 0) {
        if(longest_track > 75) {
          if(trkorigin[index_longest_track][anaBNBtree->trkpidbestplane_pandoraNuKHit[index_longest_track]] == 1) {
            // If I got till this point, it means I have selected the event
            // Fill the histograms at this point
            
            // Save relevant variables here for convenience
            double momentum_reco = anaBNBtree->trkmomrange_pandoraNuKHit[index_longest_track];
            double angle_reco;
            if(vertexatstart[index_longest_track]) angle_reco = cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]);
            if(vertexatend[index_longest_track])   angle_reco = -cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]);
            
            if(ccnc_truth == 0 && nuPDG_truth == 14) {
              
              // For the efficiency
              sel_evts_nominal ++;
              
              xsec_mom_data -> Fill(lep_mom_truth);
              xsec_mom_reco_data -> Fill(momentum_reco);
              xsec_theta_data -> Fill(lep_dcosz_truth);
              xsec_theta_reco_data -> Fill(angle_reco);
              
              // For the number of events
              pmu_numu_cc_reco_histo -> Fill(momentum_reco);
              costhetamu_numu_cc_reco_histo -> Fill(angle_reco);
              
              // Loop over the re-weighting functions and fill +-1 sigma histograms
              for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
                // For the efficiency
                sel_evts_p1[function] += (anaBNBtree->evtwgt_weight->at(function)).at(0);
                sel_evts_m1[function] += (anaBNBtree->evtwgt_weight->at(function)).at(1);
                
                xsec_mom_data_p1[function] -> Fill(lep_mom_truth, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_mom_data_m1[function] -> Fill(lep_mom_truth, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                xsec_mom_reco_data_p1[function] -> Fill(momentum_reco,(anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_mom_reco_data_m1[function] -> Fill(momentum_reco,(anaBNBtree->evtwgt_weight->at(function)).at(1));
                xsec_theta_data_p1[function] -> Fill(lep_dcosz_truth, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_theta_data_m1[function] -> Fill(lep_dcosz_truth, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                xsec_theta_reco_data_p1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_theta_reco_data_m1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                
                // For the number of events
                pmu_numu_cc_reco_histo_p1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                pmu_numu_cc_reco_histo_m1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                costhetamu_numu_cc_reco_histo_p1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                costhetamu_numu_cc_reco_histo_m1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                
                
              } // end loop over re-weighting functions
            } // end numu cc
            
            else if(ccnc_truth == 0 && nuPDG_truth == -14) {
              pmu_anumu_cc_reco_histo -> Fill(momentum_reco);
              costhetamu_anumu_cc_reco_histo -> Fill(angle_reco);
              for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
                pmu_anumu_cc_reco_histo_p1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                pmu_anumu_cc_reco_histo_m1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                costhetamu_anumu_cc_reco_histo_p1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                costhetamu_anumu_cc_reco_histo_m1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
              } // end loop over re-weighting functions
            } // end anumu cc
            
            else if(ccnc_truth == 0 && (nuPDG_truth == 12 || nuPDG_truth == -12)) {
              pmu_nue_cc_reco_histo -> Fill(momentum_reco);
              costhetamu_nue_cc_reco_histo -> Fill(angle_reco);
              for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
                pmu_nue_cc_reco_histo_p1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                pmu_nue_cc_reco_histo_m1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                costhetamu_nue_cc_reco_histo_p1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                costhetamu_nue_cc_reco_histo_m1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
              } // end loop over re-weighting functions
            } // end nue anue cc
            
            else if (ccnc_truth == 1) {
              pmu_nc_reco_histo -> Fill(momentum_reco);
              costhetamu_nc_reco_histo -> Fill(angle_reco);
              for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
                pmu_nc_reco_histo_p1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                pmu_nc_reco_histo_m1[function] -> Fill(momentum_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                costhetamu_nc_reco_histo_p1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                costhetamu_nc_reco_histo_m1[function] -> Fill(angle_reco, (anaBNBtree->evtwgt_weight->at(function)).at(1));
              } // end loop over re-weighting functions
            } // end nc
            
          }
        } // end longest_track > 75
      } // end longest_track > 0
    } // end neutrinotag == true
  } // end loop over entries
  
  
  
  //*********************************
  //
  // Calculate efficiency
  //
  //*********************************
  
  // NOMINAL
  xsec_mom_eff -> Add(xsec_mom_data);
  xsec_mom_eff -> Divide(xsec_mom_truth);
  xsec_theta_eff -> Add(xsec_theta_data);
  xsec_theta_eff -> Divide(xsec_theta_truth);
  xsec_mom_reco_eff -> Add(xsec_mom_reco_data);
  xsec_mom_reco_eff -> Divide(xsec_mom_truth);
  xsec_theta_reco_eff -> Add(xsec_theta_reco_data);
  xsec_theta_reco_eff -> Divide(xsec_theta_truth);
  
  // +- 1 SIGMA
  for (unsigned int function = 0; function < xsec_mom_eff_p1.size(); function++) {
    // e as a function of true muon momentum
    xsec_mom_eff_p1[function] -> Add(xsec_mom_data_p1[function]);
    xsec_mom_eff_p1[function] -> Divide(xsec_mom_truth_p1[function]);
    xsec_mom_eff_m1[function] -> Add(xsec_mom_data_m1[function]);
    xsec_mom_eff_m1[function] -> Divide(xsec_mom_truth_m1[function]);
    
    // e as a function of true muon angle
    xsec_theta_eff_p1[function] -> Add(xsec_theta_data_p1[function]);
    xsec_theta_eff_p1[function] -> Divide(xsec_theta_truth_p1[function]);
    xsec_theta_eff_m1[function] -> Add(xsec_theta_data_m1[function]);
    xsec_theta_eff_m1[function] -> Divide(xsec_theta_truth_m1[function]);
    
    // e as a function of reco muon momentum
    xsec_mom_reco_eff_p1[function] -> Add(xsec_mom_reco_data_p1[function]);
    xsec_mom_reco_eff_p1[function] -> Divide(xsec_mom_truth_p1[function]);
    xsec_mom_reco_eff_m1[function] -> Add(xsec_mom_reco_data_m1[function]);
    xsec_mom_reco_eff_m1[function] -> Divide(xsec_mom_truth_m1[function]);
    
    // e as a function of reco muon angle
    xsec_theta_reco_eff_p1[function] -> Add(xsec_theta_reco_data_p1[function]);
    xsec_theta_reco_eff_p1[function] -> Divide(xsec_theta_truth_p1[function]);
    xsec_theta_reco_eff_m1[function] -> Add(xsec_theta_reco_data_m1[function]);
    xsec_theta_reco_eff_m1[function] -> Divide(xsec_theta_truth_m1[function]);
  }
  
  
  
  
  //*********************************
  //
  // Save efficiency histograms to file
  //
  //*********************************
  
  TFile f = TFile("eventWeightEfficiency.root", "RECREATE");
  
  // Save nominal
  xsec_mom_eff->Write();
  xsec_theta_eff->Write();
  xsec_mom_reco_eff->Write();
  xsec_theta_reco_eff->Write();
  
  xsec_mom_data->Write();
  xsec_mom_truth->Write();
  
  
  // Save +-1 sigma
  for (unsigned int function = 0; function < xsec_mom_eff_p1.size(); function++) {
    xsec_mom_truth_p1[function]->Write();
    xsec_mom_data_p1[function]->Write();
    
    xsec_mom_eff_p1[function]->Write();
    xsec_mom_eff_m1[function]->Write();
    xsec_theta_eff_p1[function]->Write();
    xsec_theta_eff_m1[function]->Write();
    xsec_mom_reco_eff_p1[function]->Write();
    xsec_mom_reco_eff_m1[function]->Write();
    xsec_theta_reco_eff_p1[function]->Write();
    xsec_theta_reco_eff_m1[function]->Write();
  }
  
  f.Close();
  std::cout << "Efficiency histograms saved to eventWeightEfficiency.root" << std::endl;
  
  
  
  //*********************************
  //
  // Save events and backgrounds to file
  //
  //*********************************
  
  TFile ff = TFile("eventWeightBackground.root", "RECREATE");
  
  // Save nominal
  pmu_numu_cc_reco_histo -> Write();
  costhetamu_numu_cc_reco_histo -> Write();
  pmu_anumu_cc_reco_histo -> Write();
  costhetamu_anumu_cc_reco_histo -> Write();
  pmu_nue_cc_reco_histo -> Write();
  costhetamu_nue_cc_reco_histo -> Write();
  pmu_nc_reco_histo -> Write();
  costhetamu_nc_reco_histo -> Write();
  
  // Save +-1 sigma
  for (unsigned int function = 0; function < pmu_numu_cc_reco_histo_p1.size(); function++) {
    pmu_numu_cc_reco_histo_p1[function] -> Write();
    pmu_numu_cc_reco_histo_m1[function] -> Write();
    costhetamu_numu_cc_reco_histo_p1[function] -> Write();
    costhetamu_numu_cc_reco_histo_m1[function] -> Write();
    
    pmu_anumu_cc_reco_histo_p1[function] -> Write();
    pmu_anumu_cc_reco_histo_m1[function] -> Write();
    costhetamu_anumu_cc_reco_histo_p1[function] -> Write();
    costhetamu_anumu_cc_reco_histo_m1[function] -> Write();
    
    pmu_nue_cc_reco_histo_p1[function] -> Write();
    pmu_nue_cc_reco_histo_m1[function] -> Write();
    costhetamu_nue_cc_reco_histo_p1[function] -> Write();
    costhetamu_nue_cc_reco_histo_m1[function] -> Write();
    
    pmu_nc_reco_histo_p1[function] -> Write();
    pmu_nc_reco_histo_m1[function] -> Write();
    costhetamu_nc_reco_histo_p1[function] -> Write();
    costhetamu_nc_reco_histo_m1[function] -> Write();
    
  }
  
  ff.Close();
  std::cout << "Background histograms saved to eventWeightBackground.root" << std::endl;
  
  
  
  
}


//_______________________________________________________
bool EvtWgtAnalysis::inFV(double x, double y, double z) {
  
  double FVx = 256.35;
  double FVy = 233;
  double FVz = 1036.8;
  double border = 5.;
  
  if(x < (FVx - border) && x > border && y < (FVy/2. - border) && y > (-FVy/2. + border) && z < (FVz - border) && z > border) return true;
  else return false;
}





//__________________________________________________
void EvtWgtAnalysis::CalculateXSecPercDifference() {
  
  std::cout << "EvtWgtAnalysis::CalculateXSecPercDifference starts" << std::endl;
  
  std::cout << "Calling CalcEfficiency() first." << std::endl;
  CalcEfficiency();
  
  // Instantiate XSec and background histograms
  
  XSec_pmu_nominal = new TH1D("XSec_pmu_nominal", ";p_{#mu};#propto#sigma [arb.]", 10, 0, 2);
  background_pmu_nominal = new TH1D("background_pmu_nominal", ";p_{#mu}; Background events", 10, 0, 2);
  
  int nFunc = pmu_numu_cc_reco_histo_p1.size();
  
  TString nameBase;
  
  XSec_pmu_p1.resize(functionsName->size());
  XSec_pmu_m1.resize(functionsName->size());
  background_pmu_p1.resize(functionsName->size());
  background_pmu_m1.resize(functionsName->size());
  
  for ( int j=0; j<nFunc; j++) {
    nameBase = "XSec_pmu_";
    XSec_pmu_p1.at(j) = new TH1D(nameBase+functionsName->at(j)+"_Plus1Sigma",  ";p_{#mu};#propto#sigma^+ [arb.]", 10, 0, 2);
    XSec_pmu_m1.at(j) = new TH1D(nameBase+functionsName->at(j)+"_Minus1Sigma", ";p_{#mu};#propto#sigma^- [arb.]", 10, 0, 2);
    
    nameBase = "B_pmu_";
    background_pmu_p1.at(j) = new TH1D(nameBase+functionsName->at(j)+"_Plus1Sigma",  ";p_{#mu};Background events", 10, 0, 2);
    background_pmu_m1.at(j) = new TH1D(nameBase+functionsName->at(j)+"_Minus1Sigma", ";p_{#mu};Background events", 10, 0, 2);
  }
  
  // Nominal background
  
  background_pmu_nominal -> Add(pmu_anumu_cc_reco_histo);  // anumu
  background_pmu_nominal -> Add(pmu_nue_cc_reco_histo);  // anumu
  background_pmu_nominal -> Add(pmu_nc_reco_histo);  // anumu

  // Summing up the baground events for +-1 sigma parameters
  // Loop over the functions
  for ( int j=0; j<nFunc; j++) {
    
    background_pmu_p1.at(j) -> Add(pmu_anumu_cc_reco_histo_p1.at(j));  // anumu
    background_pmu_p1.at(j) -> Add(pmu_nue_cc_reco_histo_p1.at(j));    // nue
    background_pmu_p1.at(j) -> Add(pmu_nc_reco_histo_p1.at(j));        // nc
    
    background_pmu_m1.at(j) -> Add(pmu_anumu_cc_reco_histo_m1.at(j));  // anumu
    background_pmu_m1.at(j) -> Add(pmu_nue_cc_reco_histo_m1.at(j));    // nue
    background_pmu_m1.at(j) -> Add(pmu_nc_reco_histo_m1.at(j));        // nc

  } // endl loop over functions
  
  // I now have background, efficiency and nominal events
  // Background: background_pmu_nominal
  //             background_pmu_p1[]
  //             background_pmu_m1[]
  //
  // Efficiency: xsec_mom_eff
  //             xsec_mom_eff_p1[]
  //             xsec_mom_eff_m1[]
  //
  // Events:     pmu_numu_cc_reco_histo   (only in this case, we are assuming we take this value from data)
  
  // XSec: nominal
  XSec_pmu_nominal->Add(pmu_numu_cc_reco_histo);
  XSec_pmu_nominal->Add(background_pmu_nominal, -1);
  XSec_pmu_nominal->Divide(xsec_mom_eff);

  // XSec: +- 1 sigma
  std::cout << "Calculating cross section with efficiency as a function of *****TRUE***** muon momentum." << std::endl;
  for (int j=0; j<nFunc; j++) {
    XSec_pmu_p1.at(j)->Add(pmu_numu_cc_reco_histo); // Nominal (taken from data)
    XSec_pmu_p1.at(j)->Add(background_pmu_p1.at(j), -1);
    XSec_pmu_p1.at(j)->Divide(xsec_mom_eff_p1.at(j));
    
    XSec_pmu_m1.at(j)->Add(pmu_numu_cc_reco_histo); // Nominal (taken from data)
    XSec_pmu_m1.at(j)->Add(background_pmu_m1.at(j), -1);
    XSec_pmu_m1.at(j)->Divide(xsec_mom_eff_m1.at(j));
  }
  
  
  
  // Now finally calculate the percental difference between sigma and sigma^+-
  // We want (sigma-sigma^+)/sigma   and   (sigma-sigma^-)/sigma
  
  XSec_pmu_percDiff_p1.resize(functionsName->size());
  XSec_pmu_percDiff_m1.resize(functionsName->size());
  for (int j=0; j<nFunc; j++) {
    
    nameBase = "xsec_pmu_percdiff_";
    XSec_pmu_percDiff_p1.at(j) = new TH1D(nameBase+functionsName->at(j)+"_Plus1Sigma",  ";p_{#mu};Percental Difference (%)", 10, 0, 2);
    XSec_pmu_percDiff_m1.at(j) = new TH1D(nameBase+functionsName->at(j)+"_Minus1Sigma", ";p_{#mu};Percental Difference (%)", 10, 0, 2);
  }
  
  for (int j=0; j<nFunc; j++) {
    
    XSec_pmu_percDiff_p1.at(j)->Add(XSec_pmu_nominal);
    XSec_pmu_percDiff_p1.at(j)->Add(XSec_pmu_p1.at(j), -1);
    XSec_pmu_percDiff_p1.at(j)->Divide(XSec_pmu_nominal);
    XSec_pmu_percDiff_p1.at(j)->Scale(100.);
    
    XSec_pmu_percDiff_m1.at(j)->Add(XSec_pmu_nominal);
    XSec_pmu_percDiff_m1.at(j)->Add(XSec_pmu_m1.at(j), -1);
    XSec_pmu_percDiff_m1.at(j)->Divide(XSec_pmu_nominal);
    XSec_pmu_percDiff_m1.at(j)->Scale(100.);
  }
  

  
  // Loop over all the histograms and find the abs max value for each parameter
  double maxima[4];
  for (int j=0; j<nFunc; j++) {
    //XSec_pmu_percDiff_p1.at(j)->GetXaxis()->SetRange(2,8); // To calculate maximum
    //XSec_pmu_percDiff_m1.at(j)->GetXaxis()->SetRange(2,8); // in a given range only
    
    maxima[0] = abs(XSec_pmu_percDiff_p1.at(j)->GetMaximum());
    maxima[1] = abs(XSec_pmu_percDiff_p1.at(j)->GetMinimum());
    maxima[2] = abs(XSec_pmu_percDiff_m1.at(j)->GetMaximum());
    maxima[3] = abs(XSec_pmu_percDiff_m1.at(j)->GetMinimum());
    
    double thisMax = *std::max_element(maxima, maxima+4);
    std::cout << "Func: " << GetLegendName(functionsName->at(j)) << "  Max: " << thisMax << std::endl;
    
  }

  
  
  
  //*********************************
  //
  // Save XSec to file
  //
  //*********************************
  
  TFile f = TFile("eventWeightXSec.root", "RECREATE");
  
  XSec_pmu_nominal->Write();
  for (unsigned int function = 0; function < pmu_numu_cc_reco_histo_p1.size(); function++) {
    XSec_pmu_p1.at(function)->Write();
    XSec_pmu_m1.at(function)->Write();
  }

  
  for (unsigned int function = 0; function < pmu_numu_cc_reco_histo_p1.size(); function++) {
    XSec_pmu_percDiff_p1.at(function)->Write();
    XSec_pmu_percDiff_m1.at(function)->Write();
  }
  
  std::cout << "XSex percental difference histograms saved to eventWeightXSec.root" << std::endl;

}





//________________________________________
void EvtWgtAnalysis::MakeXsecDiffPlots() {

  if (xsec_mom_truth_p1.size() == 0) {
    std::cout << "Calling EvtWgtAnalysis::CalculateXSecPercDifference() first." << std::endl;
    CalculateXSecPercDifference();
  }

  system("mkdir ./EvtWgtXsecDiffPlots");

  // Open LaTeX file
  bool makeLaTeX = true;
  ofstream latexFile;
  ofstream latexFile2;
  if (makeLaTeX){
    latexFile.open("./EvtWgtXsecDiffPlots/evtwgtXsecDiffPmu.tex");
    latexFile << "\\begin{table}[]" << endl;
    latexFile << "\\caption{evtwgtXsecDiffPmu}" << endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << endl;
    latexFile << "\\label{tab:}" << endl;
    latexFile << "\\centering" << endl;
    latexFile << "\\begin{tabular}{c c c c}" << endl;
    latexFile << "\\toprule" << endl;
    latexFile << "  &  Cross-section percental   &    $\\sqrt{S+B}/S$             &   $\\sqrt{S+B}/S$ \\\\"            << endl;
    latexFile << "  &  difference (\\%)          &    @ $6.6\\times 10^{20}$ POT  &   @ $0.5\\times 10^{20}$ POT \\\\" << endl;
    latexFile << "\\midrule" << endl;

    latexFile2.open("./EvtWgtXsecDiffPlots/evtwgtXsecDiffPmu_technote.tex");
    latexFile2 << "\\begin{table}[]" << endl;
    latexFile2 << "\\caption{evtwgtXsecDiffPmu}" << endl;
    latexFile2 << "\\captionsetup{format=hang,labelfont={sf,bf}}" << endl;
    latexFile2 << "\\label{tab:}" << endl;
    latexFile2 << "\\centering" << endl;
    latexFile2 << "\\begin{tabular}{c c c}" << endl;
    latexFile2 << "\\toprule" << endl;
    latexFile2 << "  &  \\multicolumn{2}{c}{Cross-section percental difference}  \\\\"            << endl;
    latexFile2 << "  &  at maximum point  &  at cross-section peak               \\\\"            << endl;
    latexFile2 << "\\midrule" << endl;
  }

  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);

  double loopMax;
  loopMax = pmu_numu_cc_reco_histo_p1.size();


  vector<TH1D*> XSec_pmu_percDiff_max;
  XSec_pmu_percDiff_max.resize(loopMax);
  TString nameBase;
  for (unsigned int function = 0; function < loopMax; function++) {
    nameBase = "xsec_pmu_percDiff_max_";
    XSec_pmu_percDiff_max.at(function) = new TH1D(nameBase+functionsName->at(function)+"_max",  ";p_{#mu};#propto#sigma^+ [arb.]", 10, 0, 2);
  }

  vector<bool> binMax; // This vector contains bin by bin info
                       // true: in this bin Plus1 is bigger
                       // flase: in this bin Minus1 is bigger

  for (unsigned int function = 0; function < loopMax; function++) {

    binMax.resize(XSec_pmu_percDiff_p1.at(function)->GetNbinsX());

    // First set the histograms so that I am plotting absolute values only
    for (int bin = 1; bin < XSec_pmu_percDiff_p1.at(function)->GetNbinsX()+1; bin++) {

      // First set the histograms so that I am plotting absolute values only
      XSec_pmu_percDiff_p1.at(function)->SetBinContent(bin,abs(XSec_pmu_percDiff_p1.at(function)->GetBinContent(bin)));
      XSec_pmu_percDiff_m1.at(function)->SetBinContent(bin,abs(XSec_pmu_percDiff_m1.at(function)->GetBinContent(bin)));
     
      // Now fill XSec_pmu_percDiff_max with the maximum between the two prevous histograms
      if (XSec_pmu_percDiff_p1.at(function)->GetBinContent(bin) > XSec_pmu_percDiff_m1.at(function)->GetBinContent(bin)){
        binMax[bin] = true;
        XSec_pmu_percDiff_max.at(function)->SetBinContent(bin, XSec_pmu_percDiff_p1.at(function)->GetBinContent(bin));
      } else {
        binMax[bin] = false;
        XSec_pmu_percDiff_max.at(function)->SetBinContent(bin, XSec_pmu_percDiff_m1.at(function)->GetBinContent(bin));
      }
    }

    // Find the bin where the cross-section in bigger
    int xsec_max_point_bin = pmu_numu_cc_reco_histo->GetMaximumBin();

    XSec_pmu_percDiff_max.at(function)->GetXaxis()->SetRange(2,8);
    double uncert = XSec_pmu_percDiff_max.at(function)->GetMaximum();
    double uncertAtXSecPeak = XSec_pmu_percDiff_max.at(function)->GetBinContent(xsec_max_point_bin);
    double maxBin = XSec_pmu_percDiff_max.at(function)->GetMaximumBin();
    XSec_pmu_percDiff_max.at(function)->GetXaxis()->SetRange(1,10);

    // Calulate sqrt(S+B)/S

    double POTscale = 5.3e19/2.43000e20;
    POTscale = 6.6e20/2.43000e20;
    std::cout << "******************* You are using 6.6e20 to scale. ===> Scale Factor: " << POTscale << std::endl;
    double S = pmu_numu_cc_reco_histo->Integral() * POTscale;
    double B = background_pmu_nominal->Integral() * POTscale;
    
    double POTscale_2 = 0.5e20/2.43000e20;
    double S_2 = pmu_numu_cc_reco_histo->Integral() * POTscale;
    double B_2 = background_pmu_nominal->Integral() * POTscale;

    TH1D * scoreFunction   = new TH1D("scoreFunctionAT6Point6e20POT",  ";p_{#mu};#sqrt{S+B}/S", 10, 0, 2);
    TH1D * scoreFunction_2 = new TH1D("scoreFunctionAT0Point5e20POT",  ";p_{#mu};#sqrt{S+B}/S", 10, 0, 2);    

    for (int bin = 1; bin < pmu_numu_cc_reco_histo->GetNbinsX()+1; bin++) {
      S = pmu_numu_cc_reco_histo->GetBinContent(bin) * POTscale;
      B = background_pmu_nominal->GetBinContent(bin) * POTscale;
      S_2 = pmu_numu_cc_reco_histo->GetBinContent(bin) * POTscale_2;
      B_2 = background_pmu_nominal->GetBinContent(bin) * POTscale_2;
      if (S == 0) {
        scoreFunction->SetBinContent(bin,0);
      }
      else
        scoreFunction->SetBinContent(bin,sqrt(S+B)/S*100);
      if (S_2 == 0) {
        scoreFunction_2->SetBinContent(bin,0);
      }
      else
        scoreFunction_2->SetBinContent(bin,sqrt(S_2+B_2)/S_2*100);
    }
    
    TFile tempFile = TFile("temp.root", "RECREATE");
    scoreFunction->Write();
    scoreFunction_2->Write();
    //t.Close();
    

    
    
/*
    if (binMax[maxBin])
      B = background_pmu_p1->Integral() * POTscale;
    else
      B = background_pmu_m1->Integral() * POTscale;

    std::cout << GetLegendName(functionsName->at(function)) << sqrt(S+B)/S*100. << std::endl;
*/


    // Make the table now
    if (makeLaTeX){
      latexFile << " $ " << GetLegendName(functionsName->at(function)) 
                << " $ & " << uncert 
                << " &   " << scoreFunction->GetBinContent(maxBin) 
                << " &   " << scoreFunction_2->GetBinContent(maxBin) << " \\\\ " << endl;

      latexFile2 << " $ " << GetLegendName(functionsName->at(function))
                 << " $ & " << uncert
                 << " & " << uncertAtXSecPeak << " \\\\ " << endl;
    }
   




    // Make the plots now
    TCanvas *c = new TCanvas("c", "canvas", 800, 700);

    c->SetBottomMargin(0.15);
    c->SetLeftMargin(0.15);
    c->SetLogy();

    XSec_pmu_percDiff_max.at(function)->SetMinimum(0.05);
    XSec_pmu_percDiff_max.at(function)->SetMaximum(100.);

    XSec_pmu_percDiff_max.at(function)->SetStats(0);   
    XSec_pmu_percDiff_max.at(function)->SetLineColor(kRed+1);
    XSec_pmu_percDiff_max.at(function)->GetXaxis()->SetTitle("p_{#mu} [GeV]");
    XSec_pmu_percDiff_max.at(function)->GetXaxis()->CenterTitle();
    //XSec_pmu_percDiff_max.at(function)->GetXaxis()->SetTitleSize(25);
    //XSec_pmu_percDiff_max.at(function)->GetXaxis()->SetTitleFont(43);
    XSec_pmu_percDiff_max.at(function)->GetXaxis()->SetTitleOffset(1.1);
    XSec_pmu_percDiff_max.at(function)->GetYaxis()->SetTitle("Percentage [%]");//("|Percental difference| [%]");
    XSec_pmu_percDiff_max.at(function)->GetYaxis()->CenterTitle();
    //XSec_pmu_percDiff_max.at(function)->GetYaxis()->SetTitleSize(25);
    //XSec_pmu_percDiff_max.at(function)->GetYaxis()->SetTitleFont(43);
    XSec_pmu_percDiff_max.at(function)->GetYaxis()->SetTitleOffset(1.1);
    XSec_pmu_percDiff_max.at(function)->GetYaxis()->SetNdivisions(516);
    XSec_pmu_percDiff_max.at(function)->Draw();
    
    scoreFunction->SetLineColor(9);
    scoreFunction->Draw("same");

    scoreFunction_2->SetLineColor(8);
    scoreFunction_2->Draw("same");


/*
    XSec_pmu_percDiff_p1.at(function)->SetStats(0);          
    XSec_pmu_percDiff_m1.at(function)->SetStats(0);          

    XSec_pmu_percDiff_p1.at(function)->SetLineColor(kRed+1);
    XSec_pmu_percDiff_m1.at(function)->SetLineColor(kGreen+1);

    XSec_pmu_percDiff_p1.at(function)->GetXaxis()->SetTitle("p_{#mu} [GeV]");
    XSec_pmu_percDiff_p1.at(function)->GetXaxis()->CenterTitle();
    //XSec_pmu_percDiff_p1.at(function)->GetXaxis()->SetTitleSize(25);
    //XSec_pmu_percDiff_p1.at(function)->GetXaxis()->SetTitleFont(43);
    XSec_pmu_percDiff_p1.at(function)->GetXaxis()->SetTitleOffset(1.0);

    XSec_pmu_percDiff_p1.at(function)->GetYaxis()->SetTitle("|Percental difference| [%]");
    XSec_pmu_percDiff_p1.at(function)->GetYaxis()->CenterTitle();
    //XSec_pmu_percDiff_p1.at(function)->GetYaxis()->SetTitleSize(25);
    //XSec_pmu_percDiff_p1.at(function)->GetYaxis()->SetTitleFont(43);
    XSec_pmu_percDiff_p1.at(function)->GetYaxis()->SetTitleOffset(1.0);
    XSec_pmu_percDiff_p1.at(function)->GetYaxis()->SetNdivisions(516); 
    
    XSec_pmu_percDiff_p1.at(function)->Draw();
    XSec_pmu_percDiff_m1.at(function)->Draw("same");
*/

      double x = 0.9;
      double y = 0.95;
      double size = 35;//28;
      int color = kRed+1;
      int font = 43;
      int align = 32;
      TString LaTeXname = GetLegendName(functionsName->at(function));
      TLatex *latex = new TLatex( x, y, LaTeXname);
      latex->SetNDC();
      latex->SetTextSize(size);
      latex->SetTextColor(color);
      latex->SetTextFont(font);
      latex->SetTextAlign(align);

      latex->Draw();


    TLegend* leg= new TLegend(0.1829574,0.648368,0.5639098,0.8486647,NULL,"brNDC");
 
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.02967359);

    leg->SetX1(0.179198);
    leg->SetX2(0.4924812);
    leg->SetY1(0.7091988);
    leg->SetY2(0.8531157);

    //leg->SetNColumns(3);

    leg->AddEntry(XSec_pmu_percDiff_max.at(function),  "|(#sigma-#sigma^{#pm})/#sigma|");
    leg->AddEntry(scoreFunction,                       "Stat. @ 6.6e20 POT");
    leg->AddEntry(scoreFunction_2,                     "Stat. @ 0.5e20 POT");
    //leg->SetNColumns(2);


    leg->Draw();
    //leg->SetNColumns(2);

    //leg->SetX1(0.1829574);
    //leg->SetX2(0.5639098);
    //leg->SetY1(0.648368);
    //leg->SetY2(0.8486647);

    TString SaveName;
    SaveName = "Pmu_"+functionsName->at(function);

    c->SetGridy();

    c->Print("./EvtWgtXsecDiffPlots/" + SaveName + ".C");
    c->Print("./EvtWgtXsecDiffPlots/" + SaveName + ".pdf");

    if (function == loopMax-1) { // If we are at the end of the loop save stats.
      latexFile2 << "Stat. @ 0.5e20 POT & - & " << scoreFunction_2->GetBinContent(xsec_max_point_bin) << " \\\\" << std::endl;
      latexFile2 << "Stat. @ 6.6e20 POT & - & " << scoreFunction->GetBinContent(xsec_max_point_bin) << " \\\\" << std::endl;
    }

  } // end loop functions

  if (makeLaTeX) {
    latexFile << "\\bottomrule" << endl;
    latexFile << "\\end{tabular}" << endl;
    latexFile << "\\end{table}" << endl;

    latexFile2 << "\\bottomrule" << endl;
    latexFile2 << "\\end{tabular}" << endl;
    latexFile2 << "\\end{table}" << endl;
  }
}








//_____________________________________
void EvtWgtAnalysis::MakeHistograms() {
  
  std::cout << "EvtWgtAnalysis::MakeHistograms starts" << std::endl;
  fflush(0);
  //TTree *mytemp = cflux;
  //AnaBNB anaBNBtree->mytemp);
  
  anaBNBtree = new AnaBNB(cflux);
  
  // It's important to active only the branches we are interested in
  // Some branches cause segmentation violation...
  anaBNBtree->fChain->SetBranchStatus("*",0);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_funcname",1);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_nfunc",1);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_weight",1);
  anaBNBtree->fChain->SetBranchStatus("evtwgt_nweight",1);
  anaBNBtree->fChain->SetBranchStatus("genie_primaries_pdg",1);
  anaBNBtree->fChain->SetBranchStatus("genie_P",1);
  anaBNBtree->fChain->SetBranchStatus("genie_Pz",1);
  anaBNBtree->fChain->SetBranchStatus("genie_status_code",1);
  anaBNBtree->fChain->SetBranchStatus("genie_no_primaries",1);
  anaBNBtree->fChain->SetBranchStatus("ccnc_truth",1);
  anaBNBtree->fChain->SetBranchStatus("nuPDG_truth",1);
  anaBNBtree->fChain->SetBranchStatus("mode_truth",1);
  
  
  Long64_t nflux = cflux->GetEntries();
  Long64_t treeentries = cflux->GetTree()->GetEntries();
  std::cout << "Total number of entries in the first tree: " << treeentries << std::endl;
  std::cout << "Total number of entries: " << nflux << std::endl;
  for (Long64_t i=0; i < nflux; ++i ) {
    cflux->GetEntry(i);
    
    
    if(treeNumber != cflux->GetTreeNumber()) {
      treeNumber = cflux->GetTreeNumber();
      std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
    }
    if (i % 1000 == 0) cout << "On entry " << i << endl;
    if(i == 0) {
      InstantiateHistograms(anaBNBtree->evtwgt_nfunc, anaBNBtree->evtwgt_funcname);
    }
    
    if(i == 0) Nfunctions = anaBNBtree->evtwgt_nfunc;
    if(i == 0) std::cout << "Number of functions in the first tree: " << Nfunctions << std::endl;
    if(Nfunctions != anaBNBtree->evtwgt_nfunc) cout<<"############################################################"<<endl;
    
    // Select cc/nc and neutrino type
    if (anaBNBtree->ccnc_truth[0] == 0 && anaBNBtree->nuPDG_truth[0] == 14) {  // CC and numu
      Nnumucc++;
      for (int f = 0; f < anaBNBtree->evtwgt_nfunc; f++) {
        if (anaBNBtree->evtwgt_funcname->at(f).find("genie_IntraNukeNmfp_Genie"))
          weightNmfp->Fill((anaBNBtree->evtwgt_weight->at(f)).at(0));
      }
      if (anaBNBtree->mode_truth[0] == 3) Ncoh++;
      if (anaBNBtree->mode_truth[0] == 2) Ndis++;
      // Loop over GENIE particles
      for (int GENIEparticle = 0; GENIEparticle < anaBNBtree->genie_no_primaries; GENIEparticle++) {
        // Look for muon
        if (anaBNBtree->genie_primaries_pdg[GENIEparticle] == 14-1 && anaBNBtree->genie_status_code[GENIEparticle]==1) {
          // Loop over the reweighting functions
          for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
            if (anaBNBtree->evtwgt_funcname->at(function) == "genie_DISAth_Genie") temp->Fill((anaBNBtree->evtwgt_weight->at(function)).at(0));
            histoPmuVec[function]->Fill(anaBNBtree->genie_P[GENIEparticle]);
            histoPmuVec_p1[function]->Fill(anaBNBtree->genie_P[GENIEparticle], (anaBNBtree->evtwgt_weight->at(function)).at(0));
            histoPmuVec_m1[function]->Fill(anaBNBtree->genie_P[GENIEparticle], (anaBNBtree->evtwgt_weight->at(function)).at(1));
            double cosAngle = anaBNBtree->genie_Pz[GENIEparticle]/anaBNBtree->genie_P[GENIEparticle];
            histoCosThetaMuVec[function]->Fill(cosAngle);
            histoCosThetaMuVec_p1[function]->Fill(cosAngle, (anaBNBtree->evtwgt_weight->at(function)).at(0));
            histoCosThetaMuVec_m1[function]->Fill(cosAngle, (anaBNBtree->evtwgt_weight->at(function)).at(1));
          } // end loop over reweighting functions
        } // end muon if
      } // end loop GENIE particle
    } // end cc/nc and neutrino type selection
  } // end of loop over entries
  std::cout << "Loop completed. Saving to eventWeightHistograms.root now." << std::endl;
 
  std::cout << "Number of numu cc events:                               " << Nnumucc << std::endl;                  
  std::cout << "  > Number of Coherent Production interactions:         " << Ncoh    << std::endl;
  std::cout << "  > Number of Deep Inelastinc Scattering  interactions: " << Ndis    << std::endl; 

 
  // Save to file
  TFile f = TFile("eventWeightHistograms.root", "RECREATE");
  temp->Write();
  weightNmfp->Write();
  for (unsigned int j=0; j<histoPmuVec.size(); j++) {
    histoPmuVec[j]->Write();
    histoPmuVec_p1[j]->Write();
    histoPmuVec_m1[j]->Write();
    histoCosThetaMuVec[j]->Write();
    histoCosThetaMuVec_p1[j]->Write();
    histoCosThetaMuVec_m1[j]->Write();
  }
  f.Close();
  cout << "Histogram correctly written to file eventWeightHistograms.root" << endl;
  
}

//___________________________________________________________________________________________
void EvtWgtAnalysis::InstantiateHistogramsEfficiency(Int_t nFunc, vector<string> *funcName) {
  
  // Check the function names vector is of the right dimension
  unsigned int temp = nFunc;
  if (funcName->size() != temp) {
    std::cout << "FATAL: funcName->size() != nFunc" << std::endl;
    //exit(0);
  }
  
  // Save the function names
  functionsName = funcName;
  
  // Resize the vectors
  all_evts_p1.resize(funcName->size());
  all_evts_m1.resize(funcName->size());
  sel_evts_p1.resize(funcName->size());
  sel_evts_m1.resize(funcName->size());

  xsec_mom_truth_p1.resize(funcName->size());
  xsec_mom_truth_m1.resize(funcName->size());
  xsec_theta_truth_p1.resize(funcName->size());
  xsec_theta_truth_m1.resize(funcName->size());
  xsec_mom_data_p1.resize(funcName->size());
  xsec_mom_data_m1.resize(funcName->size());
  xsec_mom_reco_data_p1.resize(funcName->size());
  xsec_mom_reco_data_m1.resize(funcName->size());
  xsec_theta_data_p1.resize(funcName->size());
  xsec_theta_data_m1.resize(funcName->size());
  xsec_theta_reco_data_p1.resize(funcName->size());
  xsec_theta_reco_data_m1.resize(funcName->size());
  xsec_mom_eff_p1.resize(funcName->size());
  xsec_mom_eff_m1.resize(funcName->size());
  xsec_mom_reco_eff_p1.resize(funcName->size());
  xsec_mom_reco_eff_m1.resize(funcName->size());
  xsec_theta_eff_p1.resize(funcName->size());
  xsec_theta_eff_m1.resize(funcName->size());
  xsec_theta_reco_eff_p1.resize(funcName->size());
  xsec_theta_reco_eff_m1.resize(funcName->size());
  
  pmu_numu_cc_reco_histo_p1.resize(funcName->size());
  pmu_numu_cc_reco_histo_m1.resize(funcName->size());
  costhetamu_numu_cc_reco_histo_p1.resize(funcName->size());
  costhetamu_numu_cc_reco_histo_m1.resize(funcName->size());
  pmu_anumu_cc_reco_histo_p1.resize(funcName->size());
  pmu_anumu_cc_reco_histo_m1.resize(funcName->size());
  costhetamu_anumu_cc_reco_histo_p1.resize(funcName->size());
  costhetamu_anumu_cc_reco_histo_m1.resize(funcName->size());
  pmu_nue_cc_reco_histo_p1.resize(funcName->size());
  pmu_nue_cc_reco_histo_m1.resize(funcName->size());
  costhetamu_nue_cc_reco_histo_p1.resize(funcName->size());
  costhetamu_nue_cc_reco_histo_m1.resize(funcName->size());
  pmu_nc_reco_histo_p1.resize(funcName->size());
  pmu_nc_reco_histo_m1.resize(funcName->size());
  costhetamu_nc_reco_histo_p1.resize(funcName->size());
  costhetamu_nc_reco_histo_m1.resize(funcName->size());
  
  // Initialize counters
  for ( int j=0; j<nFunc; j++) {
    all_evts_p1[j] = 0.;
    all_evts_m1[j] = 0.;
    sel_evts_p1[j] = 0.;
    sel_evts_m1[j] = 0.;
  }
  
  // Instantiate the histograms
  TString histoTitle;
  for ( int j=0; j<nFunc; j++) {
    histoTitle = "xsec_mom_eff_";
    xsec_mom_eff_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    xsec_mom_eff_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    
    histoTitle = "xsec_theta_eff_";
    xsec_theta_eff_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    xsec_theta_eff_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    
    histoTitle = "xsec_mom_reco_eff_";
    xsec_mom_reco_eff_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    xsec_mom_reco_eff_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    
    histoTitle = "xsec_theta_reco_eff_";
    xsec_theta_reco_eff_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    xsec_theta_reco_eff_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    
    histoTitle = "xsec_mom_truth_";
    xsec_mom_truth_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    xsec_mom_truth_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    
    histoTitle = "xsec_theta_truth_";
    xsec_theta_truth_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    xsec_theta_truth_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    
    histoTitle = "xsec_mom_data_";
    xsec_mom_data_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    xsec_mom_data_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    
    histoTitle = "xsec_mom_reco_data_";
    xsec_mom_reco_data_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    xsec_mom_reco_data_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Efficiency", 10, 0, 2);
    
    histoTitle = "xsec_theta_data_";
    xsec_theta_data_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    xsec_theta_data_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    
    histoTitle = "xsec_theta_reco_data_";
    xsec_theta_reco_data_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    xsec_theta_reco_data_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Efficiency", 10, -1, 1);
    
    
    
    histoTitle = "pmu_numu_cc_reco_histo_";
    pmu_numu_cc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    pmu_numu_cc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    
    histoTitle = "costhetamu_numu_cc_reco_histo_";
    costhetamu_numu_cc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    costhetamu_numu_cc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    
    histoTitle = "pmu_anumu_cc_reco_histo_";
    pmu_anumu_cc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    pmu_anumu_cc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    
    histoTitle = "costhetamu_anumu_cc_reco_histo_";
    costhetamu_anumu_cc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    costhetamu_anumu_cc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    
    histoTitle = "pmu_nue_cc_reco_histo_";
    pmu_nue_cc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    pmu_nue_cc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    
    histoTitle = "costhetamu_nue_cc_reco_histo_";
    costhetamu_nue_cc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    costhetamu_nue_cc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    
    histoTitle = "pmu_nc_reco_histo_";
    pmu_nc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    pmu_nc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Events", 10, 0, 2);
    
    histoTitle = "costhetamu_nc_reco_histo_";
    costhetamu_nc_reco_histo_p1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    costhetamu_nc_reco_histo_m1.at(j) = new TH1D(histoTitle+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Events", 10, -1, 1);
    
  }
  
}

//_________________________________________________________________________________
void EvtWgtAnalysis::InstantiateHistograms(Int_t nFunc, vector<string> *funcName) {
  
  // Check the function names vector is of the right dimension
  unsigned int temp = nFunc;
  if (funcName->size() != temp) {
    std::cout << "FATAL: funcName->size() != nFunc" << std::endl;
    exit(0);
  }
  
  // Save the function names
  functionsName = funcName;
  
  // Instantiate the histograms
  histoPmuVec.resize(funcName->size());
  histoPmuVec_p1.resize(funcName->size());
  histoPmuVec_m1.resize(funcName->size());
  histoCosThetaMuVec.resize(funcName->size());
  histoCosThetaMuVec_p1.resize(funcName->size());
  histoCosThetaMuVec_m1.resize(funcName->size());
  
  TString histoTitleP = "Pmu_";
  TString histoTitleCos = "CosThetaMu_";
  for ( int j=0; j<nFunc; j++) {
    // Pmu
    histoPmuVec.at(j) = new TH1D(histoTitleP+funcName->at(j), ";p_{#mu} [GeV];Entries per bin", 25, 0,2);
    histoPmuVec_p1.at(j) = new TH1D(histoTitleP+funcName->at(j)+"_Plus1Sigma", ";p_{#mu} [GeV];Entries per bin", 25, 0,2);
    histoPmuVec_m1.at(j) = new TH1D(histoTitleP+funcName->at(j)+"_Minus1Sigma", ";p_{#mu} [GeV];Entries per bin", 25, 0,2);
    // CosThetaMu
    histoCosThetaMuVec.at(j) = new TH1D(histoTitleCos+funcName->at(j), ";cos#theta_{#mu};Entries per bin", 25, -1, 1);
    histoCosThetaMuVec_p1.at(j) = new TH1D(histoTitleCos+funcName->at(j)+"_Plus1Sigma", ";cos#theta_{#mu};Entries per bin", 25, -1, 1);
    histoCosThetaMuVec_m1.at(j) = new TH1D(histoTitleCos+funcName->at(j)+"_Minus1Sigma", ";cos#theta_{#mu};Entries per bin", 25, -1, 1);
    
    
    
  }
}






//______________________________________________________
void EvtWgtAnalysis::MakeBackgroundPlots(int variable) {
  
  // Pmu: variable == 0
  // CosThetaMu: variable == 1
  
  if (pmu_numu_cc_reco_histo_p1.size() == 0 || costhetamu_numu_cc_reco_histo_p1.size() == 0 ) {
    cout << "Calling CalcEfficiency() first."  << endl;
    CalcEfficiency();
  }
  
  bool makeLatex = true;
  
  // Create output directory
  system("mkdir ./EvtWgtBackgroundPlots");
  system("mkdir ./EvtWgtBackgroundPlots_reducedLegend"); 
 
  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);
  
  
  double scaleFactor = 5.3e19/1.22e20;
  
  cout << "numu:  " << pmu_numu_cc_reco_histo->Integral() * scaleFactor << endl;
  cout << "anumu: " << pmu_anumu_cc_reco_histo->Integral() * scaleFactor << endl;
  cout << "nue:   " << pmu_nue_cc_reco_histo->Integral() * scaleFactor << endl;
  cout << "nc:    " << pmu_nc_reco_histo->Integral() * scaleFactor << endl;
  
  ofstream latexFile;
  double numu_nominal = pmu_numu_cc_reco_histo->Integral();
  double anumu_nominal = pmu_anumu_cc_reco_histo->Integral();
  double nue_nominal = pmu_nue_cc_reco_histo->Integral();
  double nc_nominal = pmu_nc_reco_histo->Integral();

  if(makeLatex) {
    if (variable == 0) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundPmu.tex");
    if (variable == 1) latexFile.open("./EvtWgtBackgroundPlots/evtwgtBackgroundCosThetaMu.tex");
    latexFile << "\\begin{table}[]" << endl;
    latexFile << "\\caption{}" << endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << endl;
    latexFile << "\\label{tab:}" << endl;
    latexFile << "\\centering" << endl;
    latexFile << "\\begin{tabular}{c   c c  c c  c c  c c}" << endl;
    latexFile << "\\toprule" << endl;
    latexFile << "  &  \\multicolumn{2}{ c }{$\\nu_\\mu$ CC} &  \\multicolumn{2}{ c }{$\\bar{\\nu}_\\mu$ CC}   & \\multicolumn{2}{ c }{$\\nu_e$, $\\bar{\\nu}_e$ CC}   &  \\multicolumn{2}{ c }{NC} \\\\" << endl;
    latexFile << "  &  Events   &  Diff. (\\%) &  Events  &  Diff. (\\%) & Events  &  Diff. (\\%) &  Events &  Diff. (\\%) \\\\" << endl;
    latexFile << "\\midrule" << endl;
    latexFile << "Nominal & " << pmu_numu_cc_reco_histo->Integral() << " & 0 "
    << " & " << pmu_anumu_cc_reco_histo->Integral() << " & 0 "
    << " & " << pmu_nue_cc_reco_histo->Integral() << " & 0 "
    << " & " << pmu_nc_reco_histo->Integral() << " & 0 " << "\\\\" << endl;
  }
  
  
  
  
  
  TCanvas *c = new TCanvas("c", "canvas", 800, 800);
  
  double loopMax;
  loopMax = pmu_numu_cc_reco_histo_p1.size();
  
  for (unsigned int function = 0; function < loopMax; function++) {
    
    TString SaveName;
    if (variable == 0) SaveName = "Pmu_"+functionsName->at(function);
    if (variable == 1) SaveName = "CosThetaMu_"+functionsName->at(function);
    
    c->SetLogy();                                 // Log scale
    //c->SetGridy();                                // Horizontal grid
    
    
    TH1F* histo_numu;
    TH1F* histo_numu_p1;
    TH1F* histo_numu_m1;
    TH1F* histo_anumu;
    TH1F* histo_anumu_p1;
    TH1F* histo_anumu_m1;
    TH1F* histo_nue;
    TH1F* histo_nue_p1;
    TH1F* histo_nue_m1;
    TH1F* histo_nc;
    TH1F* histo_nc_p1;
    TH1F* histo_nc_m1;
    
    if (variable == 0) {
      histo_numu       = (TH1F*)pmu_numu_cc_reco_histo->Clone("histo_numu");
      histo_numu_p1    = (TH1F*)pmu_numu_cc_reco_histo_p1[function]->Clone("histo_numu_p1");
      histo_numu_m1    = (TH1F*)pmu_numu_cc_reco_histo_m1[function]->Clone("histo_numu_m1");
      
      histo_anumu       = (TH1F*)pmu_anumu_cc_reco_histo->Clone("histo_anumu");
      histo_anumu_p1    = (TH1F*)pmu_anumu_cc_reco_histo_p1[function]->Clone("histo_anumu_p1");
      histo_anumu_m1    = (TH1F*)pmu_anumu_cc_reco_histo_m1[function]->Clone("histo_anumu_m1");
      
      histo_nue       = (TH1F*)pmu_nue_cc_reco_histo->Clone("histo_nue");
      histo_nue_p1    = (TH1F*)pmu_nue_cc_reco_histo_p1[function]->Clone("histo_nue_p1");
      histo_nue_m1    = (TH1F*)pmu_nue_cc_reco_histo_m1[function]->Clone("histo_nue_m1");
      
      histo_nc       = (TH1F*)pmu_nc_reco_histo->Clone("histo_nc");
      histo_nc_p1    = (TH1F*)pmu_nc_reco_histo_p1[function]->Clone("histo_nc_p1");
      histo_nc_m1    = (TH1F*)pmu_nc_reco_histo_m1[function]->Clone("histo_nc_m1");
    }
    else if (variable == 1) {
      histo_numu       = (TH1F*)costhetamu_numu_cc_reco_histo->Clone("histo_numu");
      histo_numu_p1    = (TH1F*)costhetamu_numu_cc_reco_histo_p1[function]->Clone("histo_numu_p1");
      histo_numu_m1    = (TH1F*)costhetamu_numu_cc_reco_histo_m1[function]->Clone("histo_numu_m1");
      
      histo_anumu       = (TH1F*)costhetamu_anumu_cc_reco_histo->Clone("histo_anumu");
      histo_anumu_p1    = (TH1F*)costhetamu_anumu_cc_reco_histo_p1[function]->Clone("histo_anumu_p1");
      histo_anumu_m1    = (TH1F*)costhetamu_anumu_cc_reco_histo_m1[function]->Clone("histo_anumu_m1");
      
      histo_nue       = (TH1F*)costhetamu_nue_cc_reco_histo->Clone("histo_nue");
      histo_nue_p1    = (TH1F*)costhetamu_nue_cc_reco_histo_p1[function]->Clone("histo_nue_p1");
      histo_nue_m1    = (TH1F*)costhetamu_nue_cc_reco_histo_m1[function]->Clone("histo_nue_m1");
      
      histo_nc       = (TH1F*)costhetamu_nc_reco_histo->Clone("histo_nc");
      histo_nc_p1    = (TH1F*)costhetamu_nc_reco_histo_p1[function]->Clone("histo_nc_p1");
      histo_nc_m1    = (TH1F*)costhetamu_nc_reco_histo_m1[function]->Clone("histo_nc_m1");
    } else {
      cout << "Invalid option. Exit." << endl;
      exit(0);
    }
    
    
    // Settings
    histo_numu->SetStats(0);          // No statistics on upper plot
    if (variable == 1) histo_numu->SetMinimum(1);
    if (variable == 1) histo_numu->SetMaximum(1e5);

    if (variable == 0) histo_numu->GetXaxis()->SetTitle("Reconstructed p_{#mu} [GeV]");
    if (variable == 1) histo_numu->GetXaxis()->SetTitle("Reconstructed cos#theta_{#mu}");
    histo_numu->GetXaxis()->CenterTitle();
    histo_numu->GetXaxis()->SetTitleSize(25);
    histo_numu->GetXaxis()->SetTitleFont(43);
    histo_numu->GetXaxis()->SetTitleOffset(1.45);
    
    histo_numu->GetYaxis()->SetTitle("Events");
    histo_numu->GetYaxis()->CenterTitle();
    histo_numu->GetYaxis()->SetTitleSize(25);
    histo_numu->GetYaxis()->SetTitleFont(43);
    histo_numu->GetYaxis()->SetTitleOffset(1.55);
    
    double lineWidth = 1;
    
    
    // Nominal
    histo_numu  ->Draw();
    histo_anumu ->Draw("same");
    histo_nue   ->Draw("same");
    histo_nc    ->Draw("same");
    
    histo_numu  ->SetLineColor(kViolet+1);
    histo_anumu ->SetLineColor(kBlue+1);
    histo_nue   ->SetLineColor(kOrange+1);
    histo_nc    ->SetLineColor(kGreen+1);
    
    histo_numu  ->SetLineWidth(lineWidth+1);
    histo_anumu ->SetLineWidth(lineWidth+1);
    histo_nue   ->SetLineWidth(lineWidth+1);
    histo_nc    ->SetLineWidth(lineWidth+1);
    
    // +- 1 sigma
    histo_numu_p1  ->Draw("same");
    histo_anumu_p1 ->Draw("same");
    histo_nue_p1   ->Draw("same");
    histo_nc_p1    ->Draw("same");
    
    histo_numu_m1  ->Draw("same");
    histo_anumu_m1 ->Draw("same");
    histo_nue_m1   ->Draw("same");
    histo_nc_m1    ->Draw("same");
    
    histo_numu_p1  ->SetLineColor(kViolet+1);
    histo_anumu_p1 ->SetLineColor(kBlue+1);
    histo_nue_p1   ->SetLineColor(kOrange+1);
    histo_nc_p1    ->SetLineColor(kGreen+1);
    
    histo_numu_m1  ->SetLineColor(kViolet+1);
    histo_anumu_m1 ->SetLineColor(kBlue+1);
    histo_nue_m1   ->SetLineColor(kOrange+1);
    histo_nc_m1    ->SetLineColor(kGreen+1);
    
    histo_numu_p1  ->SetLineStyle(7);
    histo_anumu_p1 ->SetLineStyle(7);
    histo_nue_p1   ->SetLineStyle(7);
    histo_nc_p1    ->SetLineStyle(7);
    
    histo_numu_m1  ->SetLineStyle(3);
    histo_anumu_m1 ->SetLineStyle(3);
    histo_nue_m1   ->SetLineStyle(3);
    histo_nc_m1    ->SetLineStyle(3);
    
    histo_numu_p1  ->SetLineWidth(lineWidth);
    histo_anumu_p1 ->SetLineWidth(lineWidth);
    histo_nue_p1   ->SetLineWidth(lineWidth);
    histo_nc_p1    ->SetLineWidth(lineWidth);
    
    histo_numu_m1  ->SetLineWidth(lineWidth);
    histo_anumu_m1 ->SetLineWidth(lineWidth);
    histo_nue_m1   ->SetLineWidth(lineWidth);
    histo_nc_m1    ->SetLineWidth(lineWidth);
    
    
    // Legend
    TLegend* leg;
    if (variable == 0) leg = new TLegend(0.566416,0.5535484,0.8822055,0.8825806,NULL,"brNDC");
    if (variable == 1) leg = new TLegend(0.1679198,0.5367742,0.4837093,0.8658065,NULL,"brNDC");
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    
    leg->AddEntry(histo_numu,  "#nu_{#mu} CC (nominal)");
    leg->AddEntry(histo_anumu, "#bar{#nu}_{#mu} CC (nominal)");
    leg->AddEntry(histo_nue,   "#nu_{e} CC (nominal)");
    leg->AddEntry(histo_nc,    "NC all flavours (nominal)");
    
    leg->AddEntry(histo_numu_p1,  "#nu_{#mu} CC (" + GetLegendName(functionsName->at(function)) + " + 1 #sigma)");
    leg->AddEntry(histo_anumu_p1, "#bar{#nu}_{#mu} CC (" + GetLegendName(functionsName->at(function)) + " + 1 #sigma)");
    leg->AddEntry(histo_nue_p1,   "#nu_{e} CC (" + GetLegendName(functionsName->at(function)) + " + 1 #sigma)");
    leg->AddEntry(histo_nc_p1,    "NC all flavours (" + GetLegendName(functionsName->at(function)) + " + 1 #sigma)");
    
    leg->AddEntry(histo_numu_m1,  "#nu_{#mu} CC (" + GetLegendName(functionsName->at(function)) + " - 1 #sigma)");
    leg->AddEntry(histo_anumu_m1, "#bar{#nu}_{#mu} CC (" + GetLegendName(functionsName->at(function)) + " - 1 #sigma)");
    leg->AddEntry(histo_nue_m1,   "#nu_{e} CC (" + GetLegendName(functionsName->at(function)) + " - 1 #sigma)");
    leg->AddEntry(histo_nc_m1,    "NC all flavours (" + GetLegendName(functionsName->at(function)) + " - 1 #sigma)");
    
    leg->Draw();
    
    
    c->Print("./EvtWgtBackgroundPlots/" + SaveName + ".C");
    c->Print("./EvtWgtBackgroundPlots/" + SaveName + ".pdf");
    

    // Reduced Legend Plots
    histo_numu->GetXaxis()->SetTitleOffset(1.10);
    histo_numu->GetYaxis()->SetTitleOffset(1.30);
    histo_numu->GetXaxis()->SetTitleSize(30);
    histo_numu->GetYaxis()->SetTitleSize(30);

    histo_numu  ->Draw();
    histo_anumu ->Draw("same");
    histo_nue   ->Draw("same");
    histo_nc    ->Draw("same");
    histo_numu_p1  ->Draw("same");
    histo_anumu_p1 ->Draw("same");
    histo_nue_p1   ->Draw("same");
    histo_nc_p1    ->Draw("same");

    histo_numu_m1  ->Draw("same");
    histo_anumu_m1 ->Draw("same");
    histo_nue_m1   ->Draw("same");
    histo_nc_m1    ->Draw("same");

    TLegend* leg2;
    if (variable == 0) leg2 = new TLegend(0.566416,0.5535484,0.8822055,0.8825806,NULL,"brNDC");
    if (variable == 1) leg2 = new TLegend(0.1679198,0.5367742,0.4837093,0.8658065,NULL,"brNDC");
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);

    leg2->AddEntry(histo_numu,  "#nu_{#mu} CC");
    leg2->AddEntry(histo_anumu, "#bar{#nu}_{#mu} CC");
    leg2->AddEntry(histo_nue,   "#nu_{e}, #bar{#nu}_{e} CC");
    leg2->AddEntry(histo_nc,    "NC all flavours");

    leg2->Draw();

    c->Print("./EvtWgtBackgroundPlots_reducedLegend/" + SaveName + ".C");
    c->Print("./EvtWgtBackgroundPlots_reducedLegend/" + SaveName + ".pdf");




    double numu_p1;
    double numu_m1;
    double anumu_p1;
    double anumu_m1;
    double nue_p1;
    double nue_m1;
    double nc_p1;
    double nc_m1;
    
    if (variable == 0) {
      numu_p1  = pmu_numu_cc_reco_histo_p1[function]  ->Integral();
      numu_m1  = pmu_numu_cc_reco_histo_m1[function]  ->Integral();
      anumu_p1 = pmu_anumu_cc_reco_histo_p1[function] ->Integral();
      anumu_m1 = pmu_anumu_cc_reco_histo_m1[function] ->Integral();
      nue_p1   = pmu_nue_cc_reco_histo_p1[function]   ->Integral();
      nue_m1   = pmu_nue_cc_reco_histo_m1[function]   ->Integral();
      nc_p1    = pmu_nc_reco_histo_p1[function]       ->Integral();
      nc_m1    = pmu_nc_reco_histo_m1[function]       ->Integral();
    } else if (variable == 1) {
      numu_p1  = costhetamu_numu_cc_reco_histo_p1[function]  ->Integral();
      numu_m1  = costhetamu_numu_cc_reco_histo_m1[function]  ->Integral();
      anumu_p1 = costhetamu_anumu_cc_reco_histo_p1[function] ->Integral();
      anumu_m1 = costhetamu_anumu_cc_reco_histo_m1[function] ->Integral();
      nue_p1   = costhetamu_nue_cc_reco_histo_p1[function]   ->Integral();
      nue_m1   = costhetamu_nue_cc_reco_histo_m1[function]   ->Integral();
      nc_p1    = costhetamu_nc_reco_histo_p1[function]       ->Integral();
      nc_m1    = costhetamu_nc_reco_histo_m1[function]       ->Integral();
    }
    
    if(makeLatex) {
      double numu_diff_p1 = (numu_p1 -numu_nominal)/numu_nominal * 100.;
      double numu_diff_m1 = (numu_p1 -numu_nominal)/numu_nominal * 100.;
      
      double anumu_diff_p1 = (anumu_p1 -anumu_nominal)/anumu_nominal * 100.;
      double anumu_diff_m1 = (anumu_m1 -anumu_nominal)/anumu_nominal * 100.;

      double nue_diff_p1 = (nue_p1 -nue_nominal)/nue_nominal * 100.;
      double nue_diff_m1 = (nue_m1 -nue_nominal)/nue_nominal * 100.;
      
      double nc_diff_p1 = (nc_p1 -nc_nominal)/nc_nominal * 100.;
      double nc_diff_m1 = (nc_p1 -nc_nominal)/nc_nominal * 100.;


      latexFile << "\\midrule" << endl;
      latexFile << "$" << GetLegendName(functionsName->at(function)) << " + 1\\sigma$ "
      << " & " << std::setprecision(5) << numu_p1  << " & " << std::setprecision(2) << numu_diff_p1
      << " & " << std::setprecision(5) << anumu_p1 << " & " << std::setprecision(2) << anumu_diff_p1
      << " & " << std::setprecision(5) << nue_p1   << " & " << std::setprecision(2) << nue_diff_p1
      << " & " << std::setprecision(5) << nc_p1    << " & " << std::setprecision(2) << nc_diff_p1 << "\\\\" << endl;
      
      latexFile << "$" << GetLegendName(functionsName->at(function)) << " - 1\\sigma$ "
      << " & " << std::setprecision(5) << numu_m1  << " & " << std::setprecision(2) << numu_diff_m1
      << " & " << std::setprecision(5) << anumu_m1 << " & " << std::setprecision(2) << anumu_diff_m1
      << " & " << std::setprecision(5) << nue_m1   << " & " << std::setprecision(2) << nue_diff_m1
      << " & " << std::setprecision(5) << nc_m1    << " & " << std::setprecision(2) << nc_diff_m1 << "\\\\" << endl;
    }
  }
  
  
  
  
  if(makeLatex) {
    latexFile << "\\bottomrule" << endl;
    latexFile << "\\end{tabular}" << endl;
    latexFile << "\\end{table}" << endl;
  }
  
  
  
}









//_____________________________________________________________
void EvtWgtAnalysis::MakePlots(bool normalised, int variable) {
  
  // Variable:
  // 0: efficiency - Pmu
  // 1: efficiency - CosThetaMu
  // 2: events - Pmu
  // 3: events - CosThetaMu
  
  bool makeLatex = true;
  
  if (xsec_mom_eff_p1.size() == 0 && (variable == 0 || variable == 1)) {
    cout << "Calling CalcEfficiency() first."  << endl;
    CalcEfficiency();
  }
  if (histoPmuVec.size() == 0 && (variable == 2 || variable == 3)) {
    cout << "Calling MakeHistograms() first."  << endl;
    MakeHistograms();
  }
  
  // Make a directory to store the plots
  //if (variable == 0) system("rm -rf ./EvtWgtEfficiencyPlots/Pmu* ");
  //if (variable == 1) system("rm -rf ./EvtWgtEfficiencyPlots/Cos* ");
  //if (variable == 2) system("rm -rf ./EvtWgtEventPlots/Pmu* ");
  //if (variable == 3) system("rm -rf ./EvtWgtEventPlots/CosThetaMu* ");
  if (variable == 0 || variable == 1)system("mkdir ./EvtWgtEfficiencyPlots");
  if (variable == 2 || variable == 3)system("mkdir ./EvtWgtEventPlots");
  
  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);
  
  // Opening a text file to write the integrals
  ofstream outfile;
  if (variable == 0) outfile.open("./EvtWgtEfficiencyPlots/IntegralsPmu.txt");
  if (variable == 1) outfile.open("./EvtWgtEfficiencyPlots/IntegralsCosThetaMu.txt");
  if (variable == 2) outfile.open("./EvtWgtEventPlots/IntegralsPmu.txt");
  if (variable == 3) outfile.open("./EvtWgtEventPlots/IntegralsCosThetaMu.txt");
  
  
  
  
  // Open LaTeX file to write the table
  ofstream latexFile;
  if(makeLatex) {
    if (variable == 0) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyPmu.tex");
    if (variable == 1) latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiencyCosThetaMu.tex");
    if (variable == 2) latexFile.open("./EvtWgtEventPlots/evtwgtEventPmu.tex");
    if (variable == 3) latexFile.open("./EvtWgtEventPlots/evtwgtEventCosThetaMu.tex");
    latexFile << "\\begin{table}[]" << endl;
    latexFile << "\\caption{}" << endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << endl;
    latexFile << "\\label{tab:}" << endl;
    latexFile << "\\centering" << endl;
    latexFile << "\\begin{tabular}{ccc}" << endl;
    latexFile << "\\toprule" << endl;
    if (variable == 0 || variable == 1) latexFile << "  &  Efficiency  &  Difference (\\%) \\\\" << endl;
    else latexFile << "  &  Integral  &  Difference (\\%) \\\\" << endl;
    latexFile << "\\midrule" << endl;
  }
  
  // Looping over all the functions
  double loopMax;
  if (variable == 0 || variable == 1) loopMax = xsec_mom_eff_p1.size();
  if (variable == 2 || variable == 3) loopMax = histoPmuVec.size();
  
  for (unsigned int function = 0; function < loopMax; function++) {
    
    TH1F *histoPmu;
    TH1F *histoPmu_p1;
    TH1F *histoPmu_m1;
    
    if (variable == 0) {  // Pmu
      histoPmu    = (TH1F*)xsec_mom_reco_eff->Clone("histoPmu");
      histoPmu_p1 = (TH1F*)xsec_mom_reco_eff_p1.at(function)->Clone("histoPmu_p1");
      histoPmu_m1 = (TH1F*)xsec_mom_reco_eff_m1.at(function)->Clone("histoPmu_m1");
    }
    else if (variable == 1) { // CosThetaMu
      histoPmu    = (TH1F*)xsec_theta_eff->Clone("histoPmu");
      histoPmu_p1 = (TH1F*)xsec_theta_eff_p1.at(function)->Clone("histoPmu_p1");
      histoPmu_m1 = (TH1F*)xsec_theta_eff_m1.at(function)->Clone("histoPmu_m1");
    }
    else if (variable == 2) { // Pmu
      histoPmu    = (TH1F*)histoPmuVec.at(function)->Clone("histoPmu");
      histoPmu_p1 = (TH1F*)histoPmuVec_p1.at(function)->Clone("histoPmu_p1");
      histoPmu_m1 = (TH1F*)histoPmuVec_m1.at(function)->Clone("histoPmu_m1");
    }
    else if (variable == 3) { // CosThetaMu
      histoPmu    = (TH1F*)histoCosThetaMuVec.at(function)->Clone("histoPmu");
      histoPmu_p1 = (TH1F*)histoCosThetaMuVec_p1.at(function)->Clone("histoPmu_p1");
      histoPmu_m1 = (TH1F*)histoCosThetaMuVec_m1.at(function)->Clone("histoPmu_m1");
    }
    else {
      cout << "Invalid option. Exit." << endl;
      exit(0);
    }
    TString SaveName;
    if(variable == 0 || variable == 2) SaveName = "Pmu_"+functionsName->at(function);
    if(variable == 1 || variable == 3) SaveName = "CosThetaMu_"+functionsName->at(function);
    TString LegName  = GetLegendName(functionsName->at(function));
    
    if(normalised) SaveName += "_normalised";
    
    /*
     if(makeLatex) {
     latexFile << "\\subsection{" << functionsName->at(function) << "}" << endl;
     latexFile << "\begin{figure}" << endl;
     latexFile << "\centering" << endl;
     latexFile << "\subfloat[][caption]" << endl;
     latexFile << "{\includegraphics[width=0.45\textwidth]{" << SaveName"}} \quad
     }
     */
    double histoPmu_Int = histoPmu->Integral();
    double histoPmu_p1_Int = histoPmu_p1->Integral();
    double histoPmu_m1_Int = histoPmu_m1->Integral();
    
    
    if (normalised) {
      histoPmu->Scale(1./histoPmu_Int);
      histoPmu_p1->Scale(1./histoPmu_p1_Int);
      histoPmu_m1->Scale(1./histoPmu_m1_Int);
    }
    
    // Calculate integrals
    outfile << functionsName->at(function) << endl;
    outfile << "Integral Nominal:    " << histoPmu->Integral() << endl;
    outfile << "Integral nominal_p1: " << histoPmu_p1->Integral() << endl;
    outfile << "Integral nominal_m1: " << histoPmu_m1->Integral() << endl;
    
    outfile << "Difference w.r.t. Nominal (%):" << endl;
    outfile << "nominal_p1: " << (histoPmu_p1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
    outfile << "nominal_m1: " << (histoPmu_m1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
    outfile << "--------------------------------------" << endl << endl;
    
    
    if(makeLatex) {
      if ( (variable == 2) || (variable == 3) ) {
        if (function == 0) latexFile << "Nominal" << " & " << histoPmu->Integral() << " & 0" << "\\\\" << endl;
        latexFile << "\\midrule" << endl;
        latexFile << "$" << GetLegendName(functionsName->at(function)) << " + 1\\sigma$ & " << histoPmu_p1->Integral() << " & " << (histoPmu_p1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << "\\\\" << endl;
        latexFile << "$" << GetLegendName(functionsName->at(function)) << " - 1\\sigma$ & " << histoPmu_m1->Integral() << " & " << (histoPmu_m1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << "\\\\" << endl;
      }
      
      if (variable == 0 || variable == 1) {  // save the value of the efficiency to the LaTeX file
        
        double eff_nom = sel_evts_nominal/all_evts_nominal            * 100.;
        double eff_p1  = sel_evts_p1[function]/all_evts_p1[function]  * 100.;
        double eff_m1  = sel_evts_m1[function]/all_evts_m1[function]  * 100.;

        if (function == 0) latexFile << "Nominal" << " & " << std::setprecision(4) << eff_nom << " & 0" << "\\\\" << endl;
        latexFile << "\\midrule" << endl;
        
        latexFile << "$" << GetLegendName(functionsName->at(function)) << " + 1\\sigma$ & "
        << std::setprecision(4) << eff_p1 << " & "
        << std::setprecision(4) << (eff_p1-eff_nom)/eff_nom*100. << "\\\\" << endl;
        
        latexFile << "$" << GetLegendName(functionsName->at(function)) << " - 1\\sigma$ & "
        << std::setprecision(4) << eff_m1 << " & "
        << std::setprecision(4) << (eff_m1-eff_nom)/eff_nom*100. << "\\\\" << endl;
      }
      
    }
    
    
    
    // Define the Canvas
    TCanvas *c = new TCanvas("c", "canvas", 800, 800);
    
    
    
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    histoPmu_p1->SetMinimum(0.0001); // Otherwise 0 label overlaps
    if (variable == 0 || variable == 1) histoPmu_p1->SetMaximum(0.9);
    histoPmu_p1->SetStats(0);          // No statistics on upper plot
    histoPmu_m1->SetStats(0);          // No statistics on upper plot
    histoPmu->SetStats(0);          // No statistics on upper plot
    histoPmu_p1->Draw();               // Draw h1
    histoPmu->Draw("same");         // Draw h2 on top of h1
    histoPmu_m1->Draw("same");
    
    //uBooNESimulation();
    
    if (normalised) {
      // TLatex
      double x;
      x = 0.87;
      if (variable == 3) x = 0.46;
      double y = 0.52;
      double size = 28;
      int color = 1;
      int font = 43;
      int align = 32;
      TLatex *latex = new TLatex( x, y, "Area Normalised" );
      latex->SetNDC();
      latex->SetTextSize(size);
      latex->SetTextColor(color);
      latex->SetTextFont(font);
      latex->SetTextAlign(align);
      
      latex->Draw();
      
    }
    
    // Do not draw the Y axis label on the upper plot and redraw a small
    // axis instead, in order to avoid the first label (0) to be clipped.
    histoPmu->GetYaxis()->SetLabelSize(0.);
    TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axis->SetLabelSize(15);
    axis->Draw();
    
    // Legend for the upper plot
    TLegend* leg;
    if (variable == 0 || variable == 2) leg = new TLegend(0.65,0.6,.85,0.87);
    if (variable == 1 || variable == 3) leg = new TLegend(0.216792,0.5723502,0.4172932,0.843318,NULL,"brNDC");
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    //leg->SetHeader("");
    //leg->SetTextFont(42);
    //leg->AddEntry(histoPmu, "BBA2005 (Nominal)");
    //leg->AddEntry(histoPmu_p1, "Dipole");
    leg->AddEntry(histoPmu, "Nominal");
    leg->AddEntry(histoPmu_p1, LegName + " + 1#sigma");
    leg->AddEntry(histoPmu_m1, LegName + " - 1#sigma");
    leg->Draw();
    
    // lower plot will be in pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.5);
    pad2->SetGridx(); // vertical grid
    //pad2->SetGridy(); // orizontal grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    
    
    // Define the first ratio plot
    TH1F *ratio_p1 = (TH1F*)histoPmu->Clone("ratio_p1");
    //ratio_p1->SetMinimum(0.92);  // Define Y ..
    //ratio_p1->SetMaximum(1.08); // .. range
    //ratio_p1->Sumw2();
    ratio_p1->SetStats(0);      // No statistics on lower plot
    ratio_p1->Divide(histoPmu_p1);
    ratio_p1->SetLineWidth(2);
    ratio_p1->SetLineColor(kRed+1);
    //ratio_p1->Draw("hist");       // Draw the ratio plot
    
    // Define the second ratio plot
    TH1F *ratio_m1 = (TH1F*)histoPmu->Clone("ratio_m1");
    //ratio_m1->SetMinimum(0.9);  // Define Y ..
    //ratio_m1->SetMaximum(1.1); // .. range
    //ratio_m1->Sumw2();
    ratio_m1->SetStats(0);      // No statistics on lower plot
    ratio_m1->Divide(histoPmu_m1);
    ratio_m1->SetLineWidth(2);
    ratio_m1->SetLineColor(kGreen+2);
    //ratio_m1->Draw("hist same");       // Draw the ratio plot
    
    
    
    // Try to set the Y range for the ratio plots
    double max = ratio_p1->GetMaximum();
    double min = ratio_p1->GetMinimum();
    if (ratio_m1->GetMaximum() > max) max = ratio_m1->GetMaximum();
    if (ratio_m1->GetMinimum() < min) min = ratio_m1->GetMinimum();
    //ratio_p1->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
    //ratio_p1->SetMinimum(max);  // Define Y ..
    //ratio_p1->SetMaximum(min); // .. range
    //ratio_m1->SetMinimum(max);  // Define Y ..
    //ratio_m1->SetMaximum(min); // .. range
    cout << functionsName->at(function) <<endl;
    cout << "max: " << max << endl;
    cout << "min: " << min << endl;
    
    // Draw the ratio plot
    //ratio_p1->Draw("hist");
    //ratio_m1->Draw("hist same");
    
    THStack *hs = new THStack("hs","");
    hs->Add(ratio_p1);
    hs->Add(ratio_m1);
    hs->SetMaximum(hs->GetMaximum("nostack")+0.01/**hs->GetMaximum("nostack")*/);
    hs->SetMinimum(hs->GetMinimum("nostack")-0.01/**hs->GetMinimum("nostack")*/);
    cout << "hs->GetMinimum(): " << hs->GetMinimum("nostack") << endl;
    hs->Draw("NOSTACK");
    
    
    //ratio_p1->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
    
    
    
    
    //**********************
    //
    // Settings
    //
    //**********************
    
    // h1 settings
    histoPmu->SetLineColor(kBlack);
    histoPmu->SetLineWidth(2);
    
    // Y axis h1 plot settings
    histoPmu_p1->GetYaxis()->CenterTitle();
    histoPmu_p1->GetYaxis()->SetTitleSize(25);
    histoPmu_p1->GetYaxis()->SetTitleFont(43);
    histoPmu_p1->GetYaxis()->SetTitleOffset(1.55);
    
    // h2 settings
    histoPmu_p1->SetLineColor(kRed+1);
    histoPmu_p1->SetLineWidth(2);
    
    // h3 settings
    histoPmu_m1->SetLineColor(kGreen+2);
    histoPmu_m1->SetLineWidth(2);
    
    // Ratio plot (ratio_p1) settings
    ratio_p1->SetTitle(""); // Remove the ratio title
    
    
    hs->GetYaxis()->SetTitle("Ratio");
    hs->GetYaxis()->CenterTitle();
    hs->GetYaxis()->SetNdivisions(505);
    hs->GetYaxis()->SetTitleSize(25);
    hs->GetYaxis()->SetTitleFont(43);
    hs->GetYaxis()->SetTitleOffset(1.0);
    hs->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hs->GetYaxis()->SetLabelSize(15);
    
    
    if(variable == 0 || variable == 2) hs->GetXaxis()->SetTitle("p_{#mu} [GeV]");
    if(variable == 1 || variable == 3) hs->GetXaxis()->SetTitle("cos#theta_{#mu}");
    hs->GetXaxis()->CenterTitle();
    hs->GetXaxis()->SetTitleSize(25);
    hs->GetXaxis()->SetTitleFont(43);
    hs->GetXaxis()->SetTitleOffset(3.5);
    hs->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    hs->GetXaxis()->SetLabelSize(20);
    
    // Draw linea at 1 in ratio plot
    TLine *line;
    if (variable == 0 || variable == 2) line = new TLine(0,1,2,1);
    if (variable == 1 || variable == 3) line = new TLine(-1,1,1,1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(9); // dashed
    line->Draw();
    
    if (variable == 0 || variable == 1) {
      c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".C");
      c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".pdf");
    }
    if (variable == 2 || variable == 3) {
      c->Print("./EvtWgtEventPlots/" + SaveName + ".C");
      c->Print("./EvtWgtEventPlots/" + SaveName + ".pdf");
    }
    
  } // end loop functions
  
  if(makeLatex) {
    latexFile << "\\bottomrule" << endl;
    latexFile << "\\end{tabular}" << endl;
    latexFile << "\\end{table}" << endl;
  }
  
}









//___________________________________________________
TString EvtWgtAnalysis::GetLegendName(string fName) {
  
  TString legName = "null";
  if (fName.find("qema") != std::string::npos) legName = "M_{A}^{CCQE}";
  if (fName.find("ncelEta") != std::string::npos) legName = "#eta^{NCEL}";
  if (fName.find("ncelAxial") != std::string::npos) legName = "M_{A}^{NCEL}";
  if (fName.find("qevec") != std::string::npos) legName = "VecFF";
  if (fName.find("ccresAxial") != std::string::npos) legName = "M_{A}^{CCRES}";
  if (fName.find("ccresVector") != std::string::npos) legName = "M_{V}^{CCRES}";
  if (fName.find("ncresAxial") != std::string::npos) legName = "M_{A}^{NCRES}";
  if (fName.find("ncresVector") != std::string::npos) legName = "M_{V}^{NCRES}";
  if (fName.find("cohMA") != std::string::npos) legName = "M_{A}^{COH#pi}";
  if (fName.find("cohR0") != std::string::npos) legName = "R_{0}^{COH#pi}";
  if (fName.find("NonResRvp1pi") != std::string::npos) legName = "NonResRvp1pi";
  if (fName.find("NonResRvbarp1pi") != std::string::npos) legName = "NonResRvbarp1pi";
  if (fName.find("NonResRvp2pi") != std::string::npos) legName = "NonResRvp2pi";
  if (fName.find("NonResRvbarp2pi") != std::string::npos) legName = "NonResRvbarp2pi";
  if (fName.find("ResDecayGamma") != std::string::npos) legName = "BR(#gamma)";
  if (fName.find("ResDecayEta") != std::string::npos) legName = "BR(#eta)";
  if (fName.find("ResDecayTheta") != std::string::npos) legName = "Ang. distr.";
  if (fName.find("DISAth") != std::string::npos) legName = "A_{HT}^{BY}";
  if (fName.find("DISBth") != std::string::npos) legName = "B_{HT}^{BY}";
  if (fName.find("DISCv1u") != std::string::npos) legName = "C_{V1u}^{BY}";
  if (fName.find("DISCv2u") != std::string::npos) legName = "C_{V2u}^{BY}";
  if (fName.find("AGKYxF") != std::string::npos) legName = "x_{F}";
  if (fName.find("AGKYpT") != std::string::npos) legName = "p_{T}";
  if (fName.find("FormZone") != std::string::npos) legName = "FZ";
  if (fName.find("FermiGasModelKf") != std::string::npos) legName = "CCQE-PauliSup";
  if (fName.find("FermiGasModelSf") != std::string::npos) legName = "Fermi Gas to SF";
  if (fName.find("IntraNukeNmfp") != std::string::npos) legName = "x_{mfp}^{N}";
  if (fName.find("IntraNukeNcex") != std::string::npos) legName = "x_{cex}^{N}";
  if (fName.find("ntraNukeNel") != std::string::npos) legName = "x_{el}^{N}";
  if (fName.find("ntraNukeNinel") != std::string::npos) legName = "x_{inel}^{N}";
  if (fName.find("ntraNukeNabs") != std::string::npos) legName = "x_{abs}^{N}";
  if (fName.find("ntraNukeNpi") != std::string::npos) legName = "x_{pi}^{N}";
  if (fName.find("IntraNukePImfp") != std::string::npos) legName = "x_{mfp}^{PI}";
  if (fName.find("IntraNukePIcex") != std::string::npos) legName = "x_{cex}^{PI}";
  if (fName.find("ntraNukePIel") != std::string::npos) legName = "x_{el}^{PI}";
  if (fName.find("ntraNukePIinel") != std::string::npos) legName = "x_{inel}^{PI}";
  if (fName.find("ntraNukePIabs") != std::string::npos) legName = "x_{abs}^{PI}";
  if (fName.find("ntraNukePIpi") != std::string::npos) legName = "x_{pi}^{PI}";
  
  return legName;
}








