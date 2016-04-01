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
        for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
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
            if(ccnc_truth == 0 && nuPDG_truth == 14) {
              // If I got till this point, it means I have selected the event
              // Fill the histograms at this point
              xsec_mom_data -> Fill(lep_mom_truth);
              xsec_mom_reco_data -> Fill(anaBNBtree->trkmomrange_pandoraNuKHit[index_longest_track]);
              xsec_theta_data -> Fill(lep_dcosz_truth);
              if(vertexatstart[index_longest_track]) xsec_theta_reco_data -> Fill(cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]));
              if(vertexatend[index_longest_track]) xsec_theta_reco_data -> Fill(-cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]));
              // Loop over the re-weighting functions and fill +-1 sigma histograms
              for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
                xsec_mom_data_p1[function] -> Fill(lep_mom_truth, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_mom_data_m1[function] -> Fill(lep_mom_truth, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                xsec_mom_reco_data_p1[function] -> Fill(anaBNBtree->trkmomrange_pandoraNuKHit[index_longest_track], (anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_mom_reco_data_m1[function] -> Fill(anaBNBtree->trkmomrange_pandoraNuKHit[index_longest_track], (anaBNBtree->evtwgt_weight->at(function)).at(1));
                xsec_theta_data_p1[function] -> Fill(lep_dcosz_truth, (anaBNBtree->evtwgt_weight->at(function)).at(0));
                xsec_theta_data_m1[function] -> Fill(lep_dcosz_truth, (anaBNBtree->evtwgt_weight->at(function)).at(1));
                if(vertexatstart[index_longest_track])
                  xsec_theta_reco_data_p1[function]
                  -> Fill(cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]), (anaBNBtree->evtwgt_weight->at(function)).at(0));
                if(vertexatstart[index_longest_track])
                  xsec_theta_reco_data_m1[function]
                  -> Fill(cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]), (anaBNBtree->evtwgt_weight->at(function)).at(1));
                if(vertexatend[index_longest_track])
                  xsec_theta_reco_data_p1[function]
                  -> Fill(-cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]), (anaBNBtree->evtwgt_weight->at(function)).at(0));
                if(vertexatend[index_longest_track])
                  xsec_theta_reco_data_m1[function]
                  -> Fill(-cos(anaBNBtree->trktheta_pandoraNuKHit[index_longest_track]), (anaBNBtree->evtwgt_weight->at(function)).at(1));
              } // end loop over re-weighting functions
            }
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
  // Save to file
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
  std::cout << "Histograms saved to eventWeightEfficiency.root" << std::endl;
  
}


//_______________________________________
bool EvtWgtAnalysis::inFV(double x, double y, double z) {
  
  double FVx = 256.35;
  double FVy = 233;
  double FVz = 1036.8;
  double border = 5.;
  
  if(x < (FVx - border) && x > border && y < (FVy/2. - border) && y > (-FVy/2. + border) && z < (FVz - border) && z > border) return true;
  else return false;
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
      if (anaBNBtree->mode_truth[0] == 3) Ncoh++;
      // Loop over GENIE particles
      for (int GENIEparticle = 0; GENIEparticle < anaBNBtree->genie_no_primaries; GENIEparticle++) {
        // Look for muon
        if (anaBNBtree->genie_primaries_pdg[GENIEparticle] == 14-1 && anaBNBtree->genie_status_code[GENIEparticle]==1) {
          // Loop over the reweighting functions
          for (int function = 0; function < anaBNBtree->evtwgt_nfunc; function++) {
            if (anaBNBtree->evtwgt_funcname->at(function) == "genie_coh_Genie") temp->Fill((anaBNBtree->evtwgt_weight->at(function)).at(0));
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
    } // end of loop over entries
  } // end cc/nc and neutrino type selection
  std::cout << "Loop completed. Saving to eventWeightHistograms.root now." << std::endl;
  
  std::cout << "Number of Coherent production interactions: " << Ncoh << std::endl;
  
  // Save to file
  TFile f = TFile("eventWeightHistograms.root", "RECREATE");
  temp->Write();
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
  
  // Instantiate the histograms
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


//__________________________________________________
void EvtWgtAnalysis::MakePlotsPmu(bool normalised) {
 
  bool makeLatex = true;

  if (histoPmuVec.size() == 0) {
    cout << "Calling MakeHistograms() first."  << endl;
    MakeHistograms();
  }
  
  // Make a directory to store the plots
  system("mkdir ./EvtWgtPlots");
  
  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);
  
  // Opening a text file to write the integrals
  ofstream outfile;
  outfile.open("./EvtWgtPlots/Integrals.txt");
  
  // Open LaTeX file
  ofstream latexFile;
  if(makeLatex) {
    latexFile.open("./EvtWgtPlots/evtwgt.tex");
    latexFile << "\\begin{table}[]" << endl;
    latexFile << "\\caption{}" << endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << endl;
    latexFile << "\\label{tab:}" << endl;
    latexFile << "\\centering" << endl;
    latexFile << "\\begin{tabular}{ccc}" << endl;
    latexFile << "\\toprule" << endl;
    latexFile << "  &  Integral  &  Difference (\\%) \\\\" << endl;
    latexFile << "\\midrule" << endl;
  }

  // Looping over all the functions
  for (unsigned int function = 0; function < histoPmuVec.size(); function++) {
    
    TH1F *histoPmu    = (TH1F*)histoPmuVec.at(function)->Clone("histoPmu");
    TH1F *histoPmu_p1 = (TH1F*)histoPmuVec_p1.at(function)->Clone("histoPmu_p1");
    TH1F *histoPmu_m1 = (TH1F*)histoPmuVec_m1.at(function)->Clone("histoPmu_m1");
    
    TString SaveName = "Pmu"+functionsName->at(function);
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
    outfile << "Integral histoPmu:    " << histoPmu->Integral() << endl;
    outfile << "Integral histoPmu_p1: " << histoPmu_p1->Integral() << endl;
    outfile << "Integral histoPmu_m1: " << histoPmu_m1->Integral() << endl;
    
    outfile << "Difference w.r.t. histoPmu (%):" << endl;
    outfile << "histoPmu_p1: " << (histoPmu_p1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
    outfile << "histoPmu_p1: " << (histoPmu_m1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << endl;
    outfile << "--------------------------------------" << endl << endl;
    

    if(makeLatex) {
      if (function == 0) latexFile << "$" << "Nominal" << "$ & " << histoPmu->Integral() << " & 0" << "\\\\" << endl;
      latexFile << "\\midrule" << endl;
      latexFile << "$" << GetLegendName(functionsName->at(function)) << " + 1\\sigma$ & " << histoPmu_p1->Integral() << " & " << (histoPmu_p1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << "\\\\" << endl;
      latexFile << "$" << GetLegendName(functionsName->at(function)) << " - 1\\sigma$ & " << histoPmu_m1->Integral() << " & " << (histoPmu_m1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << "\\\\" << endl;

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
    histoPmu_p1->SetStats(0);          // No statistics on upper plot
    histoPmu_m1->SetStats(0);          // No statistics on upper plot
    histoPmu->SetStats(0);          // No statistics on upper plot
    histoPmu_p1->Draw();               // Draw h1
    histoPmu->Draw("same");         // Draw h2 on top of h1
    histoPmu_m1->Draw("same");
    
    //uBooNESimulation();
    
    if (normalised) {
      // TLatex
      double x = 0.87;
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
    TLegend* leg = new TLegend(0.65,0.6,.85,0.87);
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
    ratio_p1->SetMinimum(0.75);  // Define Y ..
    ratio_p1->SetMaximum(1.25); // .. range
    ratio_p1->Sumw2();
    ratio_p1->SetStats(0);      // No statistics on lower plot
    ratio_p1->Divide(histoPmu_p1);
    ratio_p1->SetLineWidth(2);
    ratio_p1->SetLineColor(kRed+1);
    ratio_p1->Draw("hist");       // Draw the ratio plot
    
    // Define the second ratio plot
    TH1F *ratio_m1 = (TH1F*)histoPmu->Clone("ratio_m1");
    ratio_m1->SetMinimum(0.8);  // Define Y ..
    ratio_m1->SetMaximum(1.2); // .. range
    ratio_m1->Sumw2();
    ratio_m1->SetStats(0);      // No statistics on lower plot
    ratio_m1->Divide(histoPmu_m1);
    ratio_m1->SetLineWidth(2);
    ratio_m1->SetLineColor(kGreen+2);
    ratio_m1->Draw("hist same");       // Draw the ratio plot
    
    
    
    
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
    
    // Y axis ratio plot settings
    ratio_p1->GetYaxis()->SetTitle("Ratio");
    ratio_p1->GetYaxis()->CenterTitle();
    ratio_p1->GetYaxis()->SetNdivisions(505);
    ratio_p1->GetYaxis()->SetTitleSize(25);
    ratio_p1->GetYaxis()->SetTitleFont(43);
    ratio_p1->GetYaxis()->SetTitleOffset(1.0);
    ratio_p1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_p1->GetYaxis()->SetLabelSize(15);
    
    // X axis ratio plot settings
    ratio_p1->GetXaxis()->CenterTitle();
    ratio_p1->GetXaxis()->SetTitleSize(25);
    ratio_p1->GetXaxis()->SetTitleFont(43);
    ratio_p1->GetXaxis()->SetTitleOffset(3.5);
    ratio_p1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_p1->GetXaxis()->SetLabelSize(20);
    
    // Draw linea at 1 in ratio plot
    TLine *line = new TLine(0,1,2,1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(9); // dashed
    line->Draw();
    
    
    c->Print("./EvtWgtPlots/" + SaveName + ".C");
    c->Print("./EvtWgtPlots/" + SaveName + ".pdf");
    
    
    /*
     //**********************
     //
     // Area Normalised plot
     //
     //**********************
     
     TH1F *histoPmu_copy = (TH1F*)histoPmu->Clone("histoPmu_copy");
     TH1F *histoPmu_p1_copy = (TH1F*)histoPmu_p1->Clone("histoPmu_p1_copy");
     TH1F *histoPmu_m1_copy = (TH1F*)histoPmu_m1->Clone("histoPmu_m1_copy");
     
     // Area Norm
     histoPmu_copy->Scale(1./histoPmu_copy->Integral());
     histoPmu_p1_copy->Scale(1./histoPmu_p1_copy->Integral());
     histoPmu_m1_copy->Scale(1./histoPmu_m1_copy->Integral());
     
     // Draw them in a new canvas
     TCanvas *c1 = new TCanvas("c1", "canvas", 800, 800);
     histoPmu_p1_copy->Draw();
     histoPmu_copy->Draw("same");
     histoPmu_m1_copy->Draw("same");
     
     //Settings
     histoPmu_p1_copy->GetXaxis()->CenterTitle();
     histoPmu_p1_copy->GetYaxis()->CenterTitle();
     //histoPmu_p1_copy->GetYaxis()->SetTitleOffset(1.0);
     
     leg->Draw();
     
     // TLatex
     double x = 0.87;
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
     */
    
    
  } // end loop functions

  if(makeLatex) {
    latexFile << "\\bottomrule" << endl;
    latexFile << "\\end{tabular}" << endl;
    latexFile << "\\end{table}" << endl;
  }

}










//__________________________________________________________________________
void EvtWgtAnalysis::MakeEfficiencyPlotsPmu(bool normalised, int variable) {
  
  bool makeLatex = true;
  
  if (histoPmuVec.size() == 0) {
    cout << "Calling MakeHistograms() first."  << endl;
    CalcEfficiency();
  }
  
  // Make a directory to store the plots
  system("rm -rf ./EvtWgtEfficiencyPlots/Cos* ");
  system("mkdir ./EvtWgtEfficiencyPlots");
  
  // Avoid root to dislay the canvases
  gROOT->SetBatch(kTRUE);
  
  // Opening a text file to write the integrals
  ofstream outfile;
  outfile.open("./EvtWgtEfficiencyPlots/Integrals.txt");
  
  // Open LaTeX file to write the table
  ofstream latexFile;
  if(makeLatex) {
    latexFile.open("./EvtWgtEfficiencyPlots/evtwgtEfficiency.tex");
    latexFile << "\\begin{table}[]" << endl;
    latexFile << "\\caption{}" << endl;
    latexFile << "\\captionsetup{format=hang,labelfont={sf,bf}}" << endl;
    latexFile << "\\label{tab:}" << endl;
    latexFile << "\\centering" << endl;
    latexFile << "\\begin{tabular}{ccc}" << endl;
    latexFile << "\\toprule" << endl;
    latexFile << "  &  Integral  &  Difference (\\%) \\\\" << endl;
    latexFile << "\\midrule" << endl;
  }
  
  // Looping over all the functions
  for (unsigned int function = 0; function < xsec_mom_eff_p1.size(); function++) {
    
    TH1F *histoPmu;
    TH1F *histoPmu_p1;
    TH1F *histoPmu_m1;
    
    if (variable == 0) {  // Pmu
      histoPmu    = (TH1F*)xsec_mom_eff->Clone("histoPmu");
      histoPmu_p1 = (TH1F*)xsec_mom_eff_p1.at(function)->Clone("histoPmu_p1");
      histoPmu_m1 = (TH1F*)xsec_mom_eff_m1.at(function)->Clone("histoPmu_m1");
    }
    else if (variable == 1) { // CosThetaMu
      histoPmu    = (TH1F*)xsec_theta_eff->Clone("histoPmu");
      histoPmu_p1 = (TH1F*)xsec_theta_eff_p1.at(function)->Clone("histoPmu_p1");
      histoPmu_m1 = (TH1F*)xsec_theta_eff_m1.at(function)->Clone("histoPmu_m1");
    } else {
      cout << "Invalid option. Exit." << endl;
      exit(0);
    }
    TString SaveName;
    if(variable == 0) SaveName = "Pmu"+functionsName->at(function);
    if(variable == 1) SaveName = "CosThetamu"+functionsName->at(function);
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
      if (function == 0) latexFile << "$" << "Nominal" << "$ & " << histoPmu->Integral() << " & 0" << "\\\\" << endl;
      latexFile << "\\midrule" << endl;
      latexFile << "$" << GetLegendName(functionsName->at(function)) << " + 1\\sigma$ & " << histoPmu_p1->Integral() << " & " << (histoPmu_p1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << "\\\\" << endl;
      latexFile << "$" << GetLegendName(functionsName->at(function)) << " - 1\\sigma$ & " << histoPmu_m1->Integral() << " & " << (histoPmu_m1->Integral()-histoPmu->Integral())/(histoPmu->Integral())*100. << "\\\\" << endl;
      
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
    histoPmu_p1->SetMaximum(0.5);
    histoPmu_p1->SetStats(0);          // No statistics on upper plot
    histoPmu_m1->SetStats(0);          // No statistics on upper plot
    histoPmu->SetStats(0);          // No statistics on upper plot
    histoPmu_p1->Draw();               // Draw h1
    histoPmu->Draw("same");         // Draw h2 on top of h1
    histoPmu_m1->Draw("same");
    
    //uBooNESimulation();
    
    if (normalised) {
      // TLatex
      double x = 0.87;
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
    TLegend* leg = new TLegend(0.65,0.6,.85,0.87);
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
    ratio_p1->SetMinimum(0.85);  // Define Y ..
    ratio_p1->SetMaximum(1.15); // .. range
    ratio_p1->Sumw2();
    ratio_p1->SetStats(0);      // No statistics on lower plot
    ratio_p1->Divide(histoPmu_p1);
    ratio_p1->SetLineWidth(2);
    ratio_p1->SetLineColor(kRed+1);
    ratio_p1->Draw("hist");       // Draw the ratio plot
    
    // Define the second ratio plot
    TH1F *ratio_m1 = (TH1F*)histoPmu->Clone("ratio_m1");
    ratio_m1->SetMinimum(0.9);  // Define Y ..
    ratio_m1->SetMaximum(1.1); // .. range
    ratio_m1->Sumw2();
    ratio_m1->SetStats(0);      // No statistics on lower plot
    ratio_m1->Divide(histoPmu_m1);
    ratio_m1->SetLineWidth(2);
    ratio_m1->SetLineColor(kGreen+2);
    ratio_m1->Draw("hist same");       // Draw the ratio plot
    
    
    
    
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
    
    // Y axis ratio plot settings
    ratio_p1->GetYaxis()->SetTitle("Ratio");
    ratio_p1->GetYaxis()->CenterTitle();
    ratio_p1->GetYaxis()->SetNdivisions(505);
    ratio_p1->GetYaxis()->SetTitleSize(25);
    ratio_p1->GetYaxis()->SetTitleFont(43);
    ratio_p1->GetYaxis()->SetTitleOffset(1.0);
    ratio_p1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_p1->GetYaxis()->SetLabelSize(15);
    
    // X axis ratio plot settings
    if(variable == 0) ratio_p1->GetXaxis()->SetTitle("p_{#mu} [GeV]");
    if(variable == 1) ratio_p1->GetXaxis()->SetTitle("cos#theta_{#mu}");
    ratio_p1->GetXaxis()->CenterTitle();
    ratio_p1->GetXaxis()->SetTitleSize(25);
    ratio_p1->GetXaxis()->SetTitleFont(43);
    ratio_p1->GetXaxis()->SetTitleOffset(3.5);
    ratio_p1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    ratio_p1->GetXaxis()->SetLabelSize(20);
    
    // Draw linea at 1 in ratio plot
    TLine *line;
    if (variable == 0) line = new TLine(0,1,2,1);
    if (variable == 1) line = new TLine(0,-1,1,1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(9); // dashed
    line->Draw();
    
    
    c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".C");
    c->Print("./EvtWgtEfficiencyPlots/" + SaveName + ".pdf");
    
    
    /*
     //**********************
     //
     // Area Normalised plot
     //
     //**********************
     
     TH1F *histoPmu_copy = (TH1F*)histoPmu->Clone("histoPmu_copy");
     TH1F *histoPmu_p1_copy = (TH1F*)histoPmu_p1->Clone("histoPmu_p1_copy");
     TH1F *histoPmu_m1_copy = (TH1F*)histoPmu_m1->Clone("histoPmu_m1_copy");
     
     // Area Norm
     histoPmu_copy->Scale(1./histoPmu_copy->Integral());
     histoPmu_p1_copy->Scale(1./histoPmu_p1_copy->Integral());
     histoPmu_m1_copy->Scale(1./histoPmu_m1_copy->Integral());
     
     // Draw them in a new canvas
     TCanvas *c1 = new TCanvas("c1", "canvas", 800, 800);
     histoPmu_p1_copy->Draw();
     histoPmu_copy->Draw("same");
     histoPmu_m1_copy->Draw("same");
     
     //Settings
     histoPmu_p1_copy->GetXaxis()->CenterTitle();
     histoPmu_p1_copy->GetYaxis()->CenterTitle();
     //histoPmu_p1_copy->GetYaxis()->SetTitleOffset(1.0);
     
     leg->Draw();
     
     // TLatex
     double x = 0.87;
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
     */
    
    
  } // end loop functions
  
  if(makeLatex) {
    latexFile << "\\bottomrule" << endl;
    latexFile << "\\end{tabular}" << endl;
    latexFile << "\\end{table}" << endl;
  }
  
}









//________________________________________________ 
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








