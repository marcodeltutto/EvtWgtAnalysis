void makeClass() {

  //TFile *f = new TFile("../../standard_reco_hist.root");
  //TFile *f = new TFile("/pnfs/uboone/scratch/users/mdeltutt/v04_30_03/anatree_bnb_eventWeight_all/reunion/standard_reco_hist_1.root");
  //TFile *f = new TFile("/pnfs/uboone/scratch/users/mdeltutt/v04_30_03/anatree_bnb_eventWeight_all_5/5285883_0/standard_reco_hist.root");
  TFile *f = new TFile("../standard_reco_hist_merged_all_5.root");
  f->cd("analysistree");
  TTree *v = (TTree*)f->Get("anatree");
  anatree->MakeClass("AnaBNB");
}
