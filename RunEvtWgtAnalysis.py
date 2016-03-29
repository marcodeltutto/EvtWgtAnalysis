#!/usr/bin/env python
#

import os,sys,string, time
import ROOT
ROOT.gSystem.Load("AnaTree/AnaBNB_C.so")
ROOT.gSystem.Load("EvtWgtAnalysis_cxx.so")
from ROOT import EvtWgtAnalysis
from glob import glob

fname = glob("/pnfs/uboone/scratch/users/mdeltutt/v04_30_03/anatree_bnb_eventWeight_all/5033083_*/standard_reco_hist.root")

f = EvtWgtAnalysis("/pnfs/uboone/scratch/users/mdeltutt/v04_30_03/anatree_bnb_eventWeight_all/reunion/standard_reco_hist_99*.root")
#f = EvtWgtAnalysis("../standard_reco_hist.root")
f.MakeHistograms()


raw_input("Please press enter to exit.")
