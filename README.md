# EvtWgtAnalysis
Code to study GENIE model uncertainties at MicroBooNE.


Run with:

```
source SetupEvtWgtAnalysis.sh

root -l
> gSystem->Load("AnaTree/AnaBNB_C.so");
> gSystem->Load("EvtWgtAnalysis_cxx.so");
> EvtWgtAnalysis t("generic_*_to_anatree.root");
> f.MakeBackgroundPlots();
```

or use the python script:

```
source SetupEvtWgtAnalysis.sh
python RunEvtWgtAnalysis.py
```

Compile with:

```
root -l
> gSystem->Load("AnaTree/AnaBNB_C.so");
> .L EvtWgtAnalysis.cxx+
```


Generate new AnaBNB_C.so:

```
cd AnaTree
root -l makeClass.C
root -l
.L AnaBNB.C+
```
