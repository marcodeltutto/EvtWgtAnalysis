# EvtWgtAnalysis
Code to study GENIE model uncertainties at MicroBooNE.

#Compile

```
cd AnaTree
root -l
.L AnaBNB.C+
.q
cd ..
root -l
> gSystem->Load("AnaTree/AnaBNB_C.so");
> .L EvtWgtAnalysis.cxx+
```



##Run

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

##Generate new AnaBNB_C.so 

Change file path in makeClass.C first, then:

```
cd AnaTree
root -l makeClass.C
root -l
.L AnaBNB.C+
```


