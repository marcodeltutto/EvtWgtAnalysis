# EvtWgtAnalysis
Code to study GENIE model uncertainties at MicroBooNE.

##Compile

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

It works with `root v5_34_32`.

```
source SetupEvtWgtAnalysis.sh

root -l
> gSystem->Load("AnaTree/AnaBNB_C.so");
> gSystem->Load("EvtWgtAnalysis_cxx.so");
> EvtWgtAnalysis t("generic_*_to_anatree.root");
> f.MakeBackgroundPlots();
```

To make the Number of Events or efficiency plots:
```
> f.MakePlots(false, option);
```

where `option` is an `integer`:
- `option=0`: efficiency - Pmu
- `option=1`: efficiency - CosThetaMu
- `option=2`: events - Pmu
- `option=3`: events - CosThetaMu

Change `false` to `true` to have area normalized plots.


To make the plots that show the events after the selection and the background decomposition:
```
> MakeBackgroundPlots(option)
```

where `option` is an `integer`:
- `option=0`: as a function of Pmu
- `option=1`: as a function of CosThetaMu

To make the plots that show the cross-section percental difference:
```
> MakeXsecDiffPlots()
```


##Generate new AnaBNB_C.so 

Change file path in makeClass.C first, then:

```
cd AnaTree
root -l makeClass.C
root -l
> .L AnaBNB.C+
```


