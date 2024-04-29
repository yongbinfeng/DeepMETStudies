# DeepMET studies
includes code to calculate and compare the resolution difference between different METs.

The code needs python3 and [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html). One can source the provided environment on `/cvmfs` to get started:
```bash
source setenv.sh
```

## Produce the slimmed ntuples
The slimmed ntuples are produced using the `PrepareSlimmedNtuples.py` script. The script reads the input files, applies the selection and produces the output files. This way we can reduce the size of the ntuples and subsequent IO processings in order to speed up the analysis.

First go to `RecoilResol` directory and compile the helper functions:
```bash
root [0] .L Functions.cc+
```
Then run the script:
```bash
python PrepareSlimmedNtuples.py
```

## Recoil Resolutions

Some scripts used to compare the responses and (un)scaled resolutions of different recoil estimators.
The core is the class `RecoilAnalyzer` in `RecoilResol/utils/RecoilAnalyzer.py`, which prepares the RDataFrame, 
loads the functions, prepare the responses and resolutions. The other scripts are keeping constructing this class and 
compare the performances (responses and resolutions).

To make the plots for UL data and MC, firtly compile some helper functions in `RecoilResol/Functions.cc`. This only need to be done once.
```bash
root [0] .L Functions.cc+
```

Then run the scripts:
```bash
python recoil_resol_dataMC.py
```
to make the resolution plots for data and MC on $Z\to\mu\mu$ events. The script will read the slimmed ntuples and make the plots of the responses and resolutions of different METs, with and without response corrections, before and after the calibrations.


## Recoil Corrections
(To be updated.)
Code used to study the corrections and uncertainties of DeepMET. In principle it should also work on any type of MET with minor modification.

### Run the correction code.
Then prepare the uparal(u1) and uperp (u2) distributions using `MakeU1U2.py`. One could change the flags `isData`, `isAMCNLO`,
`isMADGRAPH`, `isTTbar`, `isTTbar` to prepare the distributions for different samples.

`fitRecoilZmm.py` will read the inputs and run the multi-Gaussian fits on data or MCs. Fit results will be saved and plots will be created.

`PdfToCdf_njets_pt.py` will take the pdfs from the fit, convert it to cdfs and save all of them into the ROOT files. It also allows to do the Gaussian Kernel Smoothing from the original histograms instead of the fitting, which can be used as one source of systematic uncertainties.

`MakePlots_Z.py` is the final plotting script: plot the data-MC comparisons pre/post corrections, and do the systematic variations and plot the uncs.

`MakePlots_W.py` functions the same way as `MakePlots_Z.py`, except it runs on samples in W signal region. Thus it is very slow, since we have more than 80 million W->munu events after the lepton id and kinematic selection in 2016...


