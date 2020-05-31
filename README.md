# DeepMET studies
includes code to calculate and compare the resolution difference between different METs.

The code needs [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html), which becomes an official part
of ROOT since ROOT 6.14. Tested and confirmed at least CMSSW_10_6_0 works (with ROOT 6.14.09)

## Recoil Resolutions

Code used to compare the responses and (un)scaled resolutions of different recoil estimators

## Recoil Corrections

Code used to study the corrections and uncertainties of DeepMET. In principle it should also work on any type of MET with minor modification.

### Run the correction code.

First compile the macros in `Functions.cc`

Then prepare the uparal(u1) and uperp (u2) distributions using `MakeU1U2.py`. One could change the flags `isData`, `isAMCNLO`,
`isMADGRAPH`, `isTTbar`, `isTTbar` to prepare the distributions for different samples.

`fitRecoilZmm.py` will read the inputs and run the multi-Gaussian fits on data or MCs. Fit results will be saved and plots will be created.

`PdfToCdf_njets_pt.py` will take the pdfs from the fit, convert it to cdfs and save all of them into the ROOT files. It also allows to do the Gaussian Kernel Smoothing from the original histograms instead of the fitting, which can be used as one source of systematic uncertainties.

`MakePlots_Z.py` is the final plotting script: plot the data-MC comparisons pre/post corrections, and do the systematic variations and plot the uncs.

`MakePlots_W.py` functions the same way as `MakePlots_Z.py`, except it runs on samples in W signal region. Thus it is very slow, since we have more than 80 million W->munu events after the lepton id and kinematic selection in 2016...


