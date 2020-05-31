# Recoil Corrections

Code used to study the corrections and uncertainties of DeepMET. In principle it should work on any type of MET. 

The code needs [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html), which becomes an official part
of ROOT since ROOT 6.14. Tested and confirmed at least CMSSW_10_6_0 works (with ROOT 6.14.09)

First compile the macros in Functions.cc

Then prepare the uparal(u1) and uperp (u2) distributions using MakeU1U2.py. One could change the flags `isData`, `isAMCNLO`,
`isMADGRAPH`, `isTTbar`, `isTTbar` to prepare the distributions for different samples.


