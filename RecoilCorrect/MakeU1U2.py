'''
plot and save the u1(uparal) and u2(uperp) distributions 
in the Z SR, which will be used in the recoil corrections in 
next steps. It can be run on data/Z amc@NLO/Z Madgraph/ttbar bkg
/Dibson bkg etc. QCD scale variations can be taken into account.
'''

import ROOT
import numpy as np
import os
from collections import OrderedDict
from SampleManager import Sample
import CMS_lumi
import pickle
from Utils.utils import getPtBins, getJetBins

isData     = True
isAMCNLO   = False
isMADGRAPH = False
isTTbar    = False
isDiboson  = False

assert isData+isAMCNLO+isMADGRAPH+isTTbar+isDiboson==1, "must pick one sample from data, or amc@nlo, or madgraph"

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(18)

def main():
    
    input_data    = "inputs/inputs_Z_UL_Post/input_data.txt"
    input_dy      = "inputs/inputs_Z_UL_Post/input_zjets.txt"
    input_ttbar   = "inputs/inputs_Z_UL_Post/input_ttbar.txt"
    input_WW2L    = "inputs/inputs_Z_UL_Post/input_WWTo2L2Nu.txt"
    input_WZ2L    = "inputs/inputs_Z_UL_Post/input_WZTo2Q2L.txt"
    input_ZZ2L    = "inputs/inputs_Z_UL_Post/input_ZZTo2L2Nu.txt"
    input_ZZ2L2Q  = "inputs/inputs_Z_UL_Post/input_ZZTo2Q2L.txt"
    input_dytau   = "inputs/inputs_Z_UL_Post/input_zjets_tautau.txt"

    if isData:
        samp = Sample(input_data, isMC=False, legend="Data", name="Data", prepareVars=False, select=False)
        samp.donormalization = False
        samp.varyQCDScale = False
        samp.Define("norm", "1")
        samples = [samp]
    elif isAMCNLO:
        samp = Sample(input_dy,    xsec = 0,  color=5,  reweightzpt = False, legend="DY", name="DY", prepareVars=False, select=False)
        samp.donormalization = False
        samp.varyQCDScale = False
        samples = [samp]
    elif isMADGRAPH:
        # outdated for now
        inputfile = "inputs/inputs_Z/input_zjets_madgraph.txt"
        samp = Sample(inputfile, isMC=True,  name = "ZJets_MG")
        samp.donormalization = False
        samp.varyQCDScale = False
        samples = [samp]
    elif isTTbar:
        samp = Sample(input_ttbar, xsec = 0,  color=46, reweightzpt = False, legend="t#bar{t}", name="ttbar", prepareVars=False, select=False)
        samp.donormalization = True
        samp.varyQCDScale = False
        samples = [samp]
    elif isDiboson:
        WW2LSamp  = Sample(input_WW2L,  xsec = 0,  color=38, reweightzpt = False, legend="WW2L",     name="WW2L",  prepareVars=False, select=False)
        WZ2LSamp  = Sample(input_WZ2L,  xsec = 0,  color=39, reweightzpt = False, legend="WZ2L",     name="WZ2L",  prepareVars=False, select=False)
        ZZ2LSamp  = Sample(input_ZZ2L,  xsec = 0,  color=37, reweightzpt = False, legend="ZZ2L",     name="ZZ2L",  prepareVars=False, select=False)
        ZZ2L2QSamp  = Sample(input_ZZ2L2Q,  xsec = 0, color=36, reweightzpt = False, legend="ZZ2L2Q", name="ZZ2L2Q",  prepareVars=False, select=False) 
        DYTauSamp   = Sample(input_dytau, xsec = 0, color=8,  reweightzpt = False, legend="DY#rightarrow#tau#tau", name="DYTauTau", prepareVars=False, select=False)
        WW2LSamp.donormalization = True
        WZ2LSamp.donormalization = True
        ZZ2LSamp.donormalization = True
        ZZ2L2QSamp.donormalization = True
        DYTauSamp.donormalization = True
        WW2LSamp.varyQCDScale = False
        WZ2LSamp.varyQCDScale = False
        ZZ2LSamp.varyQCDScale = False
        ZZ2L2QSamp.varyQCDScale = False
        DYTauSamp.varyQCDScale = False
        samples = [ WW2LSamp, WZ2LSamp, ZZ2LSamp, ZZ2L2QSamp, DYTauSamp]

    ptbins = getPtBins()
    njetbins = getJetBins()

    print(ptbins)
    print(njetbins)

    # prepare u1 (uparal) and u2 (uperp) distributions 
    # in differen Z pt and jet multiplicity bins
    def prepareU1U2(rdf, postfix="", extra_weight = "1.0"):
        histos_u1 = OrderedDict() 
        histos_u2 = OrderedDict()
        for ijet in range(njetbins.size-1):
            njetmin = njetbins[ijet]
            njetmax = njetbins[ijet+1]
            histos_u1[(njetmin, njetmax)] = OrderedDict()
            histos_u2[(njetmin, njetmax)] = OrderedDict()
            for ipt in range(ptbins.size-1):
                ptmin = ptbins[ipt]
                ptmax = ptbins[ipt+1]
    
                wstring = "njetbin_{}_ptbin_{}_{}".format(ijet, ipt, postfix)
                rdf = rdf.Define(wstring, "(jet_n>={} && jet_n<={}) * (Z_pt>={} && Z_pt<{} ) * weight * norm * {}".format(njetmin, njetmax, ptmin, ptmax, extra_weight))
                hname_u1 = "hist_uparal_{}".format(wstring)
                histos_u1[(njetmin,njetmax)][(ptmin, ptmax)] = rdf.Histo1D( (hname_u1, hname_u1, 240+int(ptmax)-int(ptmin), -120.0+int(ptmin), 120.0+int(ptmax)), "u1",  wstring )
                hname_u2 = "hist_uperp_{}".format(wstring)
                histos_u2[(njetmin, njetmax)][(ptmin, ptmax)] = rdf.Histo1D( (hname_u2, hname_u2, 240, -120.0,       120.0),       "u2",  wstring )

        return histos_u1, histos_u2


    for samp in samples:
        print("\n\n\nWork on sample ", samp.name)

        weights = OrderedDict()
        if samp.varyQCDScale:
            # do the QCD scale variations
            # Skip mur_2_muf_0p5: EventWeights[5] and
            # and mur_0p5_muf_2: EventWeights[7]
            weights["weight_mur_1_muf_1"   ]  =  "EventWeights[0] / EventWeights[0]"
            #weights["weight_mur_1_muf_2"   ]  =  "EventWeights[1] / EventWeights[0]"
            #weights["weight_mur_1_muf_0p5" ]  =  "EventWeights[2] / EventWeights[0]"
            #weights["weight_mur_2_muf_1"   ]  =  "EventWeights[3] / EventWeights[0]"
            #weights["weight_mur_2_muf_2"   ]  =  "EventWeights[4] / EventWeights[0]"
            #weights["weight_mur_0p5_muf_1" ]  =  "EventWeights[6] / EventWeights[0]"
            #weights["weight_mur_0p5_muf_0p5"] =  "EventWeights[8] / EventWeights[0]"
        else:
            weights["central"] = "1.0"

        histos_u1 = OrderedDict()
        histos_u2 = OrderedDict()
        print("\n\nMaking histograms...")
        for wname, wstring in weights.items():
            # making all histograms first, such that rdataframe could handle
            # everything in one loop.
            print(wname, wstring)
            samp.rdf = samp.rdf.Define(wname, wstring)
            histos_u1[wname], histos_u2[wname] = prepareU1U2(samp.rdf, postfix="{}_{}".format(samp.name, wname), extra_weight=wname)

        #if samp.donormalization:
        #    print("\n\nScale the {} MC to the xsec with normalization factor {}".format(samp.name, samp.normfactor))
        #    for wname in histos_u1:
        #        for ijetbin in histos_u1[wname]:
        #            for iptbin in histos_u1[wname][ijetbin]:
        #                histos_u1[wname][ijetbin][iptbin].Scale(samp.normfactor)
        #                histos_u2[wname][ijetbin][iptbin].Scale(samp.normfactor)


        print("\n\nWriting to output file...")
        for wname, wstring in weights.items():
            outdir = "results/U1U2"
            if os.path.exists(outdir) == False:
                os.makedirs(outdir)
            ofilename = "results/U1U2/histos_u1u2_{}_njets_pt_{}.root".format(samp.name, wname)

            print("Write to output file ", ofilename)
            ofile = ROOT.TFile(ofilename, "recreate")
            for ijetbin in list(histos_u1[wname].keys()):
                for iptbin in list(histos_u1[wname][ijetbin].keys()):
                    histos_u1[wname][ijetbin][iptbin].GetValue().SetDirectory(ofile)
                    histos_u1[wname][ijetbin][iptbin].GetValue().Write()
                    histos_u2[wname][ijetbin][iptbin].GetValue().SetDirectory(ofile)
                    histos_u2[wname][ijetbin][iptbin].GetValue().Write()

        ofile.Close()
            
    print("Finished..")

    return 

if __name__ == "__main__":
   main()
