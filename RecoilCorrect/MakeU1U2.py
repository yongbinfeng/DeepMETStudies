import ROOT
import numpy as np
from collections import OrderedDict
from SampleManager import Sample
import CMS_lumi
import pickle
from Utils.utils import getPtBins, getJetBins

isData     = False
isAMCNLO   = True
isMADGRAPH = False
isTTbar    = False
isDiboson  = False

assert isData+isAMCNLO+isMADGRAPH+isTTbar+isDiboson==1, "must pick one sample from data, or amc@nlo, or madgraph"

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT()

def main():
    
    if isData:
        inputfile = "inputs/inputs_Z/input_data.txt"
        samp = Sample(inputfile, isMC=False, name = "Data")
        samp.donormalization = False
        samp.varyQCDScale = False
        samples = [samp]
    elif isAMCNLO:
        inputfile = "inputs/inputs_Z/input_zjets.txt"
        samp = Sample(inputfile, isMC=True,  name = "ZJets_NLO")
        samp.donormalization = False
        samp.varyQCDScale = True
        samples = [samp]
    elif isMADGRAPH:
        inputfile = "inputs/inputs_Z/input_zjets_madgraph.txt"
        samp = Sample(inputfile, isMC=True,  name = "ZJets_MG")
        samp.donormalization = False
        samp.varyQCDScale = False
        samples = [samp]
    elif isTTbar:
        inputfile = "inputs/inputs_Z/input_ttbar.txt"
        samp = Sample(inputfile, isMC=True,  name = "TTbar", xsec = 831.76*0.105*1e3)
        samp.donormalization = True
        samp.varyQCDScale = True
        samples = [samp]
    elif isDiboson:
        input_WW2L    = "inputs/inputs_Z/input_WWTo2L2Nu.txt"
        input_WZ3L    = "inputs/inputs_Z/input_WZTo3LNu.txt"
        input_ZZ2L    = "inputs/inputs_Z/input_ZZTo2L2Nu.txt"
        input_ZZ4L    = "inputs/inputs_Z/input_ZZTo4L.txt"
        WW2LSamp  = Sample(input_WW2L,  xsec = 12.178*1e3,   name="WW2L")
        WZ3LSamp  = Sample(input_WZ3L,  xsec = 5.26*1e3,     name="WZ3L")
        ZZ2LSamp  = Sample(input_ZZ2L,  xsec = 0.564*1e3,    name="ZZ2L")
        ZZ4LSamp  = Sample(input_ZZ4L,  xsec = 1.212*1e3,    name="ZZ4L")
        WW2LSamp.donormalization = True
        WZ3LSamp.donormalization = True
        ZZ2LSamp.donormalization = True
        ZZ4LSamp.donormalization = True
        WW2LSamp.varyQCDScale = False
        WZ3LSamp.varyQCDScale = False
        ZZ2LSamp.varyQCDScale = False
        ZZ4LSamp.varyQCDScale = False
        samples = [ WW2LSamp, WZ3LSamp, ZZ2LSamp, ZZ4LSamp ]

    ptbins = getPtBins()
    njetbins = getJetBins()

    print ptbins
    print njetbins

    # prepare u1 (uparal) and u2 (uperp) distributions 
    # in differen Z pt and jet multiplicity bins
    def prepareU1U2(rdf, postfix="", extra_weight = "1.0"):
        histos_u1 = OrderedDict() 
        histos_u2 = OrderedDict()
        for ijet in xrange(njetbins.size-1):
            njetmin = njetbins[ijet]
            njetmax = njetbins[ijet+1]
            histos_u1[(njetmin, njetmax)] = OrderedDict()
            histos_u2[(njetmin, njetmax)] = OrderedDict()
            for ipt in xrange(ptbins.size-1):
                ptmin = ptbins[ipt]
                ptmax = ptbins[ipt+1]
    
                wstring = "njetbin_{}_ptbin_{}_{}".format(ijet, ipt, postfix)
                rdf = rdf.Define(wstring, "(jet_n>={} && jet_n<={}) * (Z_pt>={} && Z_pt<{} ) * weight * {}".format(njetmin, njetmax, ptmin, ptmax, extra_weight))
                hname_u1 = "hist_uparal_{}".format(wstring)
                histos_u1[(njetmin,njetmax)][(ptmin, ptmax)] = rdf.Histo1D( (hname_u1, hname_u1, 240+int(ptmax)-int(ptmin), -120.0+int(ptmin), 120.0+int(ptmax)), "u1",  wstring )
                hname_u2 = "hist_uperp_{}".format(wstring)
                histos_u2[(njetmin, njetmax)][(ptmin, ptmax)] = rdf.Histo1D( (hname_u2, hname_u2, 240, -120.0,       120.0),       "u2",  wstring )

        return histos_u1, histos_u2


    for samp in samples:
        print "\n\n\nWork on sample ", samp.name

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
        print "\n\nMaking histograms..."
        for wname, wstring in weights.iteritems():
            # making all histograms first, such that rdataframe could handle
            # everything in one loop.
            print wname, wstring
            samp.rdf = samp.rdf.Define(wname, wstring)
            histos_u1[wname], histos_u2[wname] = prepareU1U2(samp.rdf, postfix="{}_{}".format(samp.name, wname), extra_weight=wname)

        if samp.donormalization:
            print "\n\nScale the {} MC to the xsec with normalization factor {}".format(samp.name, samp.normfactor)
            for wname in histos_u1:
                for ijetbin in histos_u1[wname]:
                    for iptbin in histos_u1[wname][ijetbin]:
                        histos_u1[wname][ijetbin][iptbin].Scale(samp.normfactor)
                        histos_u2[wname][ijetbin][iptbin].Scale(samp.normfactor)


        print "\n\nWriting to output file..."
        for wname, wstring in weights.iteritems():
            ofilename = "results/U1U2/histos_u1u2_{}_njets_pt_{}.root".format(samp.name, wname)

            print "Write to output file ", ofilename
            ofile = ROOT.TFile(ofilename, "recreate")
            for ijetbin in histos_u1[wname].keys():
                for iptbin in histos_u1[wname][ijetbin].keys():
                    histos_u1[wname][ijetbin][iptbin].GetValue().SetDirectory(ofile)
                    histos_u1[wname][ijetbin][iptbin].GetValue().Write()
                    histos_u2[wname][ijetbin][iptbin].GetValue().SetDirectory(ofile)
                    histos_u2[wname][ijetbin][iptbin].GetValue().Write()

        ofile.Close()
            
    print "Finished.."

    return 

if __name__ == "__main__":
   main()
