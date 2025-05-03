import ROOT
import sys
from CMSPLOTS.myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code
from utils.RecoilAnalyzer import RecoilAnalyzer
import argparse

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("Functions_cc.so")

parser = argparse.ArgumentParser()
parser.add_argument("--era", default="2016", help="Era")
args = parser.parse_args()

do2016 = (args.era == "2016")
do2017 = (args.era == "2017")
do2018 = (args.era == "2018")

era = args.era
outdir = f"plots/MC/{era}"

rdf_org = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/group/cmst3/group/wmass/w-mass-13TeV/NanoAOD/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_testPVRobustDM/240229_074036/0000/NanoV9MCPostVFP_10.root")
rdf_org1 = rdf_org.Filter("nMuon >= 1")
rdf_org1 = rdf_org1.Define("Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]")
rdf_org2 = rdf_org1.Filter("Muon_pass0")
rdf = rdf_org2

rdf = rdf.Define("muon_pt", "Muon_pt[0]").Define("muon_phi", "Muon_phi[0]")

rdf = rdf.Define("PFMET_pt", "MET_pt") \
         .Define("PUPPIMET_pt", "PuppiMET_pt") \
         .Define("DeepMETMET_pt", "DeepMETResolutionTune_pt")

rdf= rdf.Define("mT_PF", "calMT_fromMET_PtPhi(muon_pt, muon_phi, MET_pt, MET_phi)")
rdf= rdf.Define("mT_PUPPI", "calMT_fromMET_PtPhi(muon_pt, muon_phi, PuppiMET_pt, PuppiMET_phi)")
rdf= rdf.Define("mT_DeepMET", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")
rdf= rdf.Define("mT_DeepMETPVRobust", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETPVRobust_pt, DeepMETPVRobust_phi)")
rdf= rdf.Define("mT_DeepMETPVRobustNoPUPPI", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETPVRobustNoPUPPI_pt, DeepMETPVRobustNoPUPPI_phi)")

rdf_Z_org = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/group/cmst3/group/wmass/yofeng/NanoAOD/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_testPVRobustDM/240226_075056/0000/NanoV9MCPostVFP_1*.root")
rdf_Z_org1 = rdf_Z_org.Filter("nMuon > 1")
rdf_Z_org1 = rdf_Z_org1.Define("Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]")
rdf_Z_org1 = rdf_Z_org1.Define("Muon_pass1", "Muon_pt[1] > 25.0 && abs(Muon_eta[1]) < 2.4 && Muon_pfRelIso04_all[1] < 0.15 && Muon_looseId[1]")

rdf_Z_org2 = rdf_Z_org1.Filter("Muon_pass0 && Muon_pass1")
rdf_Z = rdf_Z_org2

h_pvIndex = rdf.Histo1D(("pv_Index", "pv_Index", 11, -1.5, 10.5), "PVRobustIndex")
h_pvIndex_Z = rdf_Z.Histo1D(("pv_Index", "pv_Index", 11, -1.5, 10.5), "PVRobustIndex")

h_MTs = OrderedDict()
h_METs = OrderedDict()
#for met in ["PF", "DeepMET", "DeepMETPVRobust", "DeepMETPVRobustNoPUPPI"]:
for met in ["PF", "PUPPI", "DeepMET"]:
    hname = "h_mT_{}".format(met)
    h_MTs[met] = rdf.Histo1D((hname, hname, 40, 0, 120), "mT_{}".format(met))
    hname = "h_MET_{}".format(met)
    h_METs[met] = rdf.Histo1D((hname, hname, 40, 0, 100), "{}MET_pt".format(met))

colors = {
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "DeepMET": 4,
            "DeepMETPVRobust": 6,
            "DeepMETPVRobustNoPUPPI": 7
         }

labels = {
            "TK": "TK",
            "PF": "PF",
            "PUPPI": "PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DeepMET": "DeepMET",
            "DeepMETPVRobust": "DeepMET PVRobust",
            "DeepMETPVRobustNoPUPPI": "DeepMET PVRobust NoPUPPI"
         }

DrawHistos([h_pvIndex, h_pvIndex_Z], ["W", "Z"], -1.5, 10.5, "Robust PV Index", 1e-5, 10.0, "a.u.", "W_pv_Index", donormalize=True, outdir=outdir, noLumi=True, mycolors=[1, 2])

DrawHistos(h_MTs.values(), [labels[itype] for itype in h_MTs.keys()], 0, 120, "m_{T} [GeV]", 0., 0.1, "a.u.", "MT_comp", drawashist=True, dology=False, legendPos=[0.20, 0.67, 0.38, 0.86], mycolors=[colors[itype] for itype in h_MTs.keys()], outdir=outdir, donormalize=True, MCOnly=True)

DrawHistos(h_METs.values(), [labels[itype] for itype in h_METs.keys()], 0, 100, "p^{miss}_{T} [GeV]", 0., 0.13, "a.u.", "MET_comp", drawashist=True, dology=False, legendPos=[0.62, 0.67, 0.80, 0.86], mycolors=[colors[itype] for itype in h_METs.keys()], outdir=outdir, donormalize=True, MCOnly=True)

