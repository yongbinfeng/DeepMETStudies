import argparse
import sys  # noqa
sys.path.append("../RecoilResol/CMSPLOTS")  # noqa
from utils.RecoilAnalyzer import RecoilAnalyzer
from utils.utils import getpTBins, getnVtxBins, get_response_code, doPAS
from collections import OrderedDict
from CMSPLOTS.myFunction import DrawHistos
import ROOT

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("Functions_cc.so")

parser = argparse.ArgumentParser()
parser.add_argument("--era", default="2016", help="Era")
args = parser.parse_args()

do2016 = (args.era == "2016")
do2017 = (args.era == "2017")
do2018 = (args.era == "2018")

doPAS = doPAS()
doData = True

era = args.era
outdir = f"plots/MC/{era}"

if not doData:
    fname = "/home/yongbinfeng/Desktop/DeepMET/data/wjetsMC/NanoV9MCPostVFP_294.root"
    rdf_org = ROOT.ROOT.RDataFrame("Events", fname)
    rdf_org1 = rdf_org.Filter("nMuon >= 1")
    rdf_org1 = rdf_org1.Define(
        "Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]")
    rdf_org2 = rdf_org1.Filter("Muon_pass0")
    rdf = rdf_org2

    rdf = rdf.Define("muon_pt", "Muon_pt[0]").Define("muon_phi", "Muon_phi[0]")

else:
    fname = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_wjets/Data.root"
    rdf = ROOT.ROOT.RDataFrame("Events", fname)
    # event selections are already applied in the slim
    rdf = rdf.Define("muon_pt", "leadMuon_pt") \
             .Define("muon_phi", "leadMuon_phi")

rdf = rdf.Define("PFMET_pt", "MET_pt") \
    .Define("PUPPIMET_pt", "PuppiMET_pt") \
    .Define("DeepMETMET_pt", "DeepMETResolutionTune_pt")

rdf = rdf.Define(
    "mT_PF", "calMT_fromMET_PtPhi(muon_pt, muon_phi, MET_pt, MET_phi)")
rdf = rdf.Define(
    "mT_PUPPI", "calMT_fromMET_PtPhi(muon_pt, muon_phi, PuppiMET_pt, PuppiMET_phi)")
rdf = rdf.Define(
    "mT_DeepMET", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")
rdf = rdf.Define("mT_DeepMETPVRobust",
                 "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETPVRobust_pt, DeepMETPVRobust_phi)")
rdf = rdf.Define("mT_DeepMETPVRobustNoPUPPI",
                 "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETPVRobustNoPUPPI_pt, DeepMETPVRobustNoPUPPI_phi)")

# rdf_Z_org = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/group/cmst3/group/wmass/yofeng/NanoAOD/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_testPVRobustDM/240226_075056/0000/NanoV9MCPostVFP_1*.root")
# rdf_Z_org1 = rdf_Z_org.Filter("nMuon > 1")
# rdf_Z_org1 = rdf_Z_org1.Define("Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]")
# rdf_Z_org1 = rdf_Z_org1.Define("Muon_pass1", "Muon_pt[1] > 25.0 && abs(Muon_eta[1]) < 2.4 && Muon_pfRelIso04_all[1] < 0.15 && Muon_looseId[1]")
#
# rdf_Z_org2 = rdf_Z_org1.Filter("Muon_pass0 && Muon_pass1")
# rdf_Z = rdf_Z_org2
#
# h_pvIndex = rdf.Histo1D(
#    ("pv_Index", "pv_Index", 11, -1.5, 10.5), "PVRobustIndex")
# h_pvIndex_Z = rdf_Z.Histo1D(("pv_Index", "pv_Index", 11, -1.5, 10.5), "PVRobustIndex")

h_MTs = OrderedDict()
h_MTs_FineBins = OrderedDict()
h_METs = OrderedDict()
h_METs_FineBins = OrderedDict()
# for met in ["PF", "DeepMET", "DeepMETPVRobust", "DeepMETPVRobustNoPUPPI"]:
for met in ["PF", "PUPPI", "DeepMET"]:
    hname = "h_mT_{}".format(met)
    h_MTs[met] = rdf.Histo1D((hname, hname, 60, 0, 120), "mT_{}".format(met))
    hname_fine = "h_mT_{}_fine".format(met)
    h_MTs_FineBins[met] = rdf.Histo1D(
        (hname_fine, hname_fine, 600, 0, 120), "mT_{}".format(met))
    hname = "h_MET_{}".format(met)
    h_METs[met] = rdf.Histo1D(
        (hname, hname, 50, 0, 100), "{}MET_pt".format(met))
    hname_fine = "h_MET_{}_fine".format(met)
    h_METs_FineBins[met] = rdf.Histo1D(
        (hname_fine, hname_fine, 500, 0, 100), "{}MET_pt".format(met))

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


markers = {
    "PF": 20,
    "PUPPI": 21,
    "GEN": 22,
    "DeepMET": 22,
    "DeepMETCorr": 23,
    "PUPPIUnc": 24,
    "RawPF": 25,
    "RawPUPPI": 26,
    "PFUp": 25,
    "PFDown": 26,
}

# DrawHistos([h_pvIndex], ["W", "Z"], -1.5, 10.5, "Robust PV Index", 1e-5, 10.0, "a.u.",
#           "W_pv_Index", donormalize=True, outdir=outdir, noLumi=True, mycolors=[1, 2])


def FindFWHM(histo):
    max_bin = histo.GetMaximumBin()
    max_value = histo.GetBinContent(max_bin)
    half_max = max_value / 2.0
    left_bin = max_bin
    right_bin = max_bin

    while left_bin > 1 and histo.GetBinContent(left_bin) > half_max:
        left_bin -= 1
    while right_bin < histo.GetNbinsX() and histo.GetBinContent(right_bin) > half_max:
        right_bin += 1

    # fwhm = histo.GetBinLowEdge(right_bin) - histo.GetBinLowEdge(left_bin)
    fwhmRight = histo.GetBinLowEdge(right_bin) - histo.GetBinLowEdge(max_bin)
    return fwhmRight


DrawHistos(h_MTs.values(), [labels[itype] for itype in h_MTs.keys()], 0, 120, "m_{T} [GeV]", 0., 9.9e6, "Events / 2 Gev", "MT_comp", drawashist=False, dology=False, legendPos=[
           0.20, 0.63, 0.48, 0.82], mycolors=[colors[itype] for itype in h_MTs.keys()], outdir=outdir, donormalize=False, MCOnly=not doData, doPAS=doPAS, inPaper=True, nMaxDigits=3, markerstyles=[markers[itype] for itype in h_MTs.keys()], legendoptions=['LEP'] * len(h_MTs.keys()))

DrawHistos(h_METs.values(), [labels[itype] for itype in h_METs.keys()], 0, 100, "p^{miss}_{T} [GeV]", 0., 1.5e7, "Events / 2 GeV", "MET_comp", drawashist=False, dology=False, legendPos=[
           0.62, 0.63, 0.90, 0.82], mycolors=[colors[itype] for itype in h_METs.keys()], outdir=outdir, donormalize=False, MCOnly=not doData, doPAS=doPAS, inPaper=True, markerstyles=[markers[itype] for itype in h_METs.keys()], legendoptions=['LEP'] * len(h_METs.keys()))

for met in h_MTs_FineBins.keys():
    fwhm = FindFWHM(h_MTs_FineBins[met])
    print(f"FWHM for {met}: {fwhm:.2f} GeV")

for met in h_METs_FineBins.keys():
    fwhm = FindFWHM(h_METs_FineBins[met])
    print(f"FWHM for {met}: {fwhm:.2f} GeV")

ofile = ROOT.TFile(f"{outdir}/recoil_WJets_{era}.root", "RECREATE")
for h in h_MTs_FineBins.values():
    h.Write()
for h in h_METs_FineBins.values():
    h.Write()
ofile.Close()
