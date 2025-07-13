import sys
sys.path.append("../RecoilResol/CMSPLOTS")  # noqa
import argparse
from utils.RecoilAnalyzer import RecoilAnalyzer
from utils.utils import getpTBins, getnVtxBins, get_response_code, prepRecoilVars, getqTRange, getqTLabel, getnVtxLabel, getResponseLabel, getUparalLabel, getUperpLabel, getVtxRange, getpTResponseBins, doPAS, getNPVString
from collections import OrderedDict
from CMSPLOTS.myFunction import DrawHistos
import ROOT

noLumi = False
MCOnly = False
drawashist = False
doPAS = doPAS()

nPV = getNPVString()

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("Functions_cc.so")

parser = argparse.ArgumentParser()
parser.add_argument("--era", default="2016", help="Era")
parser.add_argument("--applySc", action="store_true",
                    help="Apply response corrections and include the plots in the output")
parser.add_argument("--no-applySc", action="store_false", dest="applySc",
                    help="Do not apply response corrections and do not include the plots in the output")
parser.set_defaults(applySc=True)

args = parser.parse_args()

doTest = (args.era == "test")
do2016 = (args.era == "2016")

assert doTest + do2016 == 1, "Please specify an era: test, 2016, 2017, 2018"

era = args.era
outdir = f"plots/MC_CompPUPPI/{era}"

applySc = args.applySc
print("apply Response corrections: ", applySc)

# mcfile = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/DY.root"
datafile = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/Data.root"
chain = ROOT.TChain("Events")
chain.Add(datafile)
rdf = ROOT.ROOT.RDataFrame(chain)


bkgfiles = [
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/ttbar.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/DYTauTau.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/WW2L.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/WZ2L.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/ZZ2L.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/ZZ2L2Q.root",
]

chainBkg = ROOT.TChain("Events")
for bkgfile in bkgfiles:
    chainBkg.Add(bkgfile)
rdf_bkg = ROOT.ROOT.RDataFrame(chainBkg)


def prepareVars(rdf):
    rdf = rdf.Define("pT_muons", "Z_pt").Define("phi_muons", "Z_phi")
    rdf = rdf.Define("u_GEN_pt", "pT_muons") \
             .Define("u_GEN_x", "-pT_muons*TMath::Cos(phi_muons)") \
             .Define("u_GEN_y", "-pT_muons*TMath::Sin(phi_muons)") \
             .Define("u_PUPPI_x",  "-(pT_muons*TMath::Cos(phi_muons) + PuppiMET_pt*TMath::Cos(PuppiMET_phi))") \
             .Define("u_PUPPI_y",  "-(pT_muons*TMath::Sin(phi_muons) + PuppiMET_pt*TMath::Sin(PuppiMET_phi))") \
             .Define("u_PUPPI_pt", "TMath::Sqrt(u_PUPPI_x * u_PUPPI_x + u_PUPPI_y * u_PUPPI_y)") \
             .Define("u_PF_x",     "-(pT_muons*TMath::Cos(phi_muons) + MET_pt*TMath::Cos(MET_phi))") \
             .Define("u_PF_y",     "-(pT_muons*TMath::Sin(phi_muons) + MET_pt*TMath::Sin(MET_phi))") \
             .Define("u_PF_pt",    "TMath::Sqrt(u_PF_x * u_PF_x + u_PF_y * u_PF_y)") \
             .Define("u_DeepMET_x",  "-(pT_muons*TMath::Cos(phi_muons) + DeepMETResolutionTune_pt*TMath::Cos(DeepMETResolutionTune_phi))") \
             .Define("u_DeepMET_y",  "-(pT_muons*TMath::Sin(phi_muons) + DeepMETResolutionTune_pt*TMath::Sin(DeepMETResolutionTune_phi))") \
             .Define("u_DeepMET_pt", "TMath::Sqrt(u_DeepMET_x * u_DeepMET_x + u_DeepMET_y * u_DeepMET_y)") \
             .Define("u_DeepMETNoPUPPI_x", "-(pT_muons*TMath::Cos(phi_muons) + DeepMETPVRobustNoPUPPI_pt*TMath::Cos(DeepMETPVRobustNoPUPPI_phi) )") \
             .Define("u_DeepMETNoPUPPI_y", "-(pT_muons*TMath::Sin(phi_muons) + DeepMETPVRobustNoPUPPI_pt*TMath::Sin(DeepMETPVRobustNoPUPPI_phi) )") \
             .Define("u_DeepMETNoPUPPI_pt", "TMath::Sqrt(u_DeepMETNoPUPPI_x * u_DeepMETNoPUPPI_x + u_DeepMETNoPUPPI_y * u_DeepMETNoPUPPI_y)")

    # rdf = rdf.Define("u_DeepMETCorr_x", "-(pT_muons*TMath::Cos(phi_muons) + deepmet_pt_corr_central*TMath::Cos(deepmet_phi_corr_central) )") \
    #         .Define("u_DeepMETCorr_y", "-(pT_muons*TMath::Sin(phi_muons) + deepmet_pt_corr_central*TMath::Sin(deepmet_phi_corr_central) )") \
    #         .Define("u_DeepMETCorr_pt", "TMath::Sqrt(u_DeepMETCorr_x * u_DeepMETCorr_x + u_DeepMETCorr_y * u_DeepMETCorr_y)") \
    return rdf


rdf = prepareVars(rdf)
rdf_bkg = prepareVars(rdf_bkg)

rdf = rdf.Define("weight_dummy", "1.0").Define(
    "weight_withBtag", "weight_corr * (jet_CSVLoose_n < 1) * metFilters")
rdf_bkg = rdf_bkg.Define("weight_dummy", "1.0").Define(
    "weight_withBtag", "weight_corr * (jet_CSVLoose_n < 1) * metFilters")
weightname = "weight_withBtag"

recoils = ["PF", "PUPPI", "DeepMET", "DeepMETNoPUPPI"]

colors = {
    "PF": 1,
    "PUPPI": 2,
    "GEN": 6,
    "DeepMET": 4,
    "DeepMETCorr": 8,
    "DeepMETNoPUPPI": 8,
}

labels = {
    "TK": "TK",
    "PF": "PF",
    "PUPPI": "PUPPI",
    "GEN": "GEN",
    "TKPHO": "TK+Photon",
    "DeepMET": "DeepMET",
    "DeepMETCorr": "DeepMET",
    "DeepMETNoPUPPI": "DeepMET w/o PUPPI",
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
    "DeepMETNoPUPPI": 23,
}

linestyles = {}

for itype in recoils:
    colors[itype + "_MC"] = colors[itype]

    linestyles[itype] = 1
    linestyles[itype + "_MC"] = 2

xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins()
xbins_qT_resp = getpTResponseBins()

# loop over the different qT bins
suffix = "_testPUPPIEffect"
recoilanalyzer = RecoilAnalyzer(
    rdf, recoils, name="recoilanalyzer_" + suffix, useRMS=True, weightname=weightname, rdfBkg=rdf_bkg)
recoilanalyzer.prepareVars()
recoilanalyzer.prepareResponses('u_GEN_pt', xbins_qT)
recoilanalyzer.prepareResolutions('u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzer.prepareResponses(nPV, xbins_nVtx)
recoilanalyzer.prepareResolutions(nPV, xbins_nVtx, 400, -200, 200)

hresponses = recoilanalyzer.getResponses('u_GEN_pt')
hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions('u_GEN_pt')
hresponses_nVtx = recoilanalyzer.getResponses(nPV)
hresols_paral_diff_VS_nVtx, hresols_perp_VS_nVtx = recoilanalyzer.getResolutions(
    nPV)


# just one bin, to calculate the response correction
recoilanalyzer.prepareResponses('pT_muons', xbins_qT_resp)
hresponses_inclusive = recoilanalyzer.getResponses('pT_muons')

values_responses = OrderedDict()
for itype in recoils:
    print("hresponses_inclusive in data for ", itype, " is ",
          hresponses_inclusive[itype].GetBinContent(1))
    values_responses[itype] = hresponses_inclusive[itype].GetBinContent(1)

if applySc:
    hresolsSc_paral_diff = OrderedDict()
    hresolsSc_perp = OrderedDict()
    hresolsSc_paral_diff_VS_nVtx = OrderedDict()
    hresolsSc_perp_VS_nVtx = OrderedDict()

    for itype in recoils:
        resp = values_responses[itype]

        hresolsSc_paral_diff[itype] = hresols_paral_diff[itype].Clone(
            itype + "Sc_paral_diff")
        hresolsSc_paral_diff[itype].Scale(1.0 / resp)

        hresolsSc_perp[itype] = hresols_perp[itype].Clone(
            itype + "Sc_perp")
        hresolsSc_perp[itype].Scale(1.0 / resp)

        hresolsSc_paral_diff_VS_nVtx[itype] = hresols_paral_diff_VS_nVtx[itype].Clone(
            itype + "Sc_paral_diff_VS_nVtx")
        hresolsSc_paral_diff_VS_nVtx[itype].Scale(1.0 / resp)

        hresolsSc_perp_VS_nVtx[itype] = hresols_perp_VS_nVtx[itype].Clone(
            itype + "Sc_perp_VS_nVtx")
        hresolsSc_perp_VS_nVtx[itype].Scale(1.0 / resp)


def GetLegends(hdict):
    return [labels[itype] for itype in hdict.keys() if "_MC" not in itype]


def GetLineStyles(hdict):
    return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]


def GetMarkers(hdict):
    return [markers[itype] if "_MC" not in itype else 1 for itype in hdict.keys()]


def GetDrawOptions(hdict):
    return ["EP" if "_MC" not in itype else "HIST" for itype in hdict.keys()]


qtmin, qtmax = getqTRange()
qtlabel = getqTLabel()
nvtxlabel = getnVtxLabel()
uparallabel = getUparalLabel()
uperplabel = getUperpLabel()
responselabel = getResponseLabel()

nvtxmin, nvtxmax = getVtxRange()


hcolors = [colors[itype] for itype in hresponses.keys()]
hmarkers = GetMarkers(hresponses)
linestyles = GetLineStyles(hresponses)
drawoptions = GetDrawOptions(hresponses)
legends = GetLegends(hresponses)


args = {
    "mycolors": hcolors,
    "markerstyles": hmarkers,
    "linestyles": linestyles,
    "drawoptions": drawoptions,
    "outdir": outdir,
    "noLumi": noLumi,
    "MCOnly": MCOnly,
    "dology": False,
    "drawashist": False,
    "legendoptions": ["LEP"] * len(hresponses),
    "doPAS": doPAS,
}

DrawHistos(hresponses.values(), legends, 0, qtmax, qtlabel, 0., 1.29, responselabel,
           "reco_recoil_response" + suffix, legendPos=[0.50, 0.20, 0.80, 0.40], inPaper=True, **args)

DrawHistos(hresols_paral_diff.values(), GetLegends(hresols_paral_diff), 0, qtmax, qtlabel, 0, 39.0, uparallabel, "reco_recoil_resol_paral" + suffix, drawashist=drawashist, dology=False,
           legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresols_paral_diff), MCOnly=MCOnly)

DrawHistos(hresols_perp.values(), GetLegends(hresols_perp), 0, qtmax, qtlabel, 0, 32.0, uperplabel, "reco_recoil_resol_perp" + suffix, drawashist=drawashist, dology=False,
           legendPos=[0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresols_perp), MCOnly=MCOnly)

DrawHistos(hresponses_nVtx.values(), GetLegends(hresponses_nVtx), nvtxmin, nvtxmax, nvtxlabel, 0., 1.15, responselabel, "reco_recoil_response_VS_nVtx" + suffix, drawashist=drawashist, dology=False,
           legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresponses_nVtx), MCOnly=MCOnly)

DrawHistos(hresols_paral_diff_VS_nVtx.values(), GetLegends(hresols_paral_diff_VS_nVtx), nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uparallabel, "reco_recoil_resol_paral_VS_nVtx" + suffix, drawashist=drawashist, dology=False,
           legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresols_paral_diff_VS_nVtx), MCOnly=MCOnly)

DrawHistos(hresols_perp_VS_nVtx.values(), GetLegends(hresols_perp_VS_nVtx), nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel, "reco_recoil_resol_perp_VS_nVtx" + suffix, drawashist=drawashist, dology=False,
           legendPos=[0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresols_perp_VS_nVtx), MCOnly=MCOnly)

if applySc:
    extraToDraw = ROOT.TPaveText(0.60, 0.20, 0.90, 0.27, "NDC")
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText("Response corrected")

    args["extraToDraw"] = extraToDraw

    #
    # Scaled
    #
    DrawHistos(hresolsSc_paral_diff.values(), GetLegends(hresols_paral_diff), 0, qtmax, qtlabel, 0, 60.0, uparallabel, "reco_recoil_resol_paral_Scaled" + suffix, drawashist=drawashist, dology=False, legendPos=[
               0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresols_paral_diff), MCOnly=MCOnly, extraToDraw=extraToDraw)

    DrawHistos(hresolsSc_perp.values(), GetLegends(hresols_perp), 0, qtmax, qtlabel, 0, 50.0, uperplabel, "reco_recoil_resol_perp_Scaled" + suffix, drawashist=drawashist, dology=False, legendPos=[
               0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp.keys()], noLumi=noLumi, outdir=outdir, linestyles=GetLineStyles(hresols_perp), MCOnly=MCOnly, extraToDraw=extraToDraw)

    DrawHistos(hresolsSc_paral_diff_VS_nVtx.values(), GetLegends(hresols_paral_diff_VS_nVtx), nvtxmin, nvtxmax, nvtxlabel, 0, 50.0,
               uparallabel, "reco_recoil_resol_paral_VS_nVtx_Scaled" + suffix, legendPos=[0.20, 0.60, 0.45, 0.78], inPaper=True, **args)

    DrawHistos(hresolsSc_perp_VS_nVtx.values(), GetLegends(hresols_perp_VS_nVtx), nvtxmin, nvtxmax, nvtxlabel, 0, 50.0,
               uperplabel, "reco_recoil_resol_perp_VS_nVtx_Scaled" + suffix, legendPos=[0.20, 0.60, 0.45, 0.78], inPaper=True, **args)


# f1 = ROOT.TFile("root_h/output.root", "RECREATE")
# for h2 in h2ds_perp_VS_qT.values():
#    h2.SetDirectory(f1)
#    h2.Write()
# f1.Close()
