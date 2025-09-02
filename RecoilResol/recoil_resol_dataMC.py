import sys
sys.path.append("../RecoilResol/CMSPLOTS")  # noqa
import argparse
from utils.RecoilAnalyzer import RecoilAnalyzer
from utils.utils import getpTBins, getpTResponseBins, getnVtxBins, get_response_code, getqTLabel, getnVtxLabel, getNPVString, getqTRange, getResponseLabel, getUparalLabel, getUperpLabel, getVtxRange, doPAS
from collections import OrderedDict
from CMSPLOTS.myFunction import DrawHistos
import ROOT
import math

noLumi = False
nPV = getNPVString()
useRMS = True

# whether to plot data and MC in the same plot
combineDataMC = True

doPAS = doPAS()

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()
ROOT.gSystem.Load("Functions_cc.so")

print("start processing")

parser = argparse.ArgumentParser()
parser.add_argument("--applySc", action="store_true",
                    help="Apply response corrections and include the plots in the output")
parser.add_argument("--no-applySc", action="store_false", dest="applySc",
                    help="Do not apply response corrections and do not include the plots in the output")
parser.add_argument("--test", action="store_true",
                    help="Run the script in test mode")
parser.set_defaults(applySc=True)

args = parser.parse_args()

outdir = f"plots/Data/2016"

mcfile = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/DY.root"
datafile = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/Data.root"
bkgfiles = [
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/ttbar.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/DYTauTau.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/WW2L.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/WZ2L.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/ZZ2L.root",
    "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs/ZZ2L2Q.root",
]

applySc = args.applySc
print("apply Response corrections: ", applySc)

chain = ROOT.TChain("Events")
chain.Add(datafile)
rdf_data = ROOT.ROOT.RDataFrame(chain)
rdf_data = rdf_data.Define("MET_pt_smeared", "MET_pt").Define("MET_phi_smeared", "MET_phi") \
    .Define("MET_pt_smeared_Up", "MET_pt_smeared").Define("MET_phi_smeared_Up", "MET_phi_smeared") \
    .Define("MET_pt_smeared_Down", "MET_pt_smeared").Define("MET_phi_smeared_Down", "MET_phi_smeared") \

chainMC = ROOT.TChain("Events")
chainMC.Add(mcfile)
rdf_MC = ROOT.ROOT.RDataFrame(chainMC)

chainBkg = ROOT.TChain("Events")
for bkgfile in bkgfiles:
    chainBkg.Add(bkgfile)
rdf_bkg = ROOT.ROOT.RDataFrame(chainBkg)

rdf_data = rdf_data.Define("MET_pt_fixed", "MET_pt_smeared > 500 ? 500 : MET_pt_smeared") \
    .Define("MET_pt_fixed_Up", "MET_pt_smeared_Up > 500 ? 500 : MET_pt_smeared_Up") \
    .Define("MET_pt_fixed_Down", "MET_pt_smeared_Down > 500 ? 500 : MET_pt_smeared_Down")
rdf_MC = rdf_MC.Define("MET_pt_fixed", "MET_pt_smeared > 500 ? 500 : MET_pt_smeared") \
    .Define("MET_pt_fixed_Up", "MET_pt_smeared_Up > 500 ? 500 : MET_pt_smeared_Up") \
    .Define("MET_pt_fixed_Down", "MET_pt_smeared_Down > 500 ? 500 : MET_pt_smeared_Down")
rdf_bkg = rdf_bkg.Define("MET_pt_fixed", "MET_pt_smeared > 500 ? 500 : MET_pt_smeared") \
    .Define("MET_pt_fixed_Up", "MET_pt_smeared_Up > 500 ? 500 : MET_pt_smeared_Up") \
    .Define("MET_pt_fixed_Down", "MET_pt_smeared_Down > 500 ? 500 : MET_pt_smeared_Down")


# weight_corr includes both the pileup and Zpt corrections
weightname = "weight_corr"

uncs = ["JES", "JER", "Unclustered"]

if args.test:
    # need some temporary variables to hold the original rdf
    # otherise might crash
    rdf_data_tmp = rdf_data
    rdf_MC_tmp = rdf_MC
    rdf_bkg_tmp = rdf_bkg

    rdf_data = rdf_data_tmp.Filter("rdfentry_ < 10000")
    rdf_MC = rdf_MC_tmp.Filter("rdfentry_ < 10000")
    rdf_bkg = rdf_bkg_tmp.Filter("rdfentry_ < 10000")


def prepareVars(rdf):
    rdf = rdf.Define("pT_muons", "Z_pt").Define("phi_muons", "Z_phi")
    rdf = rdf.Define("u_GEN_pt", "pT_muons") \
             .Define("u_GEN_x", "-pT_muons*TMath::Cos(phi_muons)") \
             .Define("u_GEN_y", "-pT_muons*TMath::Sin(phi_muons)") \
             .Define("u_PUPPI_x",  "-(pT_muons*TMath::Cos(phi_muons) + PuppiMET_pt*TMath::Cos(PuppiMET_phi))") \
             .Define("u_PUPPI_y",  "-(pT_muons*TMath::Sin(phi_muons) + PuppiMET_pt*TMath::Sin(PuppiMET_phi))") \
             .Define("u_PUPPI_pt", "TMath::Sqrt(u_PUPPI_x * u_PUPPI_x + u_PUPPI_y * u_PUPPI_y)") \
             .Define("u_PF_x",     "-(pT_muons*TMath::Cos(phi_muons) + MET_pt_fixed*TMath::Cos(MET_phi_smeared))") \
             .Define("u_PF_y",     "-(pT_muons*TMath::Sin(phi_muons) + MET_pt_fixed*TMath::Sin(MET_phi_smeared))") \
             .Define("u_PF_pt",    "TMath::Sqrt(u_PF_x * u_PF_x + u_PF_y * u_PF_y)") \
             .Define("u_PFUp_x",  "-(pT_muons*TMath::Cos(phi_muons) + MET_pt_fixed_Up*TMath::Cos(MET_phi_smeared_Up))") \
             .Define("u_PFUp_y",  "-(pT_muons*TMath::Sin(phi_muons) + MET_pt_fixed_Up*TMath::Sin(MET_phi_smeared_Up))") \
             .Define("u_PFUp_pt", "TMath::Sqrt(u_PFUp_x * u_PFUp_x + u_PFUp_y * u_PFUp_y)") \
             .Define("u_PFDown_x",  "-(pT_muons*TMath::Cos(phi_muons) + MET_pt_fixed_Down*TMath::Cos(MET_phi_smeared_Down))") \
             .Define("u_PFDown_y",  "-(pT_muons*TMath::Sin(phi_muons) + MET_pt_fixed_Down*TMath::Sin(MET_phi_smeared_Down))") \
             .Define("u_PFDown_pt", "TMath::Sqrt(u_PFDown_x * u_PFDown_x + u_PFDown_y * u_PFDown_y)") \
             .Define("u_DeepMET_x",  "-(pT_muons*TMath::Cos(phi_muons) + DeepMETResolutionTune_pt*TMath::Cos(DeepMETResolutionTune_phi))") \
             .Define("u_DeepMET_y",  "-(pT_muons*TMath::Sin(phi_muons) + DeepMETResolutionTune_pt*TMath::Sin(DeepMETResolutionTune_phi))") \
             .Define("u_DeepMET_pt", "TMath::Sqrt(u_DeepMET_x * u_DeepMET_x + u_DeepMET_y * u_DeepMET_y)") \
             .Define("u_RawPF_x", "-(pT_muons*TMath::Cos(phi_muons) + RawMET_pt*TMath::Cos(RawMET_phi))") \
             .Define("u_RawPF_y", "-(pT_muons*TMath::Sin(phi_muons) + RawMET_pt*TMath::Sin(RawMET_phi))") \
             .Define("u_RawPF_pt", "TMath::Sqrt(u_RawPF_x * u_RawPF_x + u_RawPF_y * u_RawPF_y)") \
             .Define("u_RawPUPPI_x", "-(pT_muons*TMath::Cos(phi_muons) + RawPuppiMET_pt*TMath::Cos(RawPuppiMET_phi))") \
             .Define("u_RawPUPPI_y", "-(pT_muons*TMath::Sin(phi_muons) + RawPuppiMET_pt*TMath::Sin(RawPuppiMET_phi))") \
             .Define("u_RawPUPPI_pt", "TMath::Sqrt(u_RawPUPPI_x * u_RawPUPPI_x + u_RawPUPPI_y * u_RawPUPPI_y)")

    # .Define("u_DeepMETCorr_x", "-(pT_muons*TMath::Cos(phi_muons) + deepmet_pt_corr_central*TMath::Cos(deepmet_phi_corr_central) )") \
    # .Define("u_DeepMETCorr_y", "-(pT_muons*TMath::Sin(phi_muons) + deepmet_pt_corr_central*TMath::Sin(deepmet_phi_corr_central) )") \
    # .Define("u_DeepMETCorr_pt", "TMath::Sqrt(u_DeepMETCorr_x * u_DeepMETCorr_x + u_DeepMETCorr_y * u_DeepMETCorr_y)")

    # uncertainties
    for unc in uncs:
        rdf = rdf.Define(f"PuppiMET_pt{unc}UpFix", f"TMath::IsNaN(PuppiMET_pt{unc}Up) ? 0 : PuppiMET_pt{unc}Up") \
            .Define(f"PuppiMET_phi{unc}UpFix", f"TMath::IsNaN(PuppiMET_phi{unc}Up) ? 0 : PuppiMET_phi{unc}Up") \
            .Define(f"u_PUPPI{unc}_x", f"-(pT_muons*TMath::Cos(phi_muons) + PuppiMET_pt{unc}UpFix*TMath::Cos(PuppiMET_phi{unc}UpFix))") \
            .Define(f"u_PUPPI{unc}_y", f"-(pT_muons*TMath::Sin(phi_muons) + PuppiMET_pt{unc}UpFix*TMath::Sin(PuppiMET_phi{unc}UpFix))") \
            .Define(f"u_PUPPI{unc}_pt", f"TMath::Sqrt(u_PUPPI{unc}_x * u_PUPPI{unc}_x + u_PUPPI{unc}_y * u_PUPPI{unc}_y)")

    # for uncertainties
    # unc_string = "0."
    # for unc in ["JES", "JER", "Unclustered"]:
    #    rdf = rdf.Define(f"PuppiMET_pt{unc}UpFix", f"TMath::IsNaN(PuppiMET_pt{unc}Up) ? 0 : PuppiMET_pt{unc}Up") \
    #        .Define(f"PuppiMET_phi{unc}UpFix", f"TMath::IsNaN(PuppiMET_phi{unc}Up) ? 0 : PuppiMET_phi{unc}Up") \
    #        .Define(f"u_PUPPI{unc}_x", f"-(pT_muons*TMath::Cos(phi_muons) + PuppiMET_pt{unc}UpFix*TMath::Cos(PuppiMET_phi{unc}UpFix))") \
    #        .Define(f"u_PUPPI{unc}_y", f"-(pT_muons*TMath::Sin(phi_muons) + PuppiMET_pt{unc}UpFix*TMath::Sin(PuppiMET_phi{unc}UpFix))") \
    #        .Define(f"u_PUPPI{unc}_pt", f"TMath::Sqrt(u_PUPPI{unc}_x * u_PUPPI{unc}_x + u_PUPPI{unc}_y * u_PUPPI{unc}_y)")  \

    #    unc_string += f"+  TMath::Power(u_PUPPI_AXIS - u_PUPPI{unc}_AXIS, 2)"

    # quad sum of the uncertainties
    # rdf = rdf.Define("u_PUPPIUnc_x", f"u_PUPPI_x + TMath::Sqrt({unc_string.replace('AXIS', 'x')})") \
    #         .Define("u_PUPPIUnc_y", f"u_PUPPI_y + TMath::Sqrt({unc_string.replace('AXIS', 'y')})") \
    #         .Define("u_PUPPIUnc_pt", "TMath::Sqrt(u_PUPPIUnc_x * u_PUPPIUnc_x + u_PUPPIUnc_y * u_PUPPIUnc_y)")
    return rdf


rdf_data = prepareVars(rdf_data)
rdf_MC = prepareVars(rdf_MC)
rdf_bkg = prepareVars(rdf_bkg)

rdf_data = rdf_data.Define("weight_dummy", "1.0").Define(
    "weight_withBtag", "weight_corr * (jet_CSVLoose_n < 1) * metFilters")
rdf_MC = rdf_MC.Define("weight_dummy", "1.0").Define(
    "weight_withBtag", "weight_corr * (jet_CSVLoose_n < 1) * metFilters")
rdf_bkg = rdf_bkg.Define("weight_dummy", "1.0").Define(
    "weight_withBtag", "weight_corr * (jet_CSVLoose_n < 1) * metFilters")
# weightname = "weight_dummy"
weightname = "weight_withBtag"

# three bins of qT: inclusive, 0-50, 50-inf
rdf_data_qTLow = rdf_data.Filter("Z_pt < 50")
rdf_data_qTHigh = rdf_data.Filter("Z_pt >= 50")
rdf_MC_qTLow = rdf_MC.Filter("Z_pt < 50")
rdf_MC_qTHigh = rdf_MC.Filter("Z_pt >= 50")
rdf_bkg_qTLow = rdf_bkg.Filter("Z_pt < 50")
rdf_bkg_qTHigh = rdf_bkg.Filter("Z_pt >= 50")


# jet_n bins
rdf_data_njet0 = rdf_data.Filter("jet_n == 0")
rdf_data_njet1 = rdf_data.Filter("jet_n == 1")
rdf_data_njet2 = rdf_data.Filter("jet_n >= 2")
rdf_MC_njet0 = rdf_MC.Filter("jet_n == 0")
rdf_MC_njet1 = rdf_MC.Filter("jet_n == 1")
rdf_MC_njet2 = rdf_MC.Filter("jet_n >= 2")
rdf_bkg_njet0 = rdf_bkg.Filter("jet_n == 0")
rdf_bkg_njet1 = rdf_bkg.Filter("jet_n == 1")
rdf_bkg_njet2 = rdf_bkg.Filter("jet_n >= 2")


rdfs = [[rdf_data, rdf_MC, rdf_bkg]]
rdfs = [[rdf_data, rdf_MC, rdf_bkg],
        [rdf_data_qTLow, rdf_MC_qTLow, rdf_bkg_qTLow],
        [rdf_data_qTHigh, rdf_MC_qTHigh, rdf_bkg_qTHigh]]

suffixes = ["", "_qTLow", "_qTHigh"]
extraHeaders = {
    "": None,
    "_qTLow": "q_{T} < 50 GeV",
    "_qTHigh": "q_{T} > 50 GeV",
}

# rdfs = [[rdf_data, rdf_MC, rdf_bkg], [rdf_data_njet0, rdf_MC_njet0, rdf_bkg_njet0],
#        [rdf_data_njet1, rdf_MC_njet1, rdf_bkg_njet1], [rdf_data_njet2, rdf_MC_njet2, rdf_bkg_njet2]]
# suffixes = ["", "_njet0", "_njet1", "_njet2"]
# extraHeaders = {
#    "": None,
#    "_njet0": "n_{jet} = 0",
#    "_njet1": "n_{jet} = 1",
#    "_njet2": "n_{jet} >= 2",
# }


# recoils = ["PF", "PUPPI", "DeepMET", "DeepMETCorr"]
# recoils = ["PF", "PUPPI", "PUPPIUnc", "DeepMET"]
recoils = ["PF", "PUPPI", "DeepMET"]
# recoils = ["PF", "PFUp", "PFDown"]
# recoils = ["PF", "PUPPI", "DeepMET", "RawPF", "RawPUPPI"]

recoils_uncs = ["PUPPI" + unc for unc in uncs]

colors = {
    "PF": 1,
    "PUPPI": 2,
    "GEN": 6,
    "DeepMET": 4,
    "DeepMETCorr": 8,
    "PUPPIUnc": 6,
    "RawPF": 5,
    "RawPUPPI": 7,
    "PFUp": 9,
    "PFDown": 7
}

labels = {
    "TK": "TK",
    "PF": "PF",
    "PUPPI": "PUPPI",
    "GEN": "GEN",
    "TKPHO": "TK+Photon",
    "DeepMET": "DeepMET",
    "DeepMETCorr": "DeepMET",
    "PUPPIUnc": "PUPPIUnc",
    "RawPF": "RawPF",
    "RawPUPPI": "RawPUPPI",
    "PFUp": "PF JER Up",
    "PFDown": "PF JER Down",
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

xbins_qT = getpTBins()
xbins_qT_resp = getpTResponseBins()
xbins_nVtx = getnVtxBins()

hresponses_inclusive = None

hists_to_dump = []

# loop over the different qT bins
for idx, [rdf_data_tmp, rdf_MC_tmp, rdf_bkg_tmp] in enumerate(rdfs):
    suffix = suffixes[idx]
    extraHeader = extraHeaders[suffix]
    if len(suffix) > 0:
        suffix = "_" + suffix

    recoilanalyzer = RecoilAnalyzer(
        rdf_data_tmp, recoils + recoils_uncs, rdfMC=rdf_MC_tmp, rdfBkg=rdf_bkg_tmp, name="recoilanalyzer_" + suffix, useRMS=useRMS, weightname=weightname)
    recoilanalyzer.prepareVars()
    recoilanalyzer.prepareResponses('u_GEN_pt', xbins_qT)
    recoilanalyzer.prepareResolutions('u_GEN_pt', xbins_qT, 400, -200, 200)
    recoilanalyzer.prepareResponses(nPV, xbins_nVtx)
    recoilanalyzer.prepareResolutions(nPV, xbins_nVtx, 400, -200, 200)

    # just one bin, to calculate the response correction
    recoilanalyzer.prepareResponses('pT_muons', xbins_qT_resp)

    hresponses = recoilanalyzer.getResponses('u_GEN_pt')
    hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions(
        'u_GEN_pt')
    hresponses_nVtx = recoilanalyzer.getResponses(nPV)
    hresols_paral_diff_VS_nVtx, hresols_perp_VS_nVtx = recoilanalyzer.getResolutions(
        nPV)

    if hresponses_inclusive is None:
        hresponses_inclusive = recoilanalyzer.getResponses('pT_muons')

    values_responses = OrderedDict()
    values_responses_MC = OrderedDict()
    for itype in recoils + recoils_uncs:
        print("hresponses_inclusive in data for ", itype, " is ",
              hresponses_inclusive[itype].GetBinContent(1))
        print("hresponses_inclusive in MC for ", itype, " is ",
              hresponses_inclusive[itype + "_MC"].GetBinContent(1))
        values_responses[itype] = hresponses_inclusive[itype].GetBinContent(1)
        values_responses_MC[itype] = hresponses_inclusive[itype +
                                                          "_MC"].GetBinContent(1)

    if applySc:
        # if idx == 0:
        #    ROOT.gInterpreter.Declare(get_response_code)
        # create branch with the scale factors
        # for itype in recoils:
        #    # "dynamic scopes" to create a variable holding histograms
        #    # ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}{suffix}= {HNAME} ".format(
        #    #    RECOIL=itype, suffix=suffix, HNAME=hresponses[itype].GetName()))
        #    # ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}_MC{suffix}= {HNAME} ".format(
        #    #    RECOIL=itype, suffix=suffix, HNAME=hresponses[itype+"_MC"].GetName()))
        #    # rdf_data_tmp = rdf_data_tmp.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL}{suffix})'.format(RECOIL=itype, suffix=suffix)) \
        #    #    .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
        #    #    .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

        #    # rdf_MC_tmp = rdf_MC_tmp.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL}_MC{suffix})'.format(RECOIL=itype, suffix=suffix)) \
        #    #    .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
        #    #    .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

        #    rdf_data_tmp = rdf_data_tmp.Define("{RECOIL}_scale".format(RECOIL=itype), f"1.0 / {values_responses[itype]}") \
        #        .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
        #        .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))
        #
        #    # bkg: increase by the same amount of data, such that resol(data - bkg) is consistent

        #    rdf_MC_tmp = rdf_MC_tmp.Define("{RECOIL}_scale".format(RECOIL=itype), f"1.0 / {values_responses_MC[itype]}") \
        #        .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
        #        .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

        # recoilsSc = [itype + "Sc" for itype in recoils]
        # recoilanalyzerSc = RecoilAnalyzer(
        #    rdf_data_tmp, recoilsSc, rdfMC=rdf_MC_tmp, rdfBkg=rdf_bkg_tmp, name="recoilanalyzer_Scaled" + suffix, useRMS=useRMS, weightname=weightname)
        # recoilanalyzerSc.prepareVars()
        # recoilanalyzerSc.prepareResponses('u_GEN_pt', xbins_qT)
        # recoilanalyzerSc.prepareResolutions(
        #    'u_GEN_pt', xbins_qT, 400, -200, 200)
        # recoilanalyzerSc.prepareResponses(nPV, xbins_nVtx)
        # recoilanalyzerSc.prepareResolutions(nPV, xbins_nVtx, 400, -200, 200)

        # hresponsesSc = recoilanalyzerSc.getResponses('u_GEN_pt')
        # hresolsSc_paral_diff, hresolsSc_perp = recoilanalyzerSc.getResolutions(
        #    'u_GEN_pt')
        # hresponsesSc_nVtx = recoilanalyzerSc.getResponses(nPV)
        # hresolsSc_paral_diff_VS_nVtx, hresolsSc_perp_VS_nVtx = recoilanalyzerSc.getResolutions(
        #    nPV)

        # when using constant correction, it is equivalent to just scale the resolutions by constant
        hresolsSc_paral_diff = OrderedDict()
        hresolsSc_perp = OrderedDict()
        hresolsSc_paral_diff_VS_nVtx = OrderedDict()
        hresolsSc_perp_VS_nVtx = OrderedDict()

        for itype in recoils + recoils_uncs:
            resp = values_responses[itype]

            hresolsSc_paral_diff[itype] = hresols_paral_diff[itype].Clone(
                hresols_paral_diff[itype].GetName() + "Sc_paral_diff")
            hresolsSc_paral_diff[itype].Scale(1.0 / resp)

            hresolsSc_perp[itype] = hresols_perp[itype].Clone(
                hresols_perp[itype].GetName() + "Sc_perp")
            hresolsSc_perp[itype].Scale(1.0 / resp)

            hresolsSc_paral_diff_VS_nVtx[itype] = hresols_paral_diff_VS_nVtx[itype].Clone(
                hresols_paral_diff_VS_nVtx[itype].GetName() + "Sc_paral_diff_VS_nVtx")
            hresolsSc_paral_diff_VS_nVtx[itype].Scale(1.0 / resp)

            hresolsSc_perp_VS_nVtx[itype] = hresols_perp_VS_nVtx[itype].Clone(
                hresols_perp_VS_nVtx[itype].GetName() + "Sc_perp_VS_nVtx")
            hresolsSc_perp_VS_nVtx[itype].Scale(1.0 / resp)

        for itype in recoils + recoils_uncs:
            resp_MC = values_responses_MC[itype]
            hresolsSc_paral_diff[itype + "_MC"] = hresols_paral_diff[itype + "_MC"].Clone(
                hresols_paral_diff[itype + "_MC"].GetName() + "Sc_paral_diff_MC")
            hresolsSc_paral_diff[itype + "_MC"].Scale(1.0 / resp_MC)

            hresolsSc_perp[itype + "_MC"] = hresols_perp[itype + "_MC"].Clone(
                hresols_perp[itype + "_MC"].GetName() + "Sc_perp_MC")
            hresolsSc_perp[itype + "_MC"].Scale(1.0 / resp_MC)

            hresolsSc_paral_diff_VS_nVtx[itype + "_MC"] = hresols_paral_diff_VS_nVtx[itype + "_MC"].Clone(
                hresols_paral_diff_VS_nVtx[itype + "_MC"].GetName() + "Sc_paral_diff_VS_nVtx_MC")
            hresolsSc_paral_diff_VS_nVtx[itype + "_MC"].Scale(1.0 / resp_MC)

            hresolsSc_perp_VS_nVtx[itype + "_MC"] = hresols_perp_VS_nVtx[itype + "_MC"].Clone(
                hresols_perp_VS_nVtx[itype + "_MC"].GetName() + "Sc_perp_VS_nVtx_MC")
            hresolsSc_perp_VS_nVtx[itype + "_MC"].Scale(1.0 / resp_MC)

    qtmin, qtmax = getqTRange()
    qtlabel = getqTLabel()
    nvtxlabel = getnVtxLabel()
    uparallabel = getUparalLabel()
    uperplabel = getUperpLabel()
    responselabel = getResponseLabel()

    nvtxmin, nvtxmax = getVtxRange()

    def GetColors(hdict):
        hcolors = []
        for itype in hdict.keys():
            itype = itype.replace("_MC", "")
            hcolors.append(colors[itype])
        return hcolors

    def GetLegends(hdict):
        return [labels[itype] for itype in hdict.keys() if "_MC" not in itype]

    def GetLineStyles(hdict):
        linestyles = []
        for itype in hdict.keys():
            if "_MC" not in itype:
                linestyles.append(1)
            else:
                if "PUPPI" in itype:
                    linestyles.append(2)
                elif "PF" in itype:
                    linestyles.append(5)
                elif "DeepMET" in itype:
                    linestyles.append(8)
                else:
                    linestyles.append(3)
        return linestyles

    def GetMarkers(hdict):
        return [markers[itype] if "_MC" not in itype else 1 for itype in hdict.keys()]

    def GetDrawOptions(hdict):
        return ["EP" if "_MC" not in itype else "HIST" for itype in hdict.keys()]

    hresponses_todraw = OrderedDict()
    for k, v in hresponses.items():
        if k in recoils_uncs:
            continue
        if k.replace("_MC", "") in recoils_uncs:
            continue
        hresponses_todraw[k] = v

    hcolors = GetColors(hresponses_todraw)
    hmarkers = GetMarkers(hresponses_todraw)
    linestyles = GetLineStyles(hresponses_todraw)
    drawoptions = GetDrawOptions(hresponses_todraw)
    legends = GetLegends(hresponses_todraw)

    def GetRatios(hdict):
        hratios = {}
        for idx, itype in enumerate(hdict.keys()):
            if "_MC" not in itype and itype + "_MC" in hdict.keys():
                hratios[itype] = hdict[itype].Clone(itype + "_ratio")
                hratios[itype].Divide(hdict[itype + "_MC"])

                icolor = hcolors[idx]
                hratios[itype].SetLineColor(icolor)
                hratios[itype].SetMarkerColor(icolor)
                imarker = hmarkers[idx]
                hratios[itype].SetMarkerStyle(imarker)

        return hratios

    extraToDraw = ROOT.TPaveText(0.30, 0.15, 0.90, 0.25, "NDC")
    # extraToDraw.SetFillColor(0)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    if combineDataMC:
        extraToDraw.SetTextSize(0.055)
    else:
        extraToDraw.SetTextSize(0.04)

    n_todraw = int(len(hresponses_todraw) / 2)

    args = {
        "mycolors": hcolors,
        "markerstyles": hmarkers,
        "linestyles": linestyles,
        "drawoptions": drawoptions,
        "outdir": outdir,
        "noLumi": noLumi,
        "dology": False,
        "drawashist": False,
        "extraToDraw": extraToDraw,
        "legendoptions": ["LEP"] * n_todraw + ["L"] * n_todraw,
        "lheader": "MC",
        "extralheader": "Data",
        "extralabels": [""] * n_todraw,
        "yrmin": 0.76,
        "yrmax": 1.24,
        "yrlabel": "Data  / MC",
        "doPAS": doPAS,
    }
    args_DataORMC = {
        "mycolors": hcolors,
        "markerstyles": hmarkers,
        "linestyles": linestyles,
        "drawoptions": drawoptions,
        "outdir": outdir,
        "noLumi": noLumi,
        "dology": False,
        "drawashist": False,
        "extraToDraw": extraToDraw,
        "legendoptions": ["LEP"] * n_todraw + ["L"] * n_todraw,
        "doPAS": doPAS,
    }

    if extraHeader:
        extraToDraw.AddText(extraHeader)

    def DrawHistosDataMC(hdict, xmin, xmax, xlabel, ymin, ymax, ylabel, name, legendPos=[0.58, 0.15, 0.88, 0.40], DataOnly=False, MCOnly=False, inPaper=False, **kwargs):
        args_temp = args.copy()
        if DataOnly or MCOnly:
            args_temp = args_DataORMC.copy()
        args_temp.update(kwargs)
        args_temp['legendPos'] = legendPos

        if MCOnly:
            args_temp['drawashist'] = True
            args_temp['legendoptions'] = args_temp['legendoptions'][n_todraw:]
            args_temp['MCOnly'] = True

        hdict_todraw = OrderedDict()
        for k, v in hdict.items():
            # do not draw the uncs in the upper panel
            if k in recoils_uncs:
                continue
            if k.replace("_MC", "") in recoils_uncs:
                continue
            if DataOnly and "_MC" in k:
                continue
            if MCOnly and "_MC" not in k:
                continue
            hdict_todraw[k] = v

        if not DataOnly and not MCOnly:
            if "PUPPI" in hdict_todraw.keys():
                hratio_center = None
                if 'showratio' in args_temp.keys() and args_temp['showratio']:
                    # make unc for the ratio panel
                    val_center = "PUPPI"
                    hnum = hdict['PUPPI'].Clone(
                        f"{hdict['PUPPI'].GetName()}_numForUnc")
                    hratio_center = hnum.Clone(f"{hnum.GetName()}_center")
                    hratio_center.Divide(hdict[f"{val_center}_MC"])

                    hratio_uncs = []
                    for recoil_unc in recoils_uncs:
                        if recoil_unc in hdict.keys():
                            # data / MC_unc
                            # center value for data, variation for MC
                            hratio_unc = hnum.Clone(
                                f"{hnum.GetName()}_unc_for_{recoil_unc}")
                            hratio_unc.Divide(hdict[f"{recoil_unc}_MC"])
                            hratio_uncs.append(hratio_unc)

                    for ibin in range(1, hratio_center.GetNbinsX() + 1):
                        temp_center = hratio_center.GetBinContent(ibin)
                        temp = 0.
                        for hratio_unc in hratio_uncs:
                            temp += (hratio_unc.GetBinContent(ibin) -
                                     temp_center) ** 2
                        temp = math.sqrt(temp)
                        hratio_center.SetBinError(ibin, temp)
                        hratio_center.SetBinContent(ibin, 1.0)

                args_temp['hratiopanel'] = hratio_center

            hratios = GetRatios(hdict_todraw)
            args_temp['hratios'] = hratios

        DrawHistos(hdict_todraw.values(), legends, xmin, xmax,
                   xlabel, ymin, ymax, ylabel, name, inPaper=inPaper, **args_temp)
        if inPaper:
            global hists_to_dump
            hists_to_dump += hdict_todraw.values()

    DrawHistosDataMC(hresponses, qtmin, qtmax, qtlabel, 0., 1.29, responselabel,
                     "reco_recoil_response" + suffix, legendPos=[0.58, 0.20, 0.88, 0.40], inPaper=True)

    args['showratio'] = True
    DrawHistosDataMC(hresols_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
                     "reco_recoil_resol_paral" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
    DrawHistosDataMC(hresols_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
                     "reco_recoil_resol_perp" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
    DrawHistosDataMC(hresponses_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0., 1.15, responselabel,
                     "reco_recoil_response_VS_nVtx" + suffix, legendPos=[0.78, 0.20, 0.96, 0.40])
    DrawHistosDataMC(hresols_paral_diff_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0,
                     uparallabel, "reco_recoil_resol_paral_VS_nVtx" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
    DrawHistosDataMC(hresols_perp_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel,
                     "reco_recoil_resol_perp_VS_nVtx" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])

    if applySc:
        extraToDraw.Clear()
        todraw_string = "Response corrected"
        if extraHeader:
            todraw_string += f", {extraHeader}"
        extraToDraw.AddText(todraw_string)
        #
        # Scaled
        #
        if combineDataMC:
            # did not find a way to update the extraToDraw position...
            # reset the variable here
            extraToDraw = ROOT.TPaveText(0.30, 0.10, 0.90, 0.15, "NDC")
            # extraToDraw.SetFillColor(0)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.055)
            extraToDraw.AddText(todraw_string)
            args['extraToDraw'] = extraToDraw
            DrawHistosDataMC(hresponses, qtmin, qtmax, qtlabel, 0., 1.29, responselabel,
                             "reco_recoil_response_withratio" + suffix, legendPos=[0.58, 0.20, 0.88, 0.40], inPaper=True)
            # DrawHistosDataMC(hresolsSc_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
            #                 "reco_recoil_resol_paral_Scaled" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
            # DrawHistosDataMC(hresolsSc_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
            #                 "reco_recoil_resol_perp_Scaled" + suffix, legendPos=[0.33, 0.45, 0.58, 0.73])
            DrawHistosDataMC(hresolsSc_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
                             "reco_recoil_resol_paral_Scaled" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87], inPaper=True)
            DrawHistosDataMC(hresolsSc_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
                             "reco_recoil_resol_perp_Scaled" + suffix, legendPos=[0.33, 0.45, 0.58, 0.73], inPaper=True)
            DrawHistosDataMC(hresolsSc_paral_diff_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uparallabel,
                             "reco_recoil_resol_paral_VS_nVtx_Scaled" + suffix, legendPos=[0.33, 0.48, 0.58, 0.75], inPaper=True)
            DrawHistosDataMC(hresolsSc_perp_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel,
                             "reco_recoil_resol_perp_VS_nVtx_Scaled" + suffix, legendPos=[0.33, 0.48, 0.58, 0.75], inPaper=True)
        else:
            args['showratio'] = False
            suffix_data = suffix + "_Data"
            DrawHistosDataMC(hresolsSc_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
                             "reco_recoil_resol_paral_Scaled" + suffix_data, legendPos=[0.20, 0.60, 0.45, 0.78], DataOnly=True, inPaper=True)
            DrawHistosDataMC(hresolsSc_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
                             "reco_recoil_resol_perp_Scaled" + suffix_data, legendPos=[0.20, 0.60, 0.45, 0.78], DataOnly=True, inPaper=True)
            DrawHistosDataMC(hresolsSc_paral_diff_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uparallabel,
                             "reco_recoil_resol_paral_VS_nVtx_Scaled" + suffix_data, legendPos=[0.20, 0.60, 0.45, 0.78], DataOnly=True, inPaper=True)
            DrawHistosDataMC(hresolsSc_perp_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel,
                             "reco_recoil_resol_perp_VS_nVtx_Scaled" + suffix_data, legendPos=[0.20, 0.60, 0.45, 0.78], DataOnly=True, inPaper=True)

            # suffix_MC = suffix + "_MC"
            # DrawHistosDataMC(hresolsSc_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
            #                      "reco_recoil_resol_paral_Scaled" + suffix_MC, legendPos=[0.33, 0.70, 0.58, 0.87], MCOnly=True)
            # DrawHistosDataMC(hresolsSc_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
            #                      "reco_recoil_resol_perp_Scaled" + suffix_MC, legendPos=[0.33, 0.70, 0.58, 0.87], MCOnly=True)
            # DrawHistosDataMC(hresolsSc_paral_diff_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uparallabel,
            #                      "reco_recoil_resol_paral_VS_nVtx_Scaled" + suffix_MC, legendPos=[0.33, 0.70, 0.58, 0.87], MCOnly=True)
            # DrawHistosDataMC(hresolsSc_perp_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel,
            #                      "reco_recoil_resol_perp_VS_nVtx_Scaled" + suffix_MC, legendPos=[0.33, 0.70, 0.58, 0.87], MCOnly=True)

    recoilanalyzer.saveHistos(f"root/output_{suffix}.root")

# dump hists
ofile = ROOT.TFile(f"root/recoil_resol_dataMC.root", "RECREATE")
ofile.cd()
for hist in hists_to_dump:
    hist.SetDirectory(ofile)
    hist.Write()
ofile.Close()
print("Histograms dumped to root/recoil_resol_dataMC.root")
