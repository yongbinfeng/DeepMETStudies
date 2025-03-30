import sys
sys.path.append("../RecoilResol/CMSPLOTS")  # noqa
import argparse
from utils.RecoilAnalyzer import RecoilAnalyzer
from utils.utils import getpTBins, getpTResponseBins, getnVtxBins, get_response_code, getqTLabel, getnVtxLabel, getNPVString, getqTRange, getResponseLabel, getUparalLabel, getUperpLabel, getVtxRange
from collections import OrderedDict
from CMSPLOTS.myFunction import DrawHistos
import ROOT

noLumi = False
nPV = getNPVString()
useRMS = True

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

chainMC = ROOT.TChain("Events")
chainMC.Add(mcfile)
rdf_MC = ROOT.ROOT.RDataFrame(chainMC)

chainBkg = ROOT.TChain("Events")
for bkgfile in bkgfiles:
    chainBkg.Add(bkgfile)
rdf_bkg = ROOT.ROOT.RDataFrame(chainBkg)


# weight_corr includes both the pileup and Zpt corrections
weightname = "weight_corr"

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
             .Define("u_PF_x",     "-(pT_muons*TMath::Cos(phi_muons) + MET_pt*TMath::Cos(MET_phi))") \
             .Define("u_PF_y",     "-(pT_muons*TMath::Sin(phi_muons) + MET_pt*TMath::Sin(MET_phi))") \
             .Define("u_PF_pt",    "TMath::Sqrt(u_PF_x * u_PF_x + u_PF_y * u_PF_y)") \
             .Define("u_DeepMET_x",  "-(pT_muons*TMath::Cos(phi_muons) + DeepMETResolutionTune_pt*TMath::Cos(DeepMETResolutionTune_phi))") \
             .Define("u_DeepMET_y",  "-(pT_muons*TMath::Sin(phi_muons) + DeepMETResolutionTune_pt*TMath::Sin(DeepMETResolutionTune_phi))") \
             .Define("u_DeepMET_pt", "TMath::Sqrt(u_DeepMET_x * u_DeepMET_x + u_DeepMET_y * u_DeepMET_y)")
    # .Define("u_DeepMETCorr_x", "-(pT_muons*TMath::Cos(phi_muons) + deepmet_pt_corr_central*TMath::Cos(deepmet_phi_corr_central) )") \
    # .Define("u_DeepMETCorr_y", "-(pT_muons*TMath::Sin(phi_muons) + deepmet_pt_corr_central*TMath::Sin(deepmet_phi_corr_central) )") \
    # .Define("u_DeepMETCorr_pt", "TMath::Sqrt(u_DeepMETCorr_x * u_DeepMETCorr_x + u_DeepMETCorr_y * u_DeepMETCorr_y)")

    # for uncertainties
    unc_string = "0."
    for unc in ["JES", "JER", "Unclustered"]:
        rdf = rdf.Define(f"PuppiMET_pt{unc}UpFix", f"TMath::IsNaN(PuppiMET_pt{unc}Up) ? 0 : PuppiMET_pt{unc}Up") \
            .Define(f"PuppiMET_phi{unc}UpFix", f"TMath::IsNaN(PuppiMET_phi{unc}Up) ? 0 : PuppiMET_phi{unc}Up") \
            .Define(f"u_PUPPI{unc}_x", f"-(pT_muons*TMath::Cos(phi_muons) + PuppiMET_pt{unc}UpFix*TMath::Cos(PuppiMET_phi{unc}UpFix))") \
            .Define(f"u_PUPPI{unc}_y", f"-(pT_muons*TMath::Sin(phi_muons) + PuppiMET_pt{unc}UpFix*TMath::Sin(PuppiMET_phi{unc}UpFix))") \
            .Define(f"u_PUPPI{unc}_pt", f"TMath::Sqrt(u_PUPPI{unc}_x * u_PUPPI{unc}_x + u_PUPPI{unc}_y * u_PUPPI{unc}_y)")  \

        unc_string += f"+  TMath::Power(u_PUPPI_AXIS - u_PUPPI{unc}_AXIS, 2)"

    # quad sum of the uncertainties
    rdf = rdf.Define("u_PUPPIUnc_x", f"u_PUPPI_x + TMath::Sqrt({unc_string.replace('AXIS', 'x')})") \
             .Define("u_PUPPIUnc_y", f"u_PUPPI_y + TMath::Sqrt({unc_string.replace('AXIS', 'y')})") \
             .Define("u_PUPPIUnc_pt", "TMath::Sqrt(u_PUPPIUnc_x * u_PUPPIUnc_x + u_PUPPIUnc_y * u_PUPPIUnc_y)")
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


# rdfs = [[rdf_data, rdf_MC, rdf_bkg]]
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
recoils = ["PF", "PUPPI", "PUPPIUnc", "DeepMET"]

colors = {
    "PF": 1,
    "PUPPI": 2,
    "GEN": 6,
    "DeepMET": 4,
    "DeepMETCorr": 8,
    "PUPPIUnc": 6,
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
}

markers = {
    "PF": 20,
    "PUPPI": 21,
    "GEN": 22,
    "DeepMET": 22,
    "DeepMETCorr": 23,
    "PUPPIUnc": 24,
}

xbins_qT = getpTBins()
xbins_qT_resp = getpTResponseBins()
xbins_nVtx = getnVtxBins()

hresponses_inclusive = None

# loop over the different qT bins
for idx, [rdf_data_tmp, rdf_MC_tmp, rdf_bkg_tmp] in enumerate(rdfs):
    suffix = suffixes[idx]
    extraHeader = extraHeaders[suffix]
    if len(suffix) > 0:
        suffix = "_" + suffix

    recoilanalyzer = RecoilAnalyzer(
        rdf_data_tmp, recoils, rdfMC=rdf_MC_tmp, rdfBkg=rdf_bkg_tmp, name="recoilanalyzer_" + suffix, useRMS=useRMS, weightname=weightname)
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
    for itype in recoils:
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

        for itype in recoils:
            resp_MC = values_responses_MC[itype]
            hresolsSc_paral_diff[itype + "_MC"] = hresols_paral_diff[itype + "_MC"].Clone(
                itype + "Sc_paral_diff_MC")
            hresolsSc_paral_diff[itype + "_MC"].Scale(1.0 / resp_MC)

            hresolsSc_perp[itype + "_MC"] = hresols_perp[itype + "_MC"].Clone(
                itype + "Sc_perp_MC")
            hresolsSc_perp[itype + "_MC"].Scale(1.0 / resp_MC)

            hresolsSc_paral_diff_VS_nVtx[itype + "_MC"] = hresols_paral_diff_VS_nVtx[itype + "_MC"].Clone(
                itype + "Sc_paral_diff_VS_nVtx_MC")
            hresolsSc_paral_diff_VS_nVtx[itype + "_MC"].Scale(1.0 / resp_MC)

            hresolsSc_perp_VS_nVtx[itype + "_MC"] = hresols_perp_VS_nVtx[itype + "_MC"].Clone(
                itype + "Sc_perp_VS_nVtx_MC")
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
        return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]

    def GetMarkers(hdict):
        return [markers[itype] if "_MC" not in itype else 1 for itype in hdict.keys()]

    def GetDrawOptions(hdict):
        return ["EP" if "_MC" not in itype else "HIST" for itype in hdict.keys()]

    hresponses_todraw = OrderedDict()
    for k, v in hresponses.items():
        if k != "PUPPIUnc" and k != "PUPPIUnc_MC":
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

    extraToDraw = ROOT.TPaveText(0.60, 0.70, 0.90, 0.87, "NDC")
    # extraToDraw.SetFillColor(0)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.05)

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
    }

    if extraHeader:
        extraToDraw.AddText(extraHeader)

    def DrawHistosWithUncband(hdict, xmin, xmax, xlabel, ymin, ymax, ylabel, name, legendPos=[0.58, 0.15, 0.88, 0.40], **kwargs):
        args_temp = args.copy()
        args_temp.update(kwargs)
        args_temp['legendPos'] = legendPos

        hdict_todraw = OrderedDict()
        for k, v in hdict.items():
            if "PUPPIUnc" in k:
                continue
            hdict_todraw[k] = v

        if 'showratio' in args_temp.keys() and args_temp['showratio']:
            # make unc for the ratio panel
            hratiopanel = None
            if "PUPPIUnc" in hdict.keys() or "PUPPIUncSc" in hdict.keys():
                val_center = "PUPPI" if "PUPPI" in hdict.keys() else "PUPPISc"
                val_unc = "PUPPIUnc" if "PUPPIUnc" in hdict.keys() else "PUPPIUncSc"

                h_center_ratio = hdict[val_center].Clone(
                    f"{hdict[val_center].GetName()}_ratio")
                h_center_ratio.Divide(hdict[f"{val_center}_MC"])
                # h_unc_ratio = hdict["PUPPIUnc"].Clone(f"{hdict['PUPPIUnc'].GetName()}_ratio")
                h_unc_ratio = hdict[val_center].Clone(
                    f"{hdict[val_center].GetName()}_ratio_unc")
                h_unc_ratio.Divide(hdict[f"{val_unc}_MC"])
                hratiopanel = h_center_ratio.Clone(
                    f"{h_center_ratio.GetName()}_ratioForUnc")
                for ibin in range(1, h_center_ratio.GetNbinsX() + 1):
                    hratiopanel.SetBinError(ibin, abs(h_unc_ratio.GetBinContent(
                        ibin) - h_center_ratio.GetBinContent(ibin)))
                    hratiopanel.SetBinContent(ibin, 1.0)

            hratios = GetRatios(hdict_todraw)

            args_temp['hratiopanel'] = hratiopanel
            args_temp['hratios'] = hratios

        DrawHistos(hdict_todraw.values(), legends, xmin, xmax,
                   xlabel, ymin, ymax, ylabel, name, **args_temp)

    DrawHistosWithUncband(hresponses, qtmin, qtmax, qtlabel, 0., 1.19, responselabel,
                          "reco_recoil_response" + suffix, legendPos=[0.58, 0.20, 0.88, 0.40])

    args['showratio'] = True
    DrawHistosWithUncband(hresols_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
                          "reco_recoil_resol_paral" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
    DrawHistosWithUncband(hresols_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
                          "reco_recoil_resol_perp" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
    DrawHistosWithUncband(hresponses_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0., 1.15, responselabel,
                          "reco_recoil_response_VS_nVtx" + suffix, legendPos=[0.78, 0.20, 0.96, 0.40])
    DrawHistosWithUncband(hresols_paral_diff_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0,
                          uparallabel, "reco_recoil_resol_paral_VS_nVtx" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
    DrawHistosWithUncband(hresols_perp_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel,
                          "reco_recoil_resol_perp_VS_nVtx" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])

    if applySc:
        extraToDraw.Clear()
        extraToDraw.AddText("Response corrected")
        if extraHeader:
            extraToDraw.AddText(extraHeader)
        #
        # Scaled
        #
        DrawHistosWithUncband(hresolsSc_paral_diff, qtmin, qtmax, qtlabel, 0, 50.0, uparallabel,
                              "reco_recoil_resol_paral_Scaled" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
        DrawHistosWithUncband(hresolsSc_perp, qtmin, qtmax, qtlabel, 0, 50.0, uperplabel,
                              "reco_recoil_resol_perp_Scaled" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
        DrawHistosWithUncband(hresolsSc_paral_diff_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uparallabel,
                              "reco_recoil_resol_paral_VS_nVtx_Scaled" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])
        DrawHistosWithUncband(hresolsSc_perp_VS_nVtx, nvtxmin, nvtxmax, nvtxlabel, 0, 50.0, uperplabel,
                              "reco_recoil_resol_perp_VS_nVtx_Scaled" + suffix, legendPos=[0.33, 0.60, 0.58, 0.87])

    recoilanalyzer.saveHistos(f"root/output_{suffix}.root")
