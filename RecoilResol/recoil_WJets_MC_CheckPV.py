import ROOT
import sys
sys.path.append("../RecoilResol/CMSPLOTS")
from CMSPLOTS.myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code, prepRecoilVars, getqTRange, getqTLabel, getnVtxLabel, getResponseLabel, getUparalLabel, getUperpLabel, getVtxRange, getpTResponseBins
from utils.RecoilAnalyzer import RecoilAnalyzer
import argparse

noLumi = False
MCOnly = True

useRMS=False

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("Functions_cc.so")

parser = argparse.ArgumentParser()
parser.add_argument("--includeZ", action="store_true", help="Include Z events", default=False)
parser.add_argument("--applySc", action="store_true", help="Apply scale factors")
parser.set_defaults(applySc=True)
args = parser.parse_args()

includeZ = args.includeZ
applySc = args.applySc
outdir = f"plots/MC/WJets"

ifiles = "/eos/cms/store/group/cmst3/group/wmass/w-mass-13TeV/NanoAOD/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_TrackFitV722_NanoProdv6/240509_052242/0000/NanoV9MCPostVFP_2*.root"
ifiles = "/home/yongbinfeng/Desktop/DeepMET/data/wjetsMC/NanoV9MCPostVFP_294.root"
rdf_org = ROOT.ROOT.RDataFrame("Events", ifiles)
rdf_org1 = rdf_org.Filter("nMuon >= 1")
rdf_org1 = rdf_org1.Define("Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]")
rdf_org2 = rdf_org1.Filter("Muon_pass0")
rdf = rdf_org2

rdf = rdf.Define("weight", "Generator_weight > 0 ? 1 : -1")

rdf = rdf.Define("muon_pt", "Muon_pt[0]").Define("muon_phi", "Muon_phi[0]")

rdf = rdf.Define("prefsrLeps", "prefsrLeptons(GenPart_status, GenPart_statusFlags, GenPart_pdgId, GenPart_genPartIdxMother)") \
         .Define("genl", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[prefsrLeps[0]], GenPart_eta[prefsrLeps[0]], GenPart_phi[prefsrLeps[0]], GenPart_mass[prefsrLeps[0]])") \
         .Define("genlanti", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[prefsrLeps[1]], GenPart_eta[prefsrLeps[1]], GenPart_phi[prefsrLeps[1]], GenPart_mass[prefsrLeps[1]])") \
         .Define("genV", "ROOT::Math::PxPyPzEVector(genl)+ROOT::Math::PxPyPzEVector(genlanti)") \
         .Define("V_pt", "genV.Pt()") \
         .Define("V_phi", "genV.Phi()")

#rdf = rdf.Define("V_pt", "TMath::Sqrt(Muon_pt[0] * Muon_pt[0] + GenMET_pt * GenMET_pt + 2 * Muon_pt[0] * GenMET_pt * TMath::Cos(Muon_phi[0] - GenMET_phi))")
#rdf = rdf.Define("V_phi", "TMath::ATan2(Muon_pt[0] * TMath::Sin(Muon_phi[0]) + GenMET_pt * TMath::Sin(GenMET_phi), Muon_pt[0] * TMath::Cos(Muon_phi[0]) + GenMET_pt * TMath::Cos(GenMET_phi))")

rdf = rdf.Define("PFMET_pt", "MET_pt") \
         .Define("PFMET_phi", "MET_phi") \
         .Define("PUPPIMET_pt", "PuppiMET_pt") \
         .Define("PUPPIMET_phi", "PuppiMET_phi") \
         .Define("DeepMETMET_pt", "DeepMETResolutionTune_pt")

rdf= rdf.Define("mT_PF", "calMT_fromMET_PtPhi(muon_pt, muon_phi, MET_pt, MET_phi)")
rdf= rdf.Define("mT_PUPPI", "calMT_fromMET_PtPhi(muon_pt, muon_phi, PuppiMET_pt, PuppiMET_phi)")
rdf= rdf.Define("mT_DeepMET", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")
rdf= rdf.Define("mT_DeepMETPVRobust", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETPVRobust_pt, DeepMETPVRobust_phi)")
rdf= rdf.Define("mT_DeepMETPVRobustNoPUPPI", "calMT_fromMET_PtPhi(muon_pt, muon_phi, DeepMETPVRobustNoPUPPI_pt, DeepMETPVRobustNoPUPPI_phi)")

rdf = rdf.Define("PF_pt", "MET_pt") \
         .Define("PF_phi", "MET_phi") \
         .Define("PUPPI_pt", "PuppiMET_pt") \
         .Define("PUPPI_phi", "PuppiMET_phi")

#recoils = ["PF", "PUPPI", "DeepMETResolutionTune", "DeepMETPVRobust", "DeepMETPVRobustNoPUPPI"]
recoils = ["PF", "PUPPI", "DeepMETResolutionTune", "DeepMETPVRobust"]

rdf = prepRecoilVars(rdf, "muon", recoils)
rdf = rdf.Define("u_GEN_pt", "V_pt").Define("u_GEN_x", "- V_pt * TMath::Cos(V_phi)").Define("u_GEN_y", "- V_pt * TMath::Sin(V_phi)")
h_pvIndex = rdf.Histo1D(("pv_Index", "pv_Index", 11, -1.5, 10.5), "PVRobustIndex")
h_pvIndices = [h_pvIndex]
h_pvIndices_labels = ["W"]

if includeZ:
    # include Z ntuples and plot also the PV index from Z's
    rdf_Z_org = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/group/cmst3/group/wmass/yofeng/NanoAOD/DYJetsToMuMu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoV9MCPostVFP_testPVRobustDM/240226_075056/0000/NanoV9MCPostVFP_1*.root")
    rdf_Z_org1 = rdf_Z_org.Filter("nMuon > 1")
    rdf_Z_org1 = rdf_Z_org1.Define("Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]")
    rdf_Z_org1 = rdf_Z_org1.Define("Muon_pass1", "Muon_pt[1] > 25.0 && abs(Muon_eta[1]) < 2.4 && Muon_pfRelIso04_all[1] < 0.15 && Muon_looseId[1]")

    rdf_Z_org2 = rdf_Z_org1.Filter("Muon_pass0 && Muon_pass1")
    rdf_Z = rdf_Z_org2

    rdf_Z = rdf_Z.Define("V_pt", "TMath::Sqrt(Muon_pt[0] * Muon_pt[0] + Muon_pt[1] * Muon_pt[1] + 2 * Muon_pt[0] * Muon_pt[1] * TMath::Cos(Muon_phi[0] - Muon_phi[1]))")
    rdf_Z = rdf_Z.Define("V_phi", "TMath::ATan2(Muon_pt[0] * TMath::Sin(Muon_phi[0]) + Muon_pt[1] * TMath::Sin(Muon_phi[1]), Muon_pt[0] * TMath::Cos(Muon_phi[0]) + Muon_pt[1] * TMath::Cos(Muon_phi[1]))")

    rdf_Z = rdf_Z.Define("PF_pt", "MET_pt") \
                 .Define("PF_phi", "MET_phi") \
                 .Define("PUPPI_pt", "PuppiMET_pt") \
                 .Define("PUPPI_phi", "PuppiMET_phi")

    rdf_Z = prepRecoilVars(rdf_Z, "V", recoils)
    rdf_Z = rdf_Z.Define("u_GEN_pt", "V_pt").Define("u_GEN_x", "- V_pt * TMath::Cos(V_phi)").Define("u_GEN_y", "- V_pt * TMath::Sin(V_phi)")

    h_pvIndex_Z = rdf_Z.Histo1D(("pv_Index", "pv_Index", 11, -1.5, 10.5), "PVRobustIndex")
    
    h_pvIndices.append(h_pvIndex_Z)
    h_pvIndices_labels.append("Z")

h_MTs = OrderedDict()
h_METs = OrderedDict()
for met in ["PF", "PUPPI", "DeepMET"]:
    hname = "h_mT_{}".format(met)
    h_MTs[met] = rdf.Histo1D((hname, hname, 40, 0, 120), "mT_{}".format(met))
    hname = "h_MET_{}".format(met)
    h_METs[met] = rdf.Histo1D((hname, hname, 40, 0, 100), "{}MET_pt".format(met))
    
# get resolutions
xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins() 
xbins_qT_resp = getpTResponseBins()

rdf = rdf.Define("GoodPVEvent", "PVRobustIndex == 0").Define("BadPVEvent", "PVRobustIndex != 0")
rdf_goodpv = rdf.Filter("GoodPVEvent")
rdf_badpv = rdf.Filter("BadPVEvent")

recoilanalyzer_goodpv = RecoilAnalyzer(rdf_goodpv, recoils, useRMS=useRMS, name="recoilanalyzer_goodpv")
recoilanalyzer_goodpv.prepareVars()
recoilanalyzer_goodpv.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer_goodpv.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzer_goodpv.prepareResponses(   'PV_npvsGood', xbins_nVtx)
recoilanalyzer_goodpv.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

hresponses_goodpv = recoilanalyzer_goodpv.getResponses('u_GEN_pt')
hresols_paral_diff_goodpv, hresols_perp_goodpv = recoilanalyzer_goodpv.getResolutions('u_GEN_pt')
hresponses_nVtx_goodpv = recoilanalyzer_goodpv.getResponses('PV_npvsGood')
hresols_paral_diff_VS_nVtx_goodpv, hresols_perp_VS_nVtx_goodpv = recoilanalyzer_goodpv.getResolutions('PV_npvsGood')

# just one bin, to calculate the response correction
recoilanalyzer_goodpv.prepareResponses('V_pt', xbins_qT_resp)
hresponses_inclusive_goodpv = recoilanalyzer_goodpv.getResponses('V_pt')

values_responses_goodpv = OrderedDict()
for itype in recoils:
    print("hresponses_inclusive in data GoodPV for ", itype, " is ",
          hresponses_inclusive_goodpv[itype].GetBinContent(1))
    values_responses_goodpv[itype] = hresponses_inclusive_goodpv[itype].GetBinContent(1)

recoilanalyzer_badpv = RecoilAnalyzer(rdf_badpv, recoils, useRMS=useRMS, name="recoilanalyzer_badpv")
recoilanalyzer_badpv.prepareVars()
recoilanalyzer_badpv.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer_badpv.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzer_badpv.prepareResponses(   'PV_npvsGood', xbins_nVtx)
recoilanalyzer_badpv.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

hresponses_badpv = recoilanalyzer_badpv.getResponses('u_GEN_pt')
hresols_paral_diff_badpv, hresols_perp_badpv = recoilanalyzer_badpv.getResolutions('u_GEN_pt')
hresponses_nVtx_badpv = recoilanalyzer_badpv.getResponses('PV_npvsGood')
hresols_paral_diff_VS_nVtx_badpv, hresols_perp_VS_nVtx_badpv = recoilanalyzer_badpv.getResolutions('PV_npvsGood')

# just one bin, to calculate the response correction
recoilanalyzer_badpv.prepareResponses('V_pt', xbins_qT_resp)
hresponses_inclusive_badpv = recoilanalyzer_badpv.getResponses('V_pt')

values_responses_badpv = OrderedDict()
for itype in recoils:
    print("hresponses_inclusive in data BadPV for ", itype, " is ",
          hresponses_inclusive_badpv[itype].GetBinContent(1))
    values_responses_badpv[itype] = hresponses_inclusive_badpv[itype].GetBinContent(1)


if applySc:
   #ROOT.gInterpreter.Declare(get_response_code)
   ## create branch with the scale factors
   #for itype in recoils:
   #    #"dynamic scopes" to create a variable holding histograms
   #    ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}_goodpv = {HNAME} ".format(RECOIL=itype, HNAME=hresponses_goodpv[itype].GetName()))
   #    
   #    rdf_goodpv = rdf_goodpv.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL}_goodpv)'.format(RECOIL=itype)) \
   #               .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
   #               .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))
   #               
   #    ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}_badpv = {HNAME} ".format(RECOIL=itype, HNAME=hresponses_badpv[itype].GetName()))
   #     
   #    rdf_badpv = rdf_badpv.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL}_badpv)'.format(RECOIL=itype)) \
   #                .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
   #                .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))


   #recoilsSc = [itype + "Sc" for itype in recoils]
   
   #recoilanalyzer_goodpv_Sc = RecoilAnalyzer(rdf_goodpv, recoilsSc, name = "recoilanalyzer_Scaled_goodpv", useRMS=useRMS)
   #recoilanalyzer_goodpv_Sc.prepareVars()
   #recoilanalyzer_goodpv_Sc.prepareResponses(   'u_GEN_pt', xbins_qT)
   #recoilanalyzer_goodpv_Sc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
   #recoilanalyzer_goodpv_Sc.prepareResponses(   'PV_npvsGood', xbins_nVtx)
   #recoilanalyzer_goodpv_Sc.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

   #hresponses_goodpv_Sc = recoilanalyzer_goodpv_Sc.getResponses('u_GEN_pt')
   #hresols_goodpv_Sc_paral_diff, hresols_goodpv_Sc_perp = recoilanalyzer_goodpv_Sc.getResolutions('u_GEN_pt')
   #hresponses_goodpv_Sc_nVtx = recoilanalyzer_goodpv_Sc.getResponses('PV_npvsGood')
   #hresols_goodpv_Sc_paral_diff_VS_nVtx, hresols_goodpv_Sc_perp_VS_nVtx = recoilanalyzer_goodpv_Sc.getResolutions('PV_npvsGood')
   
   #recoilanalyzer_badpv_Sc = RecoilAnalyzer(rdf_badpv, recoilsSc, name = "recoilanalyzer_Scaled_badpv", useRMS=useRMS)
   #recoilanalyzer_badpv_Sc.prepareVars()
   #recoilanalyzer_badpv_Sc.prepareResponses(   'u_GEN_pt', xbins_qT)
   #recoilanalyzer_badpv_Sc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
   #recoilanalyzer_badpv_Sc.prepareResponses(   'PV_npvsGood', xbins_nVtx)
   #recoilanalyzer_badpv_Sc.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

   #hresponses_badpv_Sc = recoilanalyzer_badpv_Sc.getResponses('u_GEN_pt')
   #hresols_badpv_Sc_paral_diff, hresols_badpv_Sc_perp = recoilanalyzer_badpv_Sc.getResolutions('u_GEN_pt')
   #hresponses_badpv_Sc_nVtx = recoilanalyzer_badpv_Sc.getResponses('PV_npvsGood')
   #hresols_badpv_Sc_paral_diff_VS_nVtx, hresols_badpv_Sc_perp_VS_nVtx = recoilanalyzer_badpv_Sc.getResolutions('PV_npvsGood')
   
   hresolsSc_paral_diff_goodpv = OrderedDict()
   hresolsSc_perp_goodpv = OrderedDict()
   hresolsSc_paral_diff_VS_nVtx_goodpv = OrderedDict()
   hresolsSc_perp_VS_nVtx_goodpv = OrderedDict()

   for itype in recoils:
       resp = values_responses_goodpv[itype] + 1e-6

       hresolsSc_paral_diff_goodpv[itype] = hresols_paral_diff_goodpv[itype].Clone(
           itype + "Sc_paral_diff")
       hresolsSc_paral_diff_goodpv[itype].Scale(1.0 / resp)

       hresolsSc_perp_goodpv[itype] = hresols_perp_goodpv[itype].Clone(
           itype + "Sc_perp")
       hresolsSc_perp_goodpv[itype].Scale(1.0 / resp)

       hresolsSc_paral_diff_VS_nVtx_goodpv[itype] = hresols_paral_diff_VS_nVtx_goodpv[itype].Clone(
           itype + "Sc_paral_diff_VS_nVtx")
       hresolsSc_paral_diff_VS_nVtx_goodpv[itype].Scale(1.0 / resp)

       hresolsSc_perp_VS_nVtx_goodpv[itype] = hresols_perp_VS_nVtx_goodpv[itype].Clone(
           itype + "Sc_perp_VS_nVtx")
       hresolsSc_perp_VS_nVtx_goodpv[itype].Scale(1.0 / resp)
       

   hresolsSc_paral_diff_badpv = OrderedDict()
   hresolsSc_perp_badpv = OrderedDict()
   hresolsSc_paral_diff_VS_nVtx_badpv = OrderedDict()
   hresolsSc_perp_VS_nVtx_badpv = OrderedDict()

   for itype in recoils:
       resp = values_responses_badpv[itype] + 1e-6

       hresolsSc_paral_diff_badpv[itype] = hresols_paral_diff_badpv[itype].Clone(
           itype + "Sc_paral_diff")
       hresolsSc_paral_diff_badpv[itype].Scale(1.0 / resp)

       hresolsSc_perp_badpv[itype] = hresols_perp_badpv[itype].Clone(
           itype + "Sc_perp")
       hresolsSc_perp_badpv[itype].Scale(1.0 / resp)

       hresolsSc_paral_diff_VS_nVtx_badpv[itype] = hresols_paral_diff_VS_nVtx_badpv[itype].Clone(
           itype + "Sc_paral_diff_VS_nVtx")
       hresolsSc_paral_diff_VS_nVtx_badpv[itype].Scale(1.0 / resp)

       hresolsSc_perp_VS_nVtx_badpv[itype] = hresols_perp_VS_nVtx_badpv[itype].Clone(
           itype + "Sc_perp_VS_nVtx")
       hresolsSc_perp_VS_nVtx_badpv[itype].Scale(1.0 / resp)
       
       
   
   

h_muon_pts = OrderedDict()
h_muon_etas = OrderedDict()
h_W_qts = OrderedDict()
for pvs in ["GoodPV", "BadPV"]:
    if pvs == "GoodPV":
        tmp = rdf_goodpv
    else:
        tmp = rdf_badpv
    hname = "h_pt_muons_{}".format(pvs)
    h_muon_pts[pvs] = tmp.Histo1D((hname, hname, 40, 0, 100), "muon_pt")
    hname = "h_eta_muons_{}".format(met)
    h_muon_etas[pvs] = tmp.Histo1D((hname, hname, 40, -2.5, 2.5), "Muon_eta")
    hname = "h_qt_W_{}".format(pvs)
    h_W_qts[pvs] = tmp.Histo1D((hname, hname, xbins_qT.size - 1, xbins_qT), "V_pt")

colors = {
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "DeepMET": 4,
            "DeepMETResolutionTune": 4,
            "DeepMETPVRobust": 6,
            "DeepMETPVRobustNoPUPPI": 7
         }

for itype in list(colors.keys()):
    colors[itype + "Sc"] = colors[itype]

labels = {
            "TK": "TK",
            "PF": "PF",
            "PUPPI": "PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DeepMET": "DeepMET",
            "DeepMETResolutionTune": "DeepMET",
            "DeepMETPVRobust": "DeepMET PV-Agnostic",
            "DeepMETPVRobustNoPUPPI": "DeepMET PVRobust NoPUPPI"
         }

for itype in list(labels.keys()):
    labels[itype + "Sc"] = labels[itype]
    

DrawHistos(h_pvIndices, h_pvIndices_labels, -1.5, 10.5, "Robust PV Index", 1e-5, 10.0, "A.U.", "W_pv_Index", donormalize=True, outdir=outdir, noLumi=noLumi, mycolors=[1, 2], MCOnly=MCOnly)

# make some kinematic distributions
DrawHistos(h_MTs.values(), [labels[itype] for itype in h_MTs.keys()], 0, 120, "m_{T} [GeV]", 0., 0.1, "A.U.", "MT_comp", drawashist=True, dology=False, legendPos=[0.20, 0.67, 0.38, 0.86], mycolors=[colors[itype] for itype in h_MTs.keys()], outdir=outdir, donormalize=True, MCOnly=MCOnly)

DrawHistos(h_METs.values(), [labels[itype] for itype in h_METs.keys()], 0, 100, "p^{miss}_{T} [GeV]", 0., 0.13, "A.U.", "MET_comp", drawashist=True, dology=False, legendPos=[0.62, 0.67, 0.80, 0.86], mycolors=[colors[itype] for itype in h_METs.keys()], outdir=outdir, donormalize=True, MCOnly=MCOnly)

DrawHistos(h_muon_pts.values(), ["PV Good", "PV Bad"], 0, 100, "#mu p_{T} [GeV]", 0., 0.15, "A.U.", "muon_pt_comp", drawashist=True, dology=False, legendPos=[0.65, 0.67, 0.90, 0.86], mycolors=[1]*2, linestyles=[1,2], outdir=outdir, donormalize=True, MCOnly=MCOnly)

DrawHistos(h_muon_etas.values(), ["PV Good", "PV Bad"], -2.5, 2.5, "#mu #eta", 0., 0.05, "A.U.", "muon_eta_comp", drawashist=True, dology=False, legendPos=[0.20, 0.67, 0.38, 0.86], mycolors=[1]*2, linestyles=[1,2], outdir=outdir, donormalize=True, MCOnly=MCOnly)

DrawHistos(h_W_qts.values(), ["PV Good", "PV Bad"], 0, 150, "W q_{T} [GeV]", 0., 0.3, "A.U.", "W_qt_comp", drawashist=True, dology=False, legendPos=[0.65, 0.67, 0.90, 0.86], mycolors=[1]*2, linestyles=[1,2], outdir=outdir, donormalize=True, MCOnly=MCOnly)


# resolution distributions

qtmin, qtmax = getqTRange()
qtlabel = getqTLabel()
nvtxlabel = getnVtxLabel()
uparallabel = getUparalLabel()
uperplabel = getUperpLabel()
responselabel = getResponseLabel()

nvtxmin, nvtxmax = getVtxRange()

args = {
      "mycolors": [colors[itype] for itype in hresponses_goodpv.keys()]*2,
      "linestyles": [1] * len(hresponses_goodpv) + [2] * len(hresponses_badpv),
      "outdir": outdir,
      "noLumi": noLumi,
      "dology": False,
      "drawashist": True,
      "legendoptions": ["L"] * len(hresponses_goodpv) + ["L"] * len(hresponses_badpv),
      #"lheader": "PVIncorrect",
      #"extralheader": "PVCorrect",
      "extralabels": [""]*len(hresponses_badpv),
      "MCOnly": MCOnly,
      "exlegoffset": 0.06
}


extraToDraw = ROOT.TPaveText(0.50, 0.37, 0.90, 0.42, "NDC")
extraToDraw.SetTextAlign(11)
extraToDraw.SetFillColorAlpha(0, 0)
extraToDraw.SetBorderSize(0)
extraToDraw.SetTextFont(42)
extraToDraw.SetTextSize(0.04)
extraToDraw.AddText("PV-Incorrect")

extraToDraw2 = ROOT.TPaveText(0.33, 0.37, 0.50, 0.42, "NDC")
extraToDraw2.SetFillColorAlpha(0, 0)
extraToDraw2.SetBorderSize(0)
extraToDraw2.SetTextAlign(31)
extraToDraw2.SetTextFont(42)
extraToDraw2.SetTextSize(0.04)
extraToDraw2.AddText("PV-Correct")

args["extraToDraw"] = [extraToDraw, extraToDraw2]

linestyles = [1] * len(hresponses_goodpv) + [2] * len(hresponses_badpv)

DrawHistos(list(hresponses_goodpv.values()) + list(hresponses_badpv.values()), [labels[itype] for itype in hresponses_goodpv.keys()], 0, qtmax, qtlabel, 0., 1.15, responselabel, "reco_recoil_response", legendPos=[0.50, 0.17, 0.80, 0.37], **args)

DrawHistos(list(hresols_paral_diff_goodpv.values()) + list(hresols_paral_diff_badpv.values()), [labels[itype] for itype in hresols_paral_diff_goodpv.keys()], 0, qtmax, qtlabel, 0, 39.0, uparallabel, "reco_recoil_resol_paral", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_goodpv.keys()]*2, noLumi=noLumi, outdir=outdir, linestyles=linestyles, MCOnly=MCOnly)

DrawHistos(list(hresols_perp_goodpv.values()) + list(hresols_perp_badpv.values()), [labels[itype] for itype in hresols_perp_goodpv.keys()], 0, qtmax, qtlabel, 0, 32.0, uperplabel, "reco_recoil_resol_perp", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_goodpv.keys()]*2, noLumi=noLumi, outdir=outdir, linestyles=linestyles, MCOnly=MCOnly)

DrawHistos(list(hresponses_nVtx_goodpv.values()) + list(hresponses_nVtx_badpv.values()), [labels[itype] for itype in hresponses_nVtx_goodpv.keys()], nvtxmin, nvtxmax, nvtxlabel, 0., 1.15, responselabel, "reco_recoil_response_VS_nVtx", drawashist=True, dology=False, legendPos=[0.45, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses_goodpv.keys()]*2, noLumi=noLumi, outdir=outdir, linestyles=linestyles, MCOnly=MCOnly)

DrawHistos(list(hresols_paral_diff_VS_nVtx_goodpv.values()) + list(hresols_paral_diff_VS_nVtx_badpv.values()), [labels[itype] for itype in hresols_paral_diff_VS_nVtx_goodpv.keys()], nvtxmin, nvtxmax, nvtxlabel, 0, 60.0, uparallabel, "reco_recoil_resol_paral_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx_goodpv.keys()]*2, noLumi=noLumi, outdir=outdir, linestyles=linestyles, MCOnly=MCOnly)

DrawHistos(list(hresols_perp_VS_nVtx_goodpv.values()) + list(hresols_perp_VS_nVtx_badpv.values()), [labels[itype] for itype in hresols_perp_VS_nVtx_goodpv.keys()], nvtxmin, nvtxmax, nvtxlabel, 0, 60.0, uperplabel, "reco_recoil_resol_perp_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx_goodpv.keys()]*2, noLumi=noLumi, outdir=outdir, linestyles=linestyles, MCOnly=MCOnly)


if applySc:
    extraToDraw = ROOT.TPaveText(0.60, 0.77, 0.90, 0.82, "NDC")
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText("Response corrected")
    
    extraToDraw2 = ROOT.TPaveText(0.37, 0.85, 0.67, 0.90, "NDC")
    extraToDraw2.SetFillColorAlpha(0, 0)
    extraToDraw2.SetTextAlign(11)
    extraToDraw2.SetBorderSize(0)
    extraToDraw2.SetTextFont(42)
    extraToDraw2.SetTextSize(0.04)
    extraToDraw2.AddText("PV-Incorrect")
    
    extraToDraw3 = ROOT.TPaveText(0.20, 0.85, 0.37, 0.90, "NDC")
    extraToDraw3.SetFillColorAlpha(0, 0)
    extraToDraw3.SetTextAlign(31)
    extraToDraw3.SetBorderSize(0)
    extraToDraw3.SetTextFont(42)
    extraToDraw3.SetTextSize(0.04)
    extraToDraw3.AddText("PV-Correct")
    
    args["extraToDraw"] = [extraToDraw, extraToDraw2, extraToDraw3]
    
    DrawHistos(list(hresolsSc_paral_diff_goodpv.values()) + list(hresolsSc_paral_diff_badpv.values()), [labels[itype] for itype in hresolsSc_paral_diff_goodpv.keys()], 0, qtmax, qtlabel, 0, 39.0, uparallabel, "reco_recoil_resol_paral_Scaled", legendPos=[0.30, 0.69, 0.65, 0.92], **args)
    
    DrawHistos(list(hresolsSc_perp_goodpv.values()) + list(hresolsSc_perp_badpv.values()), [labels[itype] for itype in hresolsSc_perp_goodpv.keys()], 0, qtmax, qtlabel, 0, 32.0, uperplabel, "reco_recoil_resol_perp_Scaled", legendPos=[0.30, 0.69, 0.70, 0.92], **args)
    
    DrawHistos(list(hresolsSc_paral_diff_VS_nVtx_goodpv.values()) + list(hresolsSc_paral_diff_VS_nVtx_badpv.values()), [labels[itype] for itype in hresolsSc_paral_diff_VS_nVtx_goodpv.keys()], nvtxmin, nvtxmax, nvtxlabel, 0, 79.0, uparallabel, "reco_recoil_resol_paral_Scaled_VS_nVtx", legendPos=[0.37, 0.65, 0.67, 0.85], **args)
    
    DrawHistos(list(hresolsSc_perp_VS_nVtx_goodpv.values()) + list(hresolsSc_perp_VS_nVtx_badpv.values()), [labels[itype] for itype in hresolsSc_perp_VS_nVtx_goodpv.keys()], nvtxmin, nvtxmax, nvtxlabel, 0, 79.0, uperplabel, "reco_recoil_resol_perp_Scaled_VS_nVtx", legendPos=[0.37, 0.65, 0.67, 0.85], **args)