'''
compare the performances of PF, Puppi, Deep (w/o bias term)
'''
import ROOT
import sys
sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
from myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, prepVars, get_response_code
from utils.RecoilAnalyzer import RecoilAnalyzer

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/ANPlots/Functions_cc.so")

itree = ROOT.TChain("Events")
itree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET_DYJets/myNanoProdMc2016_NANO_1_Skim.root")
itree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET_DYJets/myNanoProdMc2016_NANO_2_Skim.root")
itree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET_DYJets/myNanoProdMc2016_NANO_3_Skim.root")
itree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET_DYJets/myNanoProdMc2016_NANO_4_Skim.root")
itree.Add("/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET_DYJets/myNanoProdMc2016_NANO_9_Skim.root")

itreef = ROOT.TChain("tree")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_1.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_2.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_3.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_4.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_9.root")

itree.AddFriend(itreef)

rdf = ROOT.ROOT.RDataFrame( itree )

rdf = rdf.Define("u_PUPPI_pt", "calUT_PtPhi(pT_muons,   phi_muons, PuppiMET_pt, PuppiMET_phi)") \
         .Define("u_PUPPI_x",  "-(pT_muons * TMath::Cos(phi_muons) + PuppiMET_pt*TMath::Cos(PuppiMET_phi))") \
         .Define("u_PUPPI_y",  "-(pT_muons * TMath::Sin(phi_muons) + PuppiMET_pt*TMath::Sin(PuppiMET_phi))") \
         .Define("u_PF_pt",    "calUT_PtPhi(pT_muons,   phi_muons, MET_pt,      MET_phi)") \
         .Define("u_PF_x",     "-(pT_muons * TMath::Cos(phi_muons) + MET_pt*TMath::Cos(MET_phi))") \
         .Define("u_PF_y",     "-(pT_muons * TMath::Sin(phi_muons) + MET_pt*TMath::Sin(MET_phi))") \
         .Define("u_GEN_x",    "-(pT_muons * TMath::Cos(phi_muons) + GenMET_pt * TMath::Cos(GenMET_phi))") \
         .Define("u_GEN_y",    "-(pT_muons * TMath::Sin(phi_muons) + GenMET_pt * TMath::Sin(GenMET_phi))") \
         .Define("u_GEN_pt",   "TMath::Sqrt(u_GEN_x * u_GEN_x + u_GEN_y * u_GEN_y)") \
         .Define("u_DEEP_x",  "u_pred_x") \
         .Define("u_DEEP_y",  "u_pred_y") \
         .Define("u_DEEP_pt",  "TMath::Sqrt(u_pred_x * u_pred_x + u_pred_y * u_pred_y)") \
         .Define("u_WoBias_x", "event_XFromPx") \
         .Define("u_WoBias_y", "event_YFromPy") \
         .Define("u_WoBias_pt", "TMath::Sqrt(u_WoBias_x * u_WoBias_x + u_WoBias_y * u_WoBias_y)")

recoils = ["DEEP", "WoBias"]

#for itype in recoils:
#    # prepare the paral, perp, response, diff variables
#    rdf = prepVars(rdf, "u_{RECOIL}".format(RECOIL=itype), "u_GEN")

colors = {
            "TK": 3,
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "TKPHO": 9,
            "DEEP": 4,
            "WoBias": 4,
         }

labels = {
            "TK": "TK",
            "PF": "Type-I PF",
            "PUPPI": "Type-I PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DEEP": "DeepMET",
            "WoBias": "DeepMET w/o bias term",
         }

xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins() 

recoilanalyzer = RecoilAnalyzer(rdf, recoils)
recoilanalyzer.prepareVars()
recoilanalyzer.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzer.prepareResponses(   'PV_npvsGood', xbins_nVtx)
recoilanalyzer.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

hresponses = recoilanalyzer.getResponses('u_GEN_pt')
hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions('u_GEN_pt')
hresponses_nVtx = recoilanalyzer.getResponses('PV_npvsGood')
hresols_paral_diff_VS_nVtx, hresols_perp_VS_nVtx = recoilanalyzer.getResolutions('PV_npvsGood')


ROOT.gInterpreter.Declare(get_response_code)
# create branch with the scale factors
for itype in recoils:
    #"dynamic scopes" to create a variable holding histograms
    ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}= {HNAME} ".format(RECOIL=itype, HNAME=hresponses[itype].GetName()))
    rdf = rdf.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL})'.format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

#for itype in recoils:
#    # prepare the paral, perp, response, diff variables
#    rdf = prepVars(rdf, "u_{RECOIL}Sc".format(RECOIL=itype), "u_GEN")

recoilsSc = [itype + "Sc" for itype in recoils]
recoilanalyzerSc = RecoilAnalyzer(rdf, recoilsSc)
recoilanalyzerSc.prepareVars()
recoilanalyzerSc.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzerSc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzerSc.prepareResponses(   'PV_npvsGood', xbins_nVtx)
recoilanalyzerSc.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

hresponsesSc = recoilanalyzerSc.getResponses('u_GEN_pt')
hresolsSc_paral_diff, hresolsSc_perp = recoilanalyzerSc.getResolutions('u_GEN_pt')
hresponsesSc_nVtx = recoilanalyzerSc.getResponses('PV_npvsGood')
hresolsSc_paral_diff_VS_nVtx, hresolsSc_perp_VS_nVtx = recoilanalyzerSc.getResolutions('PV_npvsGood')


DrawHistos(hresponses.values(), [labels[itype] for itype in hresponses.keys()], 0, 150., "q_{T} [GeV]", 0., 1.15, "Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response", drawashist=True, dology=False, legendPos=[0.60, 0.17, 0.80, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[1,2])

DrawHistos(hresols_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 150., "q_{T} [GeV]", 0, 39.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], linestyles=[1,2])

DrawHistos(hresols_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 150., "q_{T} [GeV]", 0, 32.0, "#sigma u_{#perp} [GeV]", "reco_recoil_resol_perp", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()], linestyles=[1,2])

DrawHistos(hresponses_nVtx.values(), [labels[itype] for itype in hresponses_nVtx.keys()], 0, 40., "# Vertices", 0., 1.15, "Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[1,2])

DrawHistos(hresols_paral_diff_VS_nVtx.values(), [labels[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], 0, 50., "# Vertices", 0, 40.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], linestyles=[1,2])

DrawHistos(hresols_perp_VS_nVtx.values(), [labels[itype] for itype in hresols_perp_VS_nVtx.keys()], 0, 50., "# Vertices", 0, 40.0, "#sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], linestyles=[1,2])

##
## Scaled 
##
DrawHistos(hresponsesSc.values(), [labels[itype] for itype in hresponses.keys()], 0, 150., "q_{T} [GeV]", 0., 1.15, "Scaled Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[1,2])

DrawHistos(hresolsSc_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 150., "q_{T} [GeV]", 0, 45.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_Scaled", drawashist=True, dology=False, legendPos=[0.65, 0.73, 0.83, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], linestyles=[1,2])

DrawHistos(hresolsSc_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 150., "q_{T} [GeV]", 0, 40.0, "Scaled #sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_Scaled", drawashist=True, dology=False, legendPos=[0.65, 0.73, 0.85, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()], linestyles=[1,2])

DrawHistos(hresponsesSc_nVtx.values(), [labels[itype] for itype in hresponses_nVtx.keys()], 0, 40., "# Vertices", 0., 1.15, "Scaled Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[1,2])

DrawHistos(hresolsSc_paral_diff_VS_nVtx.values(), [labels[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], 0, 30., "# Vertices", 0, 40.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], linestyles=[1,2])

DrawHistos(hresolsSc_perp_VS_nVtx.values(), [labels[itype] for itype in hresols_perp_VS_nVtx.keys()], 0, 30., "# Vertices", 0, 40.0, "Scaled #sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], linestyles=[1,2])


#f1 = ROOT.TFile("root_h/output.root", "RECREATE")
#for h2 in h2ds_perp_VS_qT.values():
#    h2.SetDirectory(f1)
#    h2.Write()
#f1.Close()

raw_input()
