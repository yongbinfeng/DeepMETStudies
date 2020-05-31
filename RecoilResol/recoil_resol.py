import ROOT
import sys
sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
from myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code
from utils.RecoilAnalyzer import RecoilAnalyzer

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/ANPlots/Functions_cc.so")

rdf = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET/myNanoProdMc2016_NANO_[1-9]_Skim.root")
#rdf = rdf_org.Filter("pT_W<30.0")

rdf = rdf.Define("u_PUPPI_pt", "calUT_PtPhi(Muon_pt[0], Muon_phi[0], PuppiMET_pt, PuppiMET_phi)") \
         .Define("u_PUPPI_x",  "-(Muon_pt[0]*TMath::Cos(Muon_phi[0]) + PuppiMET_pt*TMath::Cos(PuppiMET_phi))") \
         .Define("u_PUPPI_y",  "-(Muon_pt[0]*TMath::Sin(Muon_phi[0]) + PuppiMET_pt*TMath::Sin(PuppiMET_phi))") \
         .Define("u_PF_pt",    "calUT_PtPhi(Muon_pt[0], Muon_phi[0], MET_pt,      MET_phi)") \
         .Define("u_PF_x",     "-(Muon_pt[0]*TMath::Cos(Muon_phi[0]) + MET_pt*TMath::Cos(MET_phi))") \
         .Define("u_PF_y",     "-(Muon_pt[0]*TMath::Sin(Muon_phi[0]) + MET_pt*TMath::Sin(MET_phi))") \
         .Define("u_GEN_x",    "-(pT_muons * TMath::Cos(phi_muons) + GenMET_pt * TMath::Cos(GenMET_phi))") \
         .Define("u_GEN_y",    "-(pT_muons * TMath::Sin(phi_muons) + GenMET_pt * TMath::Sin(GenMET_phi))") \
         .Define("u_GEN_pt",   "TMath::Sqrt(u_GEN_x * u_GEN_x + u_GEN_y * u_GEN_y)") \
         .Define("u_TK_pt",    "u_tk_eta2p5_pt") \
         .Define("u_TK_x",     "u_tk_eta2p5_pt*TMath::Cos(u_tk_eta2p5_phi)") \
         .Define("u_TK_y",     "u_tk_eta2p5_pt*TMath::Sin(u_tk_eta2p5_phi)") \
         .Define("u_TKPHO_pt", "u_tkpho_eta2p5_eta3p0_pt") \
         .Define("u_TKPHO_x",  "u_tkpho_eta2p5_eta3p0_pt*TMath::Cos(u_tkpho_eta2p5_eta3p0_phi)") \
         .Define("u_TKPHO_y",  "u_tkpho_eta2p5_eta3p0_pt*TMath::Sin(u_tkpho_eta2p5_eta3p0_phi)")

recoils = ["TK", "PUPPI", "PF"]

#for itype in recoils:
#    # prepare the paral, perp, response, diff variables
#    rdf = prepVars(rdf, "u_{RECOIL}".format(RECOIL=itype), "u_GEN")

colors = {
            "TK": 3,
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "TKPHO": 9
         }

labels = {
            "TK": "TK",
            "PF": "PF",
            "PUPPI": "PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
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


DrawHistos(hresponses.values(), [labels[itype] for itype in hresponses.keys()], 0, 60., "q_{T} [GeV]", 0., 1.15, "Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()])

DrawHistos(hresols_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 60., "q_{T} [GeV]", 0, 39.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()])

DrawHistos(hresols_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 60., "q_{T} [GeV]", 0, 32.0, "#sigma u_{#perp} [GeV]", "reco_recoil_resol_perp", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()])

DrawHistos(hresponses_nVtx.values(), [labels[itype] for itype in hresponses_nVtx.keys()], 0, 50., "# Vertices", 0., 1.15, "Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()])

DrawHistos(hresols_paral_diff_VS_nVtx.values(), [labels[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], 0, 50., "# Vertices", 0, 50.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()])

DrawHistos(hresols_perp_VS_nVtx.values(), [labels[itype] for itype in hresols_perp_VS_nVtx.keys()], 0, 50., "# Vertices", 0, 50.0, "#sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()])

##
## Scaled 
##
DrawHistos(hresponsesSc.values(), [labels[itype] for itype in hresponses.keys()], 0, 60., "q_{T} [GeV]", 0., 1.15, "Scaled Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()])

DrawHistos(hresolsSc_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 60., "q_{T} [GeV]", 0, 60.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()])

DrawHistos(hresolsSc_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 60., "q_{T} [GeV]", 0, 50.0, "Scaled #sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()])

DrawHistos(hresponsesSc_nVtx.values(), [labels[itype] for itype in hresponses_nVtx.keys()], 0, 50., "# Vertices", 0., 1.15, "Scaled Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()])

DrawHistos(hresolsSc_paral_diff_VS_nVtx.values(), [labels[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], 0, 50., "# Vertices", 0, 50.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()])

DrawHistos(hresolsSc_perp_VS_nVtx.values(), [labels[itype] for itype in hresols_perp_VS_nVtx.keys()], 0, 50., "# Vertices", 0, 50.0, "Scaled #sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()])


#f1 = ROOT.TFile("root_h/output.root", "RECREATE")
#for h2 in h2ds_perp_VS_qT.values():
#    h2.SetDirectory(f1)
#    h2.Write()
#f1.Close()

raw_input()
