'''
compare the performances of PF, Puppi, Deep (w/o Puppi weights in the training).
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
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_1_NoPUPPI.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_2_NoPUPPI.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_3_NoPUPPI.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_4_NoPUPPI.root")
itreef.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_9_NoPUPPI.root")

itreef.AddFriend(itree)

rdf_WoPUPPI = ROOT.ROOT.RDataFrame( itreef )

rdf_WoPUPPI = rdf_WoPUPPI \
         .Define("u_GEN_x",    "-(pT_muons * TMath::Cos(phi_muons) + GenMET_pt * TMath::Cos(GenMET_phi))") \
         .Define("u_GEN_y",    "-(pT_muons * TMath::Sin(phi_muons) + GenMET_pt * TMath::Sin(GenMET_phi))") \
         .Define("u_GEN_pt",   "TMath::Sqrt(u_GEN_x * u_GEN_x + u_GEN_y * u_GEN_y)") \
         .Define("u_WoPUPPI_x",  "u_pred_x") \
         .Define("u_WoPUPPI_y",  "u_pred_y") \
         .Define("u_WoPUPPI_pt",  "TMath::Sqrt(u_pred_x * u_pred_x + u_pred_y * u_pred_y)")

# with PUPPI
itreef2 = ROOT.TChain("tree")
itreef2.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_1.root")
itreef2.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_2.root")
itreef2.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_3.root")
itreef2.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_4.root")
itreef2.Add("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/DeepPFMET/root/result_DYJets_2016_9.root")

itreef2.AddFriend(itree)

rdf_WPUPPI = ROOT.ROOT.RDataFrame( itreef2 )

rdf_WPUPPI = rdf_WPUPPI.Define("u_PUPPI_pt", "calUT_PtPhi(pT_muons,   phi_muons, PuppiMET_pt, PuppiMET_phi)") \
         .Define("u_PUPPI_x",  "-(pT_muons * TMath::Cos(phi_muons) + PuppiMET_pt*TMath::Cos(PuppiMET_phi))") \
         .Define("u_PUPPI_y",  "-(pT_muons * TMath::Sin(phi_muons) + PuppiMET_pt*TMath::Sin(PuppiMET_phi))") \
         .Define("u_PF_pt",    "calUT_PtPhi(pT_muons,   phi_muons, MET_pt,      MET_phi)") \
         .Define("u_PF_x",     "-(pT_muons * TMath::Cos(phi_muons) + MET_pt*TMath::Cos(MET_phi))") \
         .Define("u_PF_y",     "-(pT_muons * TMath::Sin(phi_muons) + MET_pt*TMath::Sin(MET_phi))") \
         .Define("u_GEN_x",    "-(pT_muons * TMath::Cos(phi_muons) + GenMET_pt * TMath::Cos(GenMET_phi))") \
         .Define("u_GEN_y",    "-(pT_muons * TMath::Sin(phi_muons) + GenMET_pt * TMath::Sin(GenMET_phi))") \
         .Define("u_GEN_pt",   "TMath::Sqrt(u_GEN_x * u_GEN_x + u_GEN_y * u_GEN_y)") \
         .Define("u_TK_pt",    "u_tk_eta2p5_pt") \
         .Define("u_TK_x",     "u_tk_eta2p5_pt*TMath::Cos(u_tk_eta2p5_phi)") \
         .Define("u_TK_y",     "u_tk_eta2p5_pt*TMath::Sin(u_tk_eta2p5_phi)") \
         .Define("u_DEEP_x",  "u_pred_x") \
         .Define("u_DEEP_y",  "u_pred_y") \
         .Define("u_DEEP_pt",  "TMath::Sqrt(u_pred_x * u_pred_x + u_pred_y * u_pred_y)")

recoils = ["PF", "PUPPI", "TK", "DEEP"]

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
            "WoPUPPI": 4,
         }

labels = {
            "TK": "TK",
            "PF": "Type-I PF",
            "PUPPI": "Type-I PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DEEP": "DeepMET",
            "WoPUPPI": "DeepMET w/o PUPPI"
         }

xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins() 

recoilanalyzer = RecoilAnalyzer(rdf_WPUPPI, recoils)
recoilanalyzer.prepareVars()
recoilanalyzer.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)

hresponses = recoilanalyzer.getResponses('u_GEN_pt')
hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions('u_GEN_pt')

ROOT.gInterpreter.Declare(get_response_code)
# create branch with the scale factors
for itype in recoils:
    #"dynamic scopes" to create a variable holding histograms
    ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}= {HNAME} ".format(RECOIL=itype, HNAME=hresponses[itype].GetName()))
    rdf_WPUPPI = rdf_WPUPPI.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL})'.format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

recoilsSc = [itype + "Sc" for itype in recoils]
recoilanalyzerSc = RecoilAnalyzer(rdf_WPUPPI, recoilsSc)
recoilanalyzerSc.prepareVars()
recoilanalyzerSc.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzerSc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)

hresponsesSc = recoilanalyzerSc.getResponses('u_GEN_pt')
hresolsSc_paral_diff, hresolsSc_perp = recoilanalyzerSc.getResolutions('u_GEN_pt')

# this is really stupid
# I just need the NoPUPPI columns but have to do everything once again...
recoils = ['WoPUPPI']
recoilanalyzer = RecoilAnalyzer(rdf_WoPUPPI, recoils)
recoilanalyzer.prepareVars()
recoilanalyzer.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)

hresponses_WoPUPPI = recoilanalyzer.getResponses('u_GEN_pt')
hresols_paral_diff_WoPUPPI, hresols_perp_WoPUPPI = recoilanalyzer.getResolutions('u_GEN_pt')

# create branch with the scale factors
for itype in recoils:
    #"dynamic scopes" to create a variable holding histograms
    ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}= {HNAME} ".format(RECOIL=itype, HNAME=hresponses_WoPUPPI[itype].GetName()))
    rdf_WoPUPPI = rdf_WoPUPPI.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL})'.format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

recoilsSc = [itype + "Sc" for itype in recoils]
recoilanalyzerSc = RecoilAnalyzer(rdf_WoPUPPI, recoilsSc)
recoilanalyzerSc.prepareVars()
recoilanalyzerSc.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzerSc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)

hresponsesSc_WoPUPPI = recoilanalyzerSc.getResponses('u_GEN_pt')
hresolsSc_paral_diff_WoPUPPI, hresolsSc_perp_WoPUPPI = recoilanalyzerSc.getResolutions('u_GEN_pt')

hresponses.update(hresponses_WoPUPPI)
hresols_paral_diff.update(hresols_paral_diff_WoPUPPI)
hresols_perp.update(hresols_perp_WoPUPPI)
hresponsesSc.update(hresponsesSc_WoPUPPI)
hresolsSc_paral_diff.update(hresolsSc_paral_diff_WoPUPPI)
hresolsSc_perp.update(hresolsSc_perp_WoPUPPI)

DrawHistos(hresponses.values(), [labels[itype] for itype in hresponses.keys()], 0, 150., "q_{T} [GeV]", 0., 1.15, "Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response", drawashist=True, dology=False, legendPos=[0.55, 0.17, 0.73, 0.44], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[1,1,1,1, 2])

DrawHistos(hresols_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 150., "q_{T} [GeV]", 0, 39.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], linestyles=[1,1,1,1, 2])

DrawHistos(hresols_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 150., "q_{T} [GeV]", 0, 32.0, "#sigma u_{#perp} [GeV]", "reco_recoil_resol_perp", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()], linestyles=[1,1,1,1, 2])

##
## Scaled 
##
DrawHistos(hresponsesSc.values(), [labels[itype] for itype in hresponses.keys()], 0, 150., "q_{T} [GeV]", 0., 1.15, "Scaled Reponse -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.44], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[1,1,1,1, 2])

DrawHistos(hresolsSc_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 150., "q_{T} [GeV]", 0, 50.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_Scaled", drawashist=True, dology=False, legendPos=[0.55, 0.65, 0.73, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], linestyles=[1,1,1,1,2])

DrawHistos(hresolsSc_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 150., "q_{T} [GeV]", 0, 45.0, "Scaled #sigma u_{#perp} [GeV]", "reco_recoil_resol_perp_Scaled", drawashist=True, dology=False, legendPos=[0.55, 0.65, 0.75, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()], linestyles=[1,1,1,1,2])

#f1 = ROOT.TFile("root_h/output.root", "RECREATE")
#for h2 in h2ds_perp_VS_qT.values():
#    h2.SetDirectory(f1)
#    h2.Write()
#f1.Close()

raw_input()
