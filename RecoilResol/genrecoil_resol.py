import ROOT
import sys
sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
from myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, get_response_code
from utils.RecoilAnalyzer import RecoilAnalyzer

ROOT.gROOT.SetBatch(1)

recoils = [
            "Gen_eta5p0", "Gen_pt500_eta5p0",
            "Gen_eta3p0", "Gen_pt500_eta3p0",
            "tkphoGen_eta2p5_eta3p0", "tkphoGen_pt500_eta2p5_eta3p0",
            "tkGen_eta2p5", "tkGen_pt500_eta2p5",
          ]

colors = {
            "Gen_eta5p0": 1, "Gen_pt500_eta5p0": 1,
            "Gen_eta3p0": 4, "Gen_pt500_eta3p0": 4,
            "tkphoGen_eta2p5_eta3p0": 2, "tkphoGen_pt500_eta2p5_eta3p0": 2,
            "tkGen_eta2p5": 3, "tkGen_pt500_eta2p5": 3,
         }

linestyles = {
                "Gen_eta5p0": 1,  "Gen_pt500_eta5p0": 2,
                "Gen_eta3p0": 1,  "Gen_pt500_eta3p0": 2,
                "tkphoGen_eta2p5_eta3p0": 1, "tkphoGen_pt500_eta2p5_eta3p0": 2,
                "tkGen_eta2p5": 1, "tkGen_pt500_eta2p5": 2,
             }

labels = {
            "Gen_eta5p0": "|#eta|<5.0",
            "Gen_eta3p0": "|#eta|<3.0",
            "tkphoGen_eta2p5_eta3p0": "Chg+Photons",
            "tkGen_eta2p5": "Chg",
            "Gen_pt500_eta5p0": "|#eta|<5.0, p_{T}>0.5",
            "Gen_pt500_eta3p0": "|#eta|<3.0, p_{T}>0.5",
            "tkphoGen_pt500_eta2p5_eta3p0": "Chg+Photons, p_{T}>0.5",
            "tkGen_pt500_eta2p5": "Chg, p_{T}>0.5",
         }

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("/afs/cern.ch/work/y/yofeng/public/CMSSW_10_6_0/src/ANPlots/Functions_cc.so")

rdf = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/user/yofeng/WRecoilNanoAOD_Skimmed_v10_tempMET/myNanoProdMc2016_NANO_[1-9]_Skim.root")
#rdf = rdf_org.Filter("pT_W<30.0")

rdf = rdf.Define("u_GEN_x",    "-(pT_muons * TMath::Cos(phi_muons) + GenMET_pt * TMath::Cos(GenMET_phi))") \
         .Define("u_GEN_y",    "-(pT_muons * TMath::Sin(phi_muons) + GenMET_pt * TMath::Sin(GenMET_phi))") \
         .Define("u_GEN_pt",   "TMath::Sqrt(u_GEN_x * u_GEN_x + u_GEN_y * u_GEN_y)") \

for itype in recoils:
    rdf = rdf.Define( "u_{RECOIL}_x".format(RECOIL=itype), "u_{RECOIL}_pt * TMath::Cos(u_{RECOIL}_phi)".format(RECOIL=itype) ) \
             .Define( "u_{RECOIL}_y".format(RECOIL=itype), "u_{RECOIL}_pt * TMath::Sin(u_{RECOIL}_phi)".format(RECOIL=itype) )

xbins_qT = getpTBins()

recoilanalyzer = RecoilAnalyzer(rdf, recoils)
recoilanalyzer.prepareVars()
recoilanalyzer.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzer.prepareResolutions1D(400, -200, 200)

hresponses = recoilanalyzer.getResponses('u_GEN_pt')
hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions('u_GEN_pt')
vresols_paral_diff, vresols_perp = recoilanalyzer.getResolutions1D()

ROOT.gInterpreter.Declare(get_response_code)
# create branch with the scale factors
for itype in recoils:
    #"dynamic scopes" to create a variable holding histograms
    ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}= {HNAME} ".format(RECOIL=itype, HNAME=hresponses[itype].GetName()))
    rdf = rdf.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL})'.format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
             .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

recoilsSc = [itype + "Sc" for itype in recoils]
recoilanalyzerSc = RecoilAnalyzer(rdf, recoilsSc)
recoilanalyzerSc.prepareVars()
recoilanalyzerSc.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzerSc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzerSc.prepareResolutions1D(400, -200, 200)

hresponsesSc = recoilanalyzerSc.getResponses('u_GEN_pt')
hresolsSc_paral_diff, hresolsSc_perp = recoilanalyzerSc.getResolutions('u_GEN_pt')
vresolsSc_paral_diff, vresolsSc_perp = recoilanalyzerSc.getResolutions1D()

DrawHistos(hresponses.values(), [labels[itype] for itype in hresponses.keys()], 0, 60., "q_{T} [GeV]", 0., 1.15, "Reponse -<u_{#parallel}>/<q_{T}>", "recoil_response", drawashist=True, dology=False, legendPos=[0.22, 0.17, 0.90, 0.36], legendNCols=2, mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[linestyles[itype] for itype in hresponses.keys()])

DrawHistos(hresols_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 60., "q_{T} [GeV]", 0, 30.0, "#sigma u_{#parallel} [GeV]", "recoil_resol_paral", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.88, 0.92], legendNCols=2, mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], linestyles=[linestyles[itype] for itype in hresols_paral_diff.keys()])

DrawHistos(hresols_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 60., "q_{T} [GeV]", 0, 20.0, "#sigma u_{#perp} [GeV]", "recoil_resol_perp", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.88, 0.92], legendNCols=2, mycolors=[colors[itype] for itype in hresols_perp.keys()], linestyles=[linestyles[itype] for itype in hresols_perp.keys()])

##
## Scaled 
##
DrawHistos(hresponsesSc.values(), [labels[itype] for itype in hresponses.keys()], 0, 60., "q_{T} [GeV]", 0., 1.15, "Scaled Reponse -<u_{#parallel}>/<q_{T}>", "recoil_response_scaled", drawashist=True, dology=False, legendPos=[0.22, 0.17, 0.90, 0.36], legendNCols=2, mycolors=[colors[itype] for itype in hresponses.keys()], linestyles=[linestyles[itype] for itype in hresponses.keys()])

DrawHistos(hresolsSc_paral_diff.values(), [labels[itype] for itype in hresols_paral_diff.keys()], 0, 60., "q_{T} [GeV]", 0, 50.0, "Scaled #sigma u_{#parallel} [GeV]", "recoil_resol_paral_scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.88, 0.92], legendNCols=2, mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], linestyles=[linestyles[itype] for itype in hresols_paral_diff.keys()])

DrawHistos(hresolsSc_perp.values(), [labels[itype] for itype in hresols_perp.keys()], 0, 60., "q_{T} [GeV]", 0, 30.0, "Scaled #sigma u_{#perp} [GeV]", "recoil_resol_perp_scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.88, 0.92], legendNCols=2, mycolors=[colors[itype] for itype in hresols_perp.keys()], linestyles=[linestyles[itype] for itype in hresols_perp.keys()])

for itype in vresols_paral_diff:
    print "itype \t %25s \t %f \t %f"%(itype, vresols_paral_diff[itype][1], vresols_perp[itype][1])

print "*"*50

for itype in vresolsSc_paral_diff:
    print "itype \t %25s \t %f \t %f"%(itype, vresolsSc_paral_diff[itype][1], vresolsSc_perp[itype][1])

#f1 = ROOT.TFile("root_h/output.root", "RECREATE")
#for h2 in h2ds_perp_VS_qT.values():
#    h2.SetDirectory(f1)
#    h2.Write()
#f1.Close()

raw_input()
