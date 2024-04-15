import ROOT
import sys
sys.path.append("../RecoilResol/CMSPLOTS")
from CMSPLOTS.myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code
from utils.RecoilAnalyzer import RecoilAnalyzer
import argparse

noLumi = False

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("Functions_cc.so")

parser = argparse.ArgumentParser()
parser.add_argument("--era", default="2016", help="Era")
parser.add_argument("--applySc", action="store_true", help="Apply response corrections and include the plots in the output")
parser.add_argument("--no-applySc", action="store_false", dest="applySc", help="Do not apply response corrections and do not include the plots in the output")
parser.set_defaults(applySc=True)

args = parser.parse_args()

doTest = (args.era == "test")
do2016 = (args.era == "2016")

assert doTest + do2016 == 1, "Please specify an era: test, 2016, 2017, 2018"

era = args.era
outdir = f"plots/Data/{era}"

applySc = args.applySc
print("apply Response corrections: ", applySc)

chain = ROOT.TChain("Events")
chain.Add("/afs/cern.ch/work/y/yofeng/public/outputroot/Data.root")
rdf_data = ROOT.ROOT.RDataFrame(chain)

chainMC = ROOT.TChain("Events")
chainMC.Add("/afs/cern.ch/work/y/yofeng/public/outputroot/DY.root")
rdf_MC = ROOT.ROOT.RDataFrame(chainMC)

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
            .Define("u_DeepMETCorr_x", "-(pT_muons*TMath::Cos(phi_muons) + deepmet_pt_corr_central*TMath::Cos(deepmet_phi_corr_central) )") \
            .Define("u_DeepMETCorr_y", "-(pT_muons*TMath::Sin(phi_muons) + deepmet_pt_corr_central*TMath::Sin(deepmet_phi_corr_central) )") \
            .Define("u_DeepMETCorr_pt", "TMath::Sqrt(u_DeepMETCorr_x * u_DeepMETCorr_x + u_DeepMETCorr_y * u_DeepMETCorr_y)")
   return rdf

rdf_data = prepareVars(rdf_data)
rdf_MC = prepareVars(rdf_MC)

recoils = ["PF", "PUPPI", "DeepMET", "DeepMETCorr"]

colors = {
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "DeepMET": 4,
            "DeepMETCorr": 8,
         }

labels = {
            "TK": "TK",
            "PF": "PF",
            "PUPPI": "PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DeepMET": "DeepMET (Raw)",
            "DeepMETCorr": "DeepMET",
         }

linestyles = {}

for itype in recoils:
   colors[itype + "_MC"] = colors[itype]
   
   linestyles[itype] = 1
   linestyles[itype + "_MC"] = 2

xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins() 

#
# data
#
recoilanalyzer = RecoilAnalyzer(rdf_data, recoils, rdfMC = rdf_MC)
recoilanalyzer.prepareVars()
recoilanalyzer.prepareResponses(   'u_GEN_pt', xbins_qT)
recoilanalyzer.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
recoilanalyzer.prepareResponses(   'PV_npvsGood', xbins_nVtx)
recoilanalyzer.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

hresponses = recoilanalyzer.getResponses('u_GEN_pt')
hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions('u_GEN_pt')
hresponses_nVtx = recoilanalyzer.getResponses('PV_npvsGood')
hresols_paral_diff_VS_nVtx, hresols_perp_VS_nVtx = recoilanalyzer.getResolutions('PV_npvsGood')

if applySc:
   ROOT.gInterpreter.Declare(get_response_code)
   # create branch with the scale factors
   for itype in recoils:
       #"dynamic scopes" to create a variable holding histograms
       ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}= {HNAME} ".format(RECOIL=itype, HNAME=hresponses[itype].GetName()))
       rdf_data = rdf_data.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL})'.format(RECOIL=itype)) \
                .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
                .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))
                
       ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}_MC= {HNAME} ".format(RECOIL=itype, HNAME=hresponses[itype+"_MC"].GetName()))
       rdf_MC = rdf_MC.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL}_MC)'.format(RECOIL=itype)) \
                  .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
                  .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))

   #for itype in recoils:
   #    # prepare the paral, perp, response, diff variables
   #    rdf = prepVars(rdf, "u_{RECOIL}Sc".format(RECOIL=itype), "u_GEN")

   recoilsSc = [itype + "Sc" for itype in recoils]
   recoilanalyzerSc = RecoilAnalyzer(rdf_data, recoilsSc, rdfMC = rdf_MC)
   recoilanalyzerSc.prepareVars()
   recoilanalyzerSc.prepareResponses(   'u_GEN_pt', xbins_qT)
   recoilanalyzerSc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
   recoilanalyzerSc.prepareResponses(   'PV_npvsGood', xbins_nVtx)
   recoilanalyzerSc.prepareResolutions( 'PV_npvsGood', xbins_nVtx, 400, -200, 200)

   hresponsesSc = recoilanalyzerSc.getResponses('u_GEN_pt')
   hresolsSc_paral_diff, hresolsSc_perp = recoilanalyzerSc.getResolutions('u_GEN_pt')
   hresponsesSc_nVtx = recoilanalyzerSc.getResponses('PV_npvsGood')
   hresolsSc_paral_diff_VS_nVtx, hresolsSc_perp_VS_nVtx = recoilanalyzerSc.getResolutions('PV_npvsGood')

qtmax = 150.0

def GetLegends(hdict):
    return [labels[itype] for itype in hdict.keys() if "_MC" not in itype]
 
def GetLineStyles(hdict):
    return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]

DrawHistos(hresponses.values(), GetLegends(hresponses), 0, qtmax, "q_{T} [GeV]", 0., 1.15, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles = GetLineStyles(hresponses), noLumi=noLumi, outdir=outdir)

DrawHistos(hresols_paral_diff.values(), GetLegends(hresols_paral_diff), 0, qtmax, "q_{T} [GeV]", 0, 39.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff))

DrawHistos(hresols_perp.values(), GetLegends(hresols_perp), 0, qtmax, "q_{T} [GeV]", 0, 32.0, "#sigma u_{#perp } [GeV]", "reco_recoil_resol_perp", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp))

DrawHistos(hresponses_nVtx.values(), GetLegends(hresponses_nVtx), 0, 50., "# Vertices", 0., 1.15, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses_nVtx))

DrawHistos(hresols_paral_diff_VS_nVtx.values(), GetLegends(hresols_paral_diff_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "#sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff_VS_nVtx))

DrawHistos(hresols_perp_VS_nVtx.values(), GetLegends(hresols_perp_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "#sigma u_{#perp } [GeV]", "reco_recoil_resol_perp_VS_nVtx", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp_VS_nVtx))

if applySc:
   #
   # Scaled 
   #
   DrawHistos(hresponsesSc.values(), GetLegends(hresponses), 0, qtmax, "q_{T} [GeV]", 0., 1.15, "Scaled Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses))

   DrawHistos(hresolsSc_paral_diff.values(), GetLegends(hresols_paral_diff), 0, qtmax, "q_{T} [GeV]", 0, 60.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff))

   DrawHistos(hresolsSc_perp.values(), GetLegends(hresols_perp), 0, qtmax, "q_{T} [GeV]", 0, 50.0, "Scaled #sigma u_{#perp } [GeV]", "reco_recoil_resol_perp_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp))

   DrawHistos(hresponsesSc_nVtx.values(), GetLegends(hresponses_nVtx), 0, 50., "# Vertices", 0., 1.15, "Scaled Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses_nVtx))

   DrawHistos(hresolsSc_paral_diff_VS_nVtx.values(), GetLegends(hresols_paral_diff_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "Scaled #sigma u_{#parallel} [GeV]", "reco_recoil_resol_paral_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.38, 0.92], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff_VS_nVtx))

   DrawHistos(hresolsSc_perp_VS_nVtx.values(), GetLegends(hresols_perp_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "Scaled #sigma u_{#perp } [GeV]", "reco_recoil_resol_perp_VS_nVtx_Scaled", drawashist=True, dology=False, legendPos=[0.20, 0.73, 0.40, 0.92], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp_VS_nVtx))


#f1 = ROOT.TFile("root_h/output.root", "RECREATE")
#for h2 in h2ds_perp_VS_qT.values():
#    h2.SetDirectory(f1)
#    h2.Write()
#f1.Close()
