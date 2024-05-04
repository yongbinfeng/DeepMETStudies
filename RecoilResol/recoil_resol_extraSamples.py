import ROOT
import sys
sys.path.append("../RecoilResol/CMSPLOTS")
from CMSPLOTS.myFunction import DrawHistos, getMedian, getResolution
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code
import argparse

noLumi = False
MCOnly = True

suffix = "_ttHbb"

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

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
outdir = f"plots/MC_CompPUPPI/{era}"

applySc = args.applySc
print("apply Response corrections: ", applySc)

chainMC = ROOT.TChain("Events")
chainMC.Add("/eos/cms/store/user/yofeng/DeepMETNanoAOD/ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8.root")
rdf_MC = ROOT.ROOT.RDataFrame(chainMC)

def prepareVars(rdf):
   rdf = rdf.Define("GenHT", "Sum(GenJet_pt)")  \
            .Define("PFMET_pt", "MET_pt") \
            .Define("PFMET_phi", "MET_phi") \
            .Define("PUPPIMET_pt", "PuppiMET_pt") \
            .Define("PUPPIMET_phi", "PuppiMET_phi") \
            .Define("DeepMETMET_pt", "DeepMETResolutionTune_pt") \
            .Define("DeepMETMET_phi", "DeepMETResolutionTune_phi")
   return rdf


rdf_MC = prepareVars(rdf_MC)

recoils = ["PF", "PUPPI", "DeepMET"]

for recoil in recoils:
   rdf_MC = rdf_MC.Define(f"{recoil}MET_diff",  "{RECOIL}MET_pt - GenMET_pt".format(RECOIL=recoil))
   rdf_MC = rdf_MC.Define(f"{recoil}MET_ratio", "{RECOIL}MET_pt / GenMET_pt".format(RECOIL=recoil))
   rdf_MC = rdf_MC.Define(f"{recoil}MET_ratio_diff", "{RECOIL}MET_ratio - 1".format(RECOIL=recoil))

histos_2d_diffs = OrderedDict()
histos_ratios = OrderedDict()
histos_ratios_diff = OrderedDict()
for recoil in recoils:
   hname = "h2d_diff_{RECOIL}_vs_HT".format(RECOIL=recoil)
   histos_2d_diffs[recoil] = rdf_MC.Histo2D((hname, hname, 50, 0, 1000, 400, -200, 200), "GenHT", "{RECOIL}MET_diff".format(RECOIL=recoil))
   hname = "h2d_ratio_{RECOIL}_vs_HT".format(RECOIL=recoil)
   histos_ratios[recoil] = rdf_MC.Histo2D((hname, hname, 50, 0, 1000, 400, 0.5, 1.5), "GenHT", "{RECOIL}MET_ratio".format(RECOIL=recoil))
   hname = "h2d_ratio_diff_{RECOIL}_vs_HT".format(RECOIL=recoil)
   histos_ratios_diff[recoil] = rdf_MC.Histo2D((hname, hname, 50, 0, 1000, 400, -0.5, 0.5), "GenHT", "{RECOIL}MET_ratio_diff".format(RECOIL=recoil))
   
hresols = OrderedDict()
hresols_ratio = OrderedDict()
hresponses = OrderedDict()
for recoil in recoils:
   hresponses[recoil] = getMedian(histos_ratios[recoil]) 
   hresols[recoil] = getResolution(histos_2d_diffs[recoil])
   hresols_ratio[recoil] = getResolution(histos_ratios_diff[recoil])

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

linestyles = {}

for itype in recoils:
   colors[itype + "_MC"] = colors[itype]
   
   linestyles[itype] = 1
   linestyles[itype + "_MC"] = 2

xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins() 

# loop over the different qT bins
if applySc:
   ROOT.gInterpreter.Declare(get_response_code)
   # create branch with the scale factors
   for itype in recoils:
       #"dynamic scopes" to create a variable holding histograms
       ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}{suffix}= {HNAME} ".format(RECOIL=itype, suffix = suffix, HNAME=hresponses[itype].GetName()))

       rdf_MC = rdf_MC.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(GenHT, hprof_{RECOIL}{suffix})'.format(RECOIL=itype, suffix=suffix)) \
                  .Define("{RECOIL}MET_Sc_pt".format(RECOIL=itype), "{RECOIL}_scale * {RECOIL}MET_pt".format(RECOIL=itype)) \
                  .Define("{RECOIL}MET_Sc_diff".format(RECOIL=itype), "{RECOIL}MET_Sc_pt - GenMET_pt".format(RECOIL=itype)) \
                  .Define("{RECOIL}MET_Sc_ratio".format(RECOIL=itype), "{RECOIL}MET_Sc_pt / GenMET_pt".format(RECOIL=itype)) \
                  .Define("{RECOIL}MET_Sc_ratio_diff".format(RECOIL=itype), "{RECOIL}MET_Sc_ratio - 1".format(RECOIL=itype))
                  
   histos_2d_diffs_sc = OrderedDict()
   histos_ratios_sc = OrderedDict()
   histos_ratios_diff_sc = OrderedDict()
   
   for recoil in recoils:
      hname = "h2d_diff_{RECOIL}Sc_vs_HT".format(RECOIL=recoil)
      histos_2d_diffs_sc[recoil] = rdf_MC.Histo2D((hname, hname, 50, 0, 1000, 400, -200, 200), "GenHT", "{RECOIL}MET_Sc_diff".format(RECOIL=recoil))
      hname = "h2d_ratio_{RECOIL}Sc_vs_HT".format(RECOIL=recoil)
      histos_ratios_sc[recoil] = rdf_MC.Histo2D((hname, hname, 50, 0, 1000, 400, 0.5, 1.5), "GenHT", "{RECOIL}MET_Sc_ratio".format(RECOIL=recoil))
      hname = "h2d_ratio_diff_{RECOIL}Sc_vs_HT".format(RECOIL=recoil)
      histos_ratios_diff_sc[recoil] = rdf_MC.Histo2D((hname, hname, 50, 0, 1000, 400, -0.5, 0.5), "GenHT", "{RECOIL}MET_Sc_ratio_diff".format(RECOIL=recoil))

   hresols_sc = OrderedDict()
   hresponses_sc = OrderedDict()
   hresols_ratio_sc = OrderedDict()
   for recoil in recoils:
      hresponses_sc[recoil] = getMedian(histos_ratios_sc[recoil]) 
      hresols_sc[recoil] = getResolution(histos_2d_diffs_sc[recoil])
      hresols_ratio_sc[recoil] = getResolution(histos_ratios_diff_sc[recoil])


qtmax = 1000.0

def GetLegends(hdict):
    return [labels[itype] for itype in hdict.keys() if "_MC" not in itype]

def GetLineStyles(hdict):
    return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]

DrawHistos(hresponses.values(), GetLegends(hresponses), 0, qtmax, "q_{T} [GeV]", 0., 1.20, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response" + suffix, drawashist=True, dology=False, legendPos=[0.50, 0.20, 0.80, 0.40], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles = GetLineStyles(hresponses), noLumi=noLumi, outdir=outdir, MCOnly=MCOnly)

DrawHistos(hresols.values(), GetLegends(hresols), 0, qtmax, "q_{T} [GeV]", 0, 50.0, "#sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols), MCOnly=MCOnly)

DrawHistos(hresols_ratio.values(), GetLegends(hresols_ratio), 0, qtmax, "q_{T} [GeV]", 0, 0.5, "#sigma (u_{#parallel}) / <u_{#parallel}>", "reco_recoil_resol_paral_ratio" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_ratio.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_ratio), MCOnly=MCOnly)


if applySc:
   #
   # Scaled 
   #
   DrawHistos(hresponses_sc.values(), GetLegends(hresponses), 0, qtmax, "q_{T} [GeV]", 0., 1.20, "Scaled Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses), MCOnly=MCOnly)

   DrawHistos(hresols_sc.values(), GetLegends(hresols), 0, qtmax, "q_{T} [GeV]", 0, 50.0, "Response-corrected #sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols), MCOnly=MCOnly)
   
   DrawHistos(hresols_ratio_sc.values(), GetLegends(hresols_ratio), 0, qtmax, "q_{T} [GeV]", 0, 0.5, "Response-corrected #sigma (u_{#parallel}) / <u_{#parallel}>", "reco_recoil_resol_paral_ratio_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_ratio.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_ratio), MCOnly=MCOnly)
   
   DrawHistos([hresponses["DeepMET"], hresponses_sc["DeepMET"]], ["DeepMET", "DeepMET Scaled"], 0, qtmax, "q_{T} [GeV]", 0., 1.20, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.50, 0.20, 0.80, 0.40], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, showratio=True, yrlabel="Scaled/Unscaled", yrmin=0.91, yrmax=1.09)
   
   DrawHistos([hresols["DeepMET"], hresols_sc["DeepMET"]], ["DeepMET", "DeepMET Scaled"], 0, qtmax, "q_{T} [GeV]", 0, 50.0, "#sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, showratio=True, yrlabel="Scaled/Unscaled", yrmin=0.91, yrmax=1.09)
   
   DrawHistos([hresols_ratio["DeepMET"], hresols_ratio_sc["DeepMET"]], ["DeepMET", "DeepMET Scaled"], 0, qtmax, "q_{T} [GeV]", 0, 0.5, "#sigma (u_{#parallel}) / <u_{#parallel}>", "reco_recoil_resol_paral_ratio_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, showratio=True, yrlabel="Scaled/Unscaled", yrmin=0.91, yrmax=1.09)
   
   