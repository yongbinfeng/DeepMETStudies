import ROOT
import sys
sys.path.append("../RecoilResol/CMSPLOTS")
from CMSPLOTS.myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code, prepRecoilVars, nPVString
from utils.RecoilAnalyzer import RecoilAnalyzer
import argparse

ROOT.gROOT.SetBatch(True)

ROOT.ROOT.EnableImplicitMT(10)

ROOT.gSystem.Load("Functions_cc.so")

parser = argparse.ArgumentParser()
parser.add_argument("--era", default="2016", help="Era")
args = parser.parse_args()

era = args.era
outdir = f"plots/MC_OldSample/{era}"

rdf_org = ROOT.ROOT.RDataFrame("Events", "/eos/cms/store/user/yofeng/DeepMETNanoAOD/2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_PUMoriond17_Nano02Apr2020_102X_mcRun2_asymptotic_v8_ext2-v1_RunIISummer16NanoAODv7.root")
rdf_org1 = rdf_org.Filter("nMuon >= 2").Filter("HLT_IsoMu24 || HLT_IsoTkMu24")
rdf_org1 = rdf_org1.Define("Muon_pass0", "Muon_pt[0] > 25.0 && abs(Muon_eta[0]) < 2.4 && Muon_pfRelIso04_all[0] < 0.15 && Muon_looseId[0]").Define("Muon_pass1", "Muon_pt[1] > 25.0 && abs(Muon_eta[1]) < 2.4 && Muon_pfRelIso04_all[1] < 0.15 && Muon_looseId[1]")
rdf_org2 = rdf_org1.Filter("Muon_pass0").Filter("Muon_pass1")
rdf = rdf_org2

rdf = rdf.Define("V_pt", "TMath::Sqrt(Muon_pt[0] * Muon_pt[0] + Muon_pt[1] * Muon_pt[1] + 2 * Muon_pt[0] * Muon_pt[1] * TMath::Cos(Muon_phi[0] - Muon_phi[1]))")
rdf = rdf.Define("V_phi", "TMath::ATan2(Muon_pt[0] * TMath::Sin(Muon_phi[0]) + Muon_pt[1] * TMath::Sin(Muon_phi[1]), Muon_pt[0] * TMath::Cos(Muon_phi[0]) + Muon_pt[1] * TMath::Cos(Muon_phi[1]))")

rdf = rdf.Define("PF_pt", "MET_pt") \
         .Define("PF_phi", "MET_phi") \
         .Define("PUPPI_pt", "PuppiMET_pt") \
         .Define("PUPPI_phi", "PuppiMET_phi")

recoils = ["PF", "PUPPI"]
rdf = prepRecoilVars(rdf, "V", recoils)
rdf = rdf.Define("u_GEN_pt", "V_pt").Define("u_GEN_x", "- V_pt * TMath::Cos(V_phi)").Define("u_GEN_y", "- V_pt * TMath::Sin(V_phi)")

rdf_qTLow = rdf.Filter("V_pt < 50.0")
rdf_qTHigh = rdf.Filter("V_pt >= 50.0")

suffixes = ['_OldSample', '_OldSample_qTLow', '_OldSample_qTHigh']

xbins_qT = getpTBins()
xbins_nVtx = getnVtxBins()
nPV = nPVString()

colors = {
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "DeepMET": 4,
         }

labels = {
            "TK": "TK",
            "PF": "PF",
            "PUPPI": "PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DeepMET": "DeepMET",
         }

applySc = True 
noLumi = False
extraText = "Simulation (2016 Pre-UL)"
MCOnly = True

for idx, rdf in enumerate([rdf, rdf_qTLow, rdf_qTHigh]):
   suffix = suffixes[idx]
   recoilanalyzer = RecoilAnalyzer(rdf, recoils, name = "recoilanalyzer" + suffix)
   recoilanalyzer.prepareVars()
   recoilanalyzer.prepareResponses(   'u_GEN_pt', xbins_qT)
   recoilanalyzer.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
   recoilanalyzer.prepareResponses(   nPV, xbins_nVtx)
   recoilanalyzer.prepareResolutions( nPV, xbins_nVtx, 400, -200, 200)

   hresponses = recoilanalyzer.getResponses('u_GEN_pt')
   hresols_paral_diff, hresols_perp = recoilanalyzer.getResolutions('u_GEN_pt')
   hresponses_nVtx = recoilanalyzer.getResponses(nPV)
   hresols_paral_diff_VS_nVtx, hresols_perp_VS_nVtx = recoilanalyzer.getResolutions(nPV)

   if applySc:
      if idx == 0:
         ROOT.gInterpreter.Declare(get_response_code)
      # create branch with the scale factors
      for itype in recoils:
          #"dynamic scopes" to create a variable holding histograms
          ROOT.gInterpreter.ProcessLine("auto hprof_{RECOIL}{suffix}= {HNAME} ".format(RECOIL=itype, suffix = suffix, HNAME=hresponses[itype].GetName()))
          rdf = rdf.Define("{RECOIL}_scale".format(RECOIL=itype), '1.0/get_response(u_GEN_pt, hprof_{RECOIL}{suffix})'.format(RECOIL=itype, suffix=suffix)) \
                   .Define("u_{RECOIL}Sc_x".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_x".format(RECOIL=itype)) \
                   .Define("u_{RECOIL}Sc_y".format(RECOIL=itype), "{RECOIL}_scale * u_{RECOIL}_y".format(RECOIL=itype))


      recoilsSc = [itype + "Sc" for itype in recoils]
      recoilanalyzerSc = RecoilAnalyzer(rdf, recoilsSc, name = "recoilanalyzer_Scaled" + suffix)
      recoilanalyzerSc.prepareVars()
      recoilanalyzerSc.prepareResponses(   'u_GEN_pt', xbins_qT)
      recoilanalyzerSc.prepareResolutions( 'u_GEN_pt', xbins_qT, 400, -200, 200)
      recoilanalyzerSc.prepareResponses(   nPV, xbins_nVtx)
      recoilanalyzerSc.prepareResolutions( nPV, xbins_nVtx, 400, -200, 200)

      hresponsesSc = recoilanalyzerSc.getResponses('u_GEN_pt')
      hresolsSc_paral_diff, hresolsSc_perp = recoilanalyzerSc.getResolutions('u_GEN_pt')
      hresponsesSc_nVtx = recoilanalyzerSc.getResponses(nPV)
      hresolsSc_paral_diff_VS_nVtx, hresolsSc_perp_VS_nVtx = recoilanalyzerSc.getResolutions(nPV)


   qtmax = 150.0

   def GetLegends(hdict):
       return [labels[itype] for itype in hdict.keys() if "_MC" not in itype]

   def GetLineStyles(hdict):
       return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]

   DrawHistos(hresponses.values(), GetLegends(hresponses), 0, qtmax, "q_{T} [GeV]", 0., 1.15, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response" + suffix, drawashist=True, dology=False, legendPos=[0.50, 0.20, 0.80, 0.40], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles = GetLineStyles(hresponses), noLumi=noLumi, outdir=outdir, extraText=extraText, MCOnly=MCOnly)

   DrawHistos(hresols_paral_diff.values(), GetLegends(hresols_paral_diff), 0, qtmax, "q_{T} [GeV]", 0, 39.0, "#sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff), extraText=extraText, MCOnly=MCOnly)

   DrawHistos(hresols_perp.values(), GetLegends(hresols_perp), 0, qtmax, "q_{T} [GeV]", 0, 32.0, "#sigma (u_{#perp } ) [GeV]", "reco_recoil_resol_perp" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp), extraText=extraText, MCOnly=MCOnly)

   DrawHistos(hresponses_nVtx.values(), GetLegends(hresponses_nVtx), 0, 50., "# Vertices", 0., 1.15, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx" + suffix, drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses_nVtx), extraText=extraText, MCOnly=MCOnly)

   DrawHistos(hresols_paral_diff_VS_nVtx.values(), GetLegends(hresols_paral_diff_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "#sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral_VS_nVtx" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff_VS_nVtx), extraText=extraText, MCOnly=MCOnly)

   DrawHistos(hresols_perp_VS_nVtx.values(), GetLegends(hresols_perp_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "#sigma (u_{#perp } ) [GeV]", "reco_recoil_resol_perp_VS_nVtx" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp_VS_nVtx), extraText=extraText, MCOnly=MCOnly)

   if applySc:
      #
      # Scaled 
      #
      DrawHistos(hresponsesSc.values(), GetLegends(hresponses), 0, qtmax, "q_{T} [GeV]", 0., 1.15, "Scaled Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses), extraText=extraText, MCOnly=MCOnly)

      DrawHistos(hresolsSc_paral_diff.values(), GetLegends(hresols_paral_diff), 0, qtmax, "q_{T} [GeV]", 0, 60.0, "Response-corrected #sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff), extraText=extraText, MCOnly=MCOnly)

      DrawHistos(hresolsSc_perp.values(), GetLegends(hresols_perp), 0, qtmax, "q_{T} [GeV]", 0, 50.0, "Response-corrected #sigma (u_{#perp } ) [GeV]", "reco_recoil_resol_perp_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp), extraText=extraText, MCOnly=MCOnly)

      DrawHistos(hresponsesSc_nVtx.values(), GetLegends(hresponses_nVtx), 0, 50., "# Vertices", 0., 1.15, "Response -<u_{#parallel}>/<q_{T}>", "reco_recoil_response_VS_nVtx_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses_nVtx), extraText=extraText, MCOnly=MCOnly)

      DrawHistos(hresolsSc_paral_diff_VS_nVtx.values(), GetLegends(hresols_paral_diff_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "Response-corrected #sigma (u_{#parallel}) [GeV]", "reco_recoil_resol_paral_VS_nVtx_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_paral_diff_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_paral_diff_VS_nVtx), extraText=extraText, MCOnly=MCOnly)

      DrawHistos(hresolsSc_perp_VS_nVtx.values(), GetLegends(hresols_perp_VS_nVtx), 0, 50., "# Vertices", 0, 50.0, "Response-corrected #sigma (u_{#perp } ) [GeV]", "reco_recoil_resol_perp_VS_nVtx_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.40, 0.88], mycolors=[colors[itype] for itype in hresols_perp_VS_nVtx.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_perp_VS_nVtx), extraText=extraText, MCOnly=MCOnly)

   recoilanalyzer.saveHistos(f"root/output_{suffix}.root")
   recoilanalyzerSc.saveHistos(f"root/outputSc_{suffix}.root")