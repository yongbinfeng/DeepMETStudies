import ROOT
import sys, os
import numpy as np
sys.path.append("../RecoilResol/CMSPLOTS")
from CMSPLOTS.myFunction import DrawHistos, getMedian, getResolution, getQuantileResolution
from collections import OrderedDict
from utils.utils import getpTBins, getnVtxBins, get_response_code
import argparse

noLumi = False
MCOnly = True

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)

outdir = f"plots/MC_CompPUPPI/2016"

ROOT.gInterpreter.Declare(get_response_code)

def prepareVars(rdf):
   rdf = rdf.Define("GenHT", "Sum(GenJet_pt)")  \
            .Define("PFMET_pt", "MET_pt") \
            .Define("PFMET_phi", "MET_phi") \
            .Define("PUPPIMET_pt", "PuppiMET_pt") \
            .Define("PUPPIMET_phi", "PuppiMET_phi") \
            .Define("DeepMETMET_pt", "DeepMETResolutionTune_pt") \
            .Define("DeepMETMET_phi", "DeepMETResolutionTune_phi")
   return rdf

def makePlots(fname, suffix, lheader, htbins, ymin, ymax, applySc = True):
   #fdir = "/eos/cms/store/user/yofeng/DeepMETNanoAOD/2016/"
   fdir = "/home/yongbinfeng/Desktop/DeepMET/data/NanoAOD2016/"
   fpath = fdir + fname
   if not os.path.isfile(fpath):
      sys.exit(f"File {fpath} does not exist")
   
   chainMC = ROOT.TChain("Events")
   chainMC.Add(fpath)
   rdf_MC = ROOT.ROOT.RDataFrame(chainMC)
   
   rdf_MC = prepareVars(rdf_MC)

   recoils = ["PF", "PUPPI", "DeepMET"]

   for recoil in recoils:
      rdf_MC = rdf_MC.Define(f"{recoil}MET_diff",  "{RECOIL}MET_pt - GenMET_pt".format(RECOIL=recoil))
      rdf_MC = rdf_MC.Define(f"{recoil}MET_ratio", "{RECOIL}MET_pt / GenMET_pt".format(RECOIL=recoil))
      rdf_MC = rdf_MC.Define(f"{recoil}MET_ratio_diff", "{RECOIL}MET_ratio - 1".format(RECOIL=recoil))
      
   
   nbins_ht = htbins.size - 1
   htmax = htbins[-1]
   htmin = htbins[0]

   histos_2d_diffs = OrderedDict()
   histos_ratios = OrderedDict()
   histos_ratios_diff = OrderedDict()
   for recoil in recoils:
      hname = "h2d_diff_{RECOIL}_vs_HT".format(RECOIL=recoil)
      histos_2d_diffs[recoil] = rdf_MC.Histo2D((hname, hname, nbins_ht, htbins, 400, -200, 200), "GenHT", "{RECOIL}MET_diff".format(RECOIL=recoil))
      hname = "h2d_ratio_{RECOIL}_vs_HT".format(RECOIL=recoil)
      histos_ratios[recoil] = rdf_MC.Histo2D((hname, hname, nbins_ht, htbins, 400, 0.5, 1.5), "GenHT", "{RECOIL}MET_ratio".format(RECOIL=recoil))
      hname = "h2d_ratio_diff_{RECOIL}_vs_HT".format(RECOIL=recoil)
      histos_ratios_diff[recoil] = rdf_MC.Histo2D((hname, hname, nbins_ht, htbins, 400, -0.5, 0.5), "GenHT", "{RECOIL}MET_ratio_diff".format(RECOIL=recoil))

   hresols = OrderedDict()
   hresolsRMS = OrderedDict()
   hresols_ratio = OrderedDict()
   hresponses = OrderedDict()
   for recoil in recoils:
      hresponses[recoil] = getMedian(histos_ratios[recoil]) 
      hresols[recoil] = getQuantileResolution(histos_2d_diffs[recoil])
      hresolsRMS[recoil] = getResolution(histos_2d_diffs[recoil], useRMS=True)
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

   # loop over the different qT bins
   if applySc:
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
         histos_2d_diffs_sc[recoil] = rdf_MC.Histo2D((hname, hname, nbins_ht, htbins, 400, -200, 200), "GenHT", "{RECOIL}MET_Sc_diff".format(RECOIL=recoil))
         hname = "h2d_ratio_{RECOIL}Sc_vs_HT".format(RECOIL=recoil)
         histos_ratios_sc[recoil] = rdf_MC.Histo2D((hname, hname, nbins_ht, htbins, 400, 0.5, 1.5), "GenHT", "{RECOIL}MET_Sc_ratio".format(RECOIL=recoil))
         hname = "h2d_ratio_diff_{RECOIL}Sc_vs_HT".format(RECOIL=recoil)
         histos_ratios_diff_sc[recoil] = rdf_MC.Histo2D((hname, hname, nbins_ht, htbins, 400, -0.5, 0.5), "GenHT", "{RECOIL}MET_Sc_ratio_diff".format(RECOIL=recoil))

      hresponses_sc = OrderedDict()
      hresols_sc = OrderedDict()
      hresolsRMS_sc = OrderedDict()
      hresols_ratio_sc = OrderedDict()
      for recoil in recoils:
         hresponses_sc[recoil] = getMedian(histos_ratios_sc[recoil]) 
         hresols_sc[recoil] = getQuantileResolution(histos_2d_diffs_sc[recoil])
         hresolsRMS_sc[recoil] = getResolution(histos_2d_diffs_sc[recoil], useRMS=True)
         hresols_ratio_sc[recoil] = getResolution(histos_ratios_diff_sc[recoil])

   def GetLegends(hdict):
       return [labels[itype] for itype in hdict.keys() if "_MC" not in itype]

   def GetLineStyles(hdict):
       return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]

   DrawHistos(hresponses.values(), GetLegends(hresponses), htmin, htmax, "Gen H_{T} [GeV]", 0., 1.20, "Response -<p^{miss}_{T}>/<Gen p^{miss}_{T}>", "met_response" + suffix, drawashist=True, dology=False, legendPos=[0.50, 0.20, 0.80, 0.40], mycolors=[colors[itype] for itype in hresponses.keys()], linestyles = GetLineStyles(hresponses), noLumi=noLumi, outdir=outdir, MCOnly=MCOnly, lheader=lheader)

   DrawHistos(hresols.values(), GetLegends(hresols), htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols), MCOnly=MCOnly, lheader=lheader)
   
   DrawHistos(hresolsRMS.values(), GetLegends(hresolsRMS), htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol_RMS" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresolsRMS.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresolsRMS), MCOnly=MCOnly, lheader=lheader)
   
   # difference between quantile and RMS resol
   DrawHistos([hresols["DeepMET"], hresolsRMS["DeepMET"]], ["DeepMET", "DeepMET RMS"], htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol_RMSVsQuantitle_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, lheader=lheader, showratio=True, yrlabel="RMS/Quantile", yrmin=0.91, yrmax=1.09)
   

   DrawHistos(hresols_ratio.values(), GetLegends(hresols_ratio), htmin, htmax, "Gen H_{T} [GeV]", 0, 0.5, "#sigma (p^{miss}_{T}) / <p^{miss}_{T}>", "met_resol_ratio" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_ratio.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_ratio), MCOnly=MCOnly, lheader=lheader)


   if applySc:
      
      extraToDraw = ROOT.TPaveText(0.60, 0.77, 0.90, 0.82, "NDC")
      extraToDraw.SetFillColorAlpha(0, 0)
      extraToDraw.SetBorderSize(0)
      extraToDraw.SetTextFont(42)
      extraToDraw.SetTextSize(0.04)
      extraToDraw.AddText("Response corrected")
      #
      # Scaled 
      #
      DrawHistos(hresponses_sc.values(), GetLegends(hresponses), htmin, htmax, "Gen H_{T} [GeV]", 0., 1.20, "Scaled Response -<p^{miss}_{T}>/<Gen p^{miss}_{T}>", "met_response_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.70, 0.17, 0.88, 0.36], mycolors=[colors[itype] for itype in hresponses.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresponses), MCOnly=MCOnly, lheader=lheader)

      DrawHistos(hresols_sc.values(), GetLegends(hresols), htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols), MCOnly=MCOnly, lheader=lheader, extraToDraw=extraToDraw)
      
      DrawHistos(hresolsRMS_sc.values(), GetLegends(hresolsRMS), htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol_RMS_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresolsRMS.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresolsRMS), MCOnly=MCOnly, lheader=lheader, extraToDraw=extraToDraw)
      
      DrawHistos([hresols_sc["DeepMET"], hresolsRMS_sc["DeepMET"]], ["DeepMET", "DeepMET RMS"], htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol_Scaled_RMSVsQuantitle_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, lheader=lheader, showratio=True, yrlabel="RMS/Quantile", yrmin=0.91, yrmax=1.09)

      DrawHistos(hresols_ratio_sc.values(), GetLegends(hresols_ratio), htmin, htmax, "Gen H_{T} [GeV]", 0, 0.5, "Response-corrected #sigma (p^{miss}_{T}) / <p^{miss}_{T}>", "met_resol_ratio_Scaled" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors[itype] for itype in hresols_ratio.keys()], noLumi=noLumi, outdir=outdir, linestyles = GetLineStyles(hresols_ratio), MCOnly=MCOnly, lheader=lheader)

      DrawHistos([hresponses["DeepMET"], hresponses_sc["DeepMET"]], ["DeepMET", "DeepMET Scaled"], htmin, htmax, "Gen H_{T} [GeV]", 0., 1.20, "Response -<p^{miss}_{T}>/<Gen p^{miss}_{T}>", "met_response_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.50, 0.20, 0.80, 0.40], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, showratio=True, yrlabel="Scaled/Unscaled", yrmin=0.91, yrmax=1.09, lheader=lheader)

      DrawHistos([hresols["DeepMET"], hresols_sc["DeepMET"]], ["DeepMET", "DeepMET Scaled"], htmin, htmax, "Gen H_{T} [GeV]", ymin, ymax, "#sigma (p^{miss}_{T}) [GeV]", "met_resol_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, showratio=True, yrlabel="Scaled/Unscaled", yrmin=0.91, yrmax=1.09, lheader=lheader)

      DrawHistos([hresols_ratio["DeepMET"], hresols_ratio_sc["DeepMET"]], ["DeepMET", "DeepMET Scaled"], htmin, htmax, "Gen H_{T} [GeV]", 0, 0.5, "#sigma (p^{miss}_{T}) / <p^{miss}_{T}>", "met_resol_ratio_DeepMET" + suffix, drawashist=True, dology=False, legendPos=[0.20, 0.69, 0.38, 0.88], mycolors=[colors["DeepMET"], colors["DeepMET"]], noLumi=noLumi, outdir=outdir, linestyles = [1, 2], MCOnly=MCOnly, showratio=True, yrlabel="Scaled/Unscaled", yrmin=0.91, yrmax=1.09, lheader=lheader)
   
if __name__ == "__main__":
   
   parser = argparse.ArgumentParser()
   parser.add_argument("--applySc", action="store_true", help="Apply response corrections and include the plots in the output")
   parser.set_defaults(applySc=True)
   args = parser.parse_args()
   
   fnames = [
      ("GluGluToHHTo2B2Tau_TuneCP5_PSWeights_node_SM_13TeV-madgraph-pythia8.root", "GGToHHTo2B2Tau", "HH#rightarrowb#bar{b}#tau^{+}#tau^{-}", np.array([0.0, 200.0, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000]), 10, 45.0),
      ("SMS-T2tt-4bd_genMET-100_genHT200_mStop-350_mLSP-335_TuneCP5_LLStop_13TeV-madgraphMLM-pythia8.root", "T2tt_4bd", "SMS T2tt 4bd", np.array([200.0, 250, 300, 350.0, 400., 450., 500, 600, 700, 800, 1000]), 15.0, 50.0),
      ("SMS-TChiZZ-mNLSP400_mLSP1_TuneCP5_13TeV-madgraphMLM-pythia8.root", "TChiZZ", "SMS TChiZZ", np.array([0, 100.0, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 1000]), 5.0, 45.0),
      ("ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8.root", "ttHTobb", "ttH(H#rightarrowb#bar{b}) dilepton", np.array([100.0, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000]), 10.0, 45.0),
      ("ttHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8.root", "ttHToMuMu", "ttH(H#rightarrow#mu^{+}#mu^{-})", np.array([200.0, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000]), 10.0, 45.0),
      ("VBF_HToInvisible_M125_TuneCP5_withDipoleRecoil_13TeV_powheg_pythia8.root", "VBF_HToInvisible", "VBF H#rightarrow invisible", np.array([0.0, 100.0, 240, 280, 320, 360, 400, 450, 500, 600]), 10.0, 50.0),
   ]
   for fname, suffix, lheader, htbins, ymin, ymax in fnames:
      makePlots(fname, "_" + suffix, lheader, htbins, ymin, ymax, args.applySc)
   
