'''
plot the data-MC comparisons pre/post DeepMET corrections.
and the systematic uncs.
'''
import ROOT
import numpy as np
from collections import OrderedDict
import sys
sys.path.append("../RecoilResol/CMSPLOTS")
from tdrstyle import setTDRStyle
from myFunction import DrawHistos, THStack2TH1
import CMS_lumi
import pickle

from SampleManager import DrawConfig, Sample, SampleManager

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(15)

dotest = 0
VetoB = False

def main():
    print("Program start...")

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

    if dotest:
        input_data    = "inputs/inputs_Z_UL/input_data_test.txt"
        input_dy      = "inputs/inputs_Z_UL/input_zjets.txt"
        input_ttbar   = "inputs/inputs_Z_UL/input_ttbar_test.txt"
    else:
        input_data    = "inputs/inputs_Z_UL/input_data.txt"
        input_dy      = "inputs/inputs_Z_UL/input_zjets_all.txt"
        input_ttbar   = "inputs/inputs_Z_UL/input_ttbar.txt"
        input_WW2L    = "inputs/inputs_Z_UL/input_WWTo2L2Nu.txt"
        input_WZ2L    = "inputs/inputs_Z_UL/input_WZTo2Q2L.txt"
        input_ZZ2L    = "inputs/inputs_Z_UL/input_ZZTo2L2Nu.txt"
        input_ZZ2L2Q  = "inputs/inputs_Z_UL/input_ZZTo2Q2L.txt"
        input_dytau   = "inputs/inputs_Z_UL/input_zjets_tautau.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", bjetVeto=VetoB)
    DYSamp    = Sample(input_dy,    xsec = 2025.74*1e3,      color=5,  reweightzpt = False, legend="DY", name="DY", bjetVeto=VetoB)
    TTbarSamp = Sample(input_ttbar, xsec = 831.76*0.105*1e3, color=46, reweightzpt = False, legend="t#bar{t}", name="ttbar", bjetVeto=VetoB)
    if not dotest:
        WW2LSamp  = Sample(input_WW2L,  xsec = 12.178*1e3,       color=38, reweightzpt = False, legend="WW2L",     name="WW2L",  bjetVeto=VetoB)
        WZ2LSamp  = Sample(input_WZ2L,  xsec = 5.26*1e3,         color=39, reweightzpt = False, legend="WZ2L",     name="WZ2L",  bjetVeto=VetoB)
        ZZ2LSamp  = Sample(input_ZZ2L,  xsec = 0.564*1e3,        color=37, reweightzpt = False, legend="ZZ2L",     name="ZZ2L",  bjetVeto=VetoB)
        ZZ2L2QSamp  = Sample(input_ZZ2L2Q,  xsec = 1.212*1e3,        color=36, reweightzpt = False, legend="ZZ2L2Q",     name="ZZ2L2Q",  bjetVeto=VetoB) 
        DYTauSamp   = Sample(input_dytau, xsec = 2025.74*1e3,      color=8,  reweightzpt = False, legend="DY#rightarrow#tau#tau", name="DYTauTau", bjetVeto=VetoB)

    if not dotest:
        sampMan = SampleManager(DataSamp, [DYSamp, DYTauSamp, WW2LSamp, WZ2LSamp, TTbarSamp, ZZ2LSamp, ZZ2L2QSamp])
    else:
        sampMan = SampleManager(DataSamp, [DYSamp, TTbarSamp])
    sampMan.groupMCs(["WW2L", "WZ2L", "ZZ2L", "ZZ2L2Q"], "Dibosons", 38, "Dibosons")

    sampMan.DefineAll("zpt", "Z_pt")

    # define u1 and u2 for PUPPI and PF
    # u1 and u2 for DeepMET has already been defined in SampleManager
    #sampMan.DefineAll("Uvec_PF", "UVec(Z_pt, Z_phi, met_pt, met_phi)") 
    #sampMan.DefineAll("u_PF_pt",  "Uvec_PF.Mod()") 
    #sampMan.DefineAll("u_PF_phi", "Uvec_PF.Phi()")
    #sampMan.DefineAll("u1_PF",    "u_PF_pt * TMath::Cos(u_PF_phi + TMath::Pi() - Z_phi)")
    #sampMan.DefineAll("u2_PF",    "u_PF_pt * TMath::Sin(u_PF_phi + TMath::Pi() - Z_phi)")

    #sampMan.DefineAll("Uvec_PUPPI", "UVec(Z_pt, Z_phi, puppimet_pt, puppimet_phi)")
    #sampMan.DefineAll("u_PUPPI_pt",  "Uvec_PUPPI.Mod()")
    #sampMan.DefineAll("u_PUPPI_phi", "Uvec_PUPPI.Phi()")
    #sampMan.DefineAll("u1_PUPPI",    "u_PUPPI_pt * TMath::Cos(u_PUPPI_phi + TMath::Pi() - Z_phi)")
    #sampMan.DefineAll("u2_PUPPI",    "u_PUPPI_pt * TMath::Sin(u_PUPPI_phi + TMath::Pi() - Z_phi)") 

    # Default correction:
    ROOT.gROOT.ProcessLine('TFile* fitfunctions_Data = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TH1F* h1_ptbins_Data_central   = (TH1F*)fitfunctions_Data->Get("h1_ptbins_Data_central")')
    ROOT.gROOT.ProcessLine('TH1F* h1_njetbins_Data_central = (TH1F*)fitfunctions_Data->Get("h1_njetbins_Data_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_Data_u1_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_Data_u2_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u2_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TFile* fitfunctions_ZJets_NLO_central = TFile::Open("results/Fit/fitfunctions_ZJets_NLO_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TList* tfs_MC_u1_njets_pt_central = (TList*)fitfunctions_ZJets_NLO_central->Get("tfs_ZJets_NLO_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_MC_u2_njets_pt_central = (TList*)fitfunctions_ZJets_NLO_central->Get("tfs_ZJets_NLO_u2_njets_pt_central")')

    DYSamp.Define("u1_corr_central",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_central, tfs_MC_u1_njets_pt_central, 0.00001)")
    DYSamp.Define("u2_corr_central",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_central, tfs_MC_u2_njets_pt_central, 0.00001)")
    DYSamp.Define("u_pt_corr_central", "TMath::Sqrt(u1_corr_central*u1_corr_central + u2_corr_central*u2_corr_central)")
    sampMan.DefineAll("u1_corr_central",     "u1"  , excludes=['DY'])
    sampMan.DefineAll("u2_corr_central",     "u2"  , excludes=['DY'])
    sampMan.DefineAll("u_pt_corr_central",   "u_pt", excludes=['DY'])
    sampMan.DefineAll("deepmet_corr_central", "METVec(Z_pt, Z_phi, u1_corr_central, u2_corr_central)")
    sampMan.DefineAll("deepmet_pt_corr_central", "deepmet_corr_central.Mod()") 
    sampMan.DefineAll("deepmet_phi_corr_central", "TVector2::Phi_mpi_pi(deepmet_corr_central.Phi())")


    # systematic variations
    # corrections from MadGraph
    #ROOT.gROOT.ProcessLine('TFile* fitfunctions_ZJets_MG_central = TFile::Open("results/Fit/fitfunctions_ZJets_MG_njets_pt_central.root")')
    #ROOT.gROOT.ProcessLine('TList* tfs_MC_MG_u1_njets_pt_central = (TList*)fitfunctions_ZJets_MG_central->Get("tfs_ZJets_MG_u1_njets_pt_central")')
    #ROOT.gROOT.ProcessLine('TList* tfs_MC_MG_u2_njets_pt_central = (TList*)fitfunctions_ZJets_MG_central->Get("tfs_ZJets_MG_u2_njets_pt_central")')

    #DYSamp.Define("u1_corr_central_MG",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_central, tfs_MC_MG_u1_njets_pt_central, 0.00001)")
    #DYSamp.Define("u2_corr_central_MG",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_central, tfs_MC_MG_u2_njets_pt_central, 0.00001)")
    #DYSamp.Define("u_pt_corr_central_MG", "TMath::Sqrt(u1_corr_central_MG*u1_corr_central_MG + u2_corr_central_MG*u2_corr_central_MG)")
    #sampMan.DefineAll("u1_corr_central_MG",     "u1"  , excludes=['DY'])
    #sampMan.DefineAll("u2_corr_central_MG",     "u2"  , excludes=['DY'])
    #sampMan.DefineAll("u_pt_corr_central_MG",   "u_pt", excludes=['DY'])
    #sampMan.DefineAll("deepmet_corr_central_MG", "METVec(Z_pt, Z_phi, u1_corr_central_MG, u2_corr_central_MG)")
    #sampMan.DefineAll("deepmet_pt_corr_central_MG", "deepmet_corr_central_MG.Mod()")
    #sampMan.DefineAll("deepmet_phi_corr_central_MG", "TVector2::Phi_mpi_pi(deepmet_corr_central_MG.Phi())")

    ## corrections from Gaussian Smearing
    #ROOT.gROOT.ProcessLine('TFile* gaussSmoother_Data_njets_pt_central = TFile::Open("results/GaussSmoother/gaussSmoother_Data_njets_pt_central.root")')
    #ROOT.gROOT.ProcessLine('TList* cdfs_Data_u1_njets_pt_central = (TList*)gaussSmoother_Data_njets_pt_central->Get("cdfs_Data_u1_njets_pt_central")')
    #ROOT.gROOT.ProcessLine('TList* cdfs_Data_u2_njets_pt_central = (TList*)gaussSmoother_Data_njets_pt_central->Get("cdfs_Data_u2_njets_pt_central")')
    #ROOT.gROOT.ProcessLine('TFile* gaussSmoother_ZJets_NLO_njets_pt_central = TFile::Open("results/GaussSmoother/gaussSmoother_ZJets_NLO_njets_pt_central.root")')
    #ROOT.gROOT.ProcessLine('TList* cdfs_ZJets_NLO_u1_njets_pt_central = (TList*)gaussSmoother_ZJets_NLO_njets_pt_central->Get("cdfs_ZJets_NLO_u1_njets_pt_central")')
    #ROOT.gROOT.ProcessLine('TList* cdfs_ZJets_NLO_u2_njets_pt_central = (TList*)gaussSmoother_ZJets_NLO_njets_pt_central->Get("cdfs_ZJets_NLO_u2_njets_pt_central")')

    #DYSamp.Define("u1_corr_central_GKS",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, cdfs_Data_u1_njets_pt_central, cdfs_ZJets_NLO_u1_njets_pt_central, 0.00001)")
    #DYSamp.Define("u2_corr_central_GKS",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, cdfs_Data_u2_njets_pt_central, cdfs_ZJets_NLO_u2_njets_pt_central, 0.00001)")
    #DYSamp.Define("u_pt_corr_central_GKS", "TMath::Sqrt(u1_corr_central_GKS*u1_corr_central_GKS + u2_corr_central_GKS*u2_corr_central_GKS)")
    #sampMan.DefineAll("u1_corr_central_GKS",     "u1"  , excludes=['DY'])
    #sampMan.DefineAll("u2_corr_central_GKS",     "u2"  , excludes=['DY'])
    #sampMan.DefineAll("u_pt_corr_central_GKS",   "u_pt", excludes=['DY'])
    #sampMan.DefineAll("deepmet_corr_central_GKS", "METVec(Z_pt, Z_phi, u1_corr_central_GKS, u2_corr_central_GKS)")
    #sampMan.DefineAll("deepmet_pt_corr_central_GKS", "deepmet_corr_central_GKS.Mod()") 
    #sampMan.DefineAll("deepmet_phi_corr_central_GKS", "TVector2::Phi_mpi_pi(deepmet_corr_central_GKS.Phi())")

    # correction from ttbar qcd scale:
    def varyQCDScale(wstring):
        ROOT.gROOT.ProcessLine('TFile* fitfunctions_Data_{WS} = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_{WS}.root")'.format(WS=wstring))
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u1_njets_pt_{WS}  = (TList*)fitfunctions_Data_{WS}->Get("tfs_Data_u1_njets_pt_{WS}")'.format(WS=wstring))
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u2_njets_pt_{WS}  = (TList*)fitfunctions_Data_{WS}->Get("tfs_Data_u2_njets_pt_{WS}")'.format(WS=wstring))

        DYSamp.Define("u1_corr_{}".format(wstring),   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_{WS}, tfs_MC_u1_njets_pt_central, 0.00001)".format(WS=wstring))
        DYSamp.Define("u2_corr_{}".format(wstring),   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_{WS}, tfs_MC_u2_njets_pt_central, 0.00001)".format(WS=wstring))
        DYSamp.Define("u_pt_corr_{}".format(wstring), "TMath::Sqrt(u1_corr_{WS}*u1_corr_{WS} + u2_corr_{WS}*u2_corr_{WS})".format(WS=wstring))
        sampMan.DefineAll("u1_corr_{}".format(wstring),     "u1"  , excludes=['DY'])
        sampMan.DefineAll("u2_corr_{}".format(wstring),     "u2"  , excludes=['DY'])
        sampMan.DefineAll("u_pt_corr_{}".format(wstring),   "u_pt", excludes=['DY'])
        sampMan.DefineAll("deepmet_corr_{}".format(wstring), "METVec(Z_pt, Z_phi, u1_corr_{WS}, u2_corr_{WS})".format(WS=wstring))
        sampMan.DefineAll("deepmet_pt_corr_{}".format(wstring), "deepmet_corr_{}.Mod()".format(wstring))
        sampMan.DefineAll("deepmet_phi_corr_{}".format(wstring), "TVector2::Phi_mpi_pi(deepmet_corr_{}.Phi())".format(wstring))

    #varyQCDScale("mur_1_muf_0p5")
    #varyQCDScale("mur_1_muf_2")
    #varyQCDScale("mur_0p5_muf_1")
    #varyQCDScale("mur_2_muf_1")
    #varyQCDScale("mur_0p5_muf_0p5")
    #varyQCDScale("mur_2_muf_2")

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150, 165, 180, 200])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    # z pt befor and after pt reweighting
    sampMan.cacheDraw("zpt", "histo_zjets_zpt_WoZptWeight", 30, 0, 60, DrawConfig(xmin=0, xmax=60, xlabel='p^{ll}_{T} [GeV]'), weightname="weight_WoVpt")
    sampMan.cacheDraw("zpt", "histo_zjets_zpt", 30, 0, 60, DrawConfig(xmin=0, xmax=60, xlabel='p^{ll}_{T} [GeV]'))

    sampMan.cacheDraw("MET_pt", "histo_zjets_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='PF MET [GeV]'))
    sampMan.cacheDraw("MET_phi", "histo_zjets_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF #phi', ymax=1e10))
    sampMan.cacheDraw("PuppiMET_pt", "histo_zjets_puppimet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='PUPPI MET [GeV]'))
    sampMan.cacheDraw("PuppiMET_phi", "histo_zjets_puppimet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PUPPI #phi', ymax=1e10))

    # uncorrected deepmet
    sampMan.cacheDraw("DeepMETResolutionTune_pt", "histo_zjets_deepmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET [GeV]'))
    sampMan.cacheDraw("DeepMETResolutionTune_phi", "histo_zjets_deepmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET #phi', ymax=1e10))
    sampMan.cacheDraw("u1", "histo_zjets_u1", u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} [GeV]'))
    sampMan.cacheDraw("u2", "histo_zjets_u2", u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} [GeV]'))
    sampMan.cacheDraw("u_pt",    "histo_zjets_u_pt", u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u [GeV]'))

    # corrected deepmet
    def DrawCorrection(postfix):
        sampMan.cacheDraw("deepmet_pt_corr_"+postfix, "histo_zjets_deepmet_pt_corr_"+postfix, met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET Corrected [GeV]'))
        sampMan.cacheDraw("deepmet_phi_corr_"+postfix, "histo_zjets_deepmet_phi_corr_"+postfix, 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET Corrected #phi'))
        sampMan.cacheDraw("u1_corr_"+postfix, "histo_zjets_u1_corr_"+postfix, u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]'))
        sampMan.cacheDraw("u2_corr_"+postfix, "histo_zjets_u2_corr_"+postfix, u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} Corrected [GeV]'))
        sampMan.cacheDraw("u_pt_corr_"+postfix,    "histo_zjets_u_pt_corr_"+postfix, u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u Corrected [GeV]'))

    DrawCorrection("central")
    #DrawCorrection("central_MG")
    #DrawCorrection("central_GKS")
    #DrawCorrection("mur_1_muf_0p5")
    #DrawCorrection("mur_1_muf_2")
    #DrawCorrection("mur_0p5_muf_1")
    #DrawCorrection("mur_2_muf_1")
    #DrawCorrection("mur_0p5_muf_0p5")
    #DrawCorrection("mur_2_muf_2")

    ## sys: corrected deepmet from MG
    #sampMan.cacheDraw("deepmet_pt_corr_central_MG", "histo_zjets_deepmet_pt_corr_central_MG", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET Corrected [GeV]'))
    #sampMan.cacheDraw("u1_corr_central_MG", "histo_zjets_u1_corr_central_MG", u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]'))
    #sampMan.cacheDraw("u2_corr_central_MG", "histo_zjets_u2_corr_central_MG", u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} Corrected [GeV]'))
    #sampMan.cacheDraw("u_pt_corr_central_MG",    "histo_zjets_u_pt_corr_central_MG", u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u Corrected [GeV]')) 

    ## sys: corrected deepmet from Gauss Smearing
    #sampMan.cacheDraw("deepmet_pt_corr_central_GKS", "histo_zjets_deepmet_pt_corr_central_GKS", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET Corrected [GeV]'))
    #sampMan.cacheDraw("u1_corr_central_GKS", "histo_zjets_u1_corr_central_GKS", u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]'))
    #sampMan.cacheDraw("u2_corr_central_GKS", "histo_zjets_u2_corr_central_GKS", u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} Corrected [GeV]'))
    #sampMan.cacheDraw("u_pt_corr_central_GKS",    "histo_zjets_u_pt_corr_central_GKS", u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u Corrected [GeV]'))  

    sampMan.launchDraw()

    # compare some ratio distributions for systematics
    rmin = 0.7
    rmax = 1.3
    colors = [1,2, 3]
    linestyles = [1,3,3]

    def compRatios(hnames, colors, linestyles, legends, outputtag):
        histos_met =  [list(sampMan.hratios["histo_zjets_deepmet_pt_{}".format(hname)].values())[0] for hname in hnames]
        histos_u_pt = [list(sampMan.hratios["histo_zjets_u_pt_{}".format(hname)      ].values())[0] for hname in hnames]
        histos_u1   = [list(sampMan.hratios["histo_zjets_u1_{}".format(hname)        ].values())[0] for hname in hnames]
        histos_u2   = [list(sampMan.hratios["histo_zjets_u2_{}".format(hname)        ].values())[0] for hname in hnames]

        if len(hnames)>5:
            legendNCols = 2
            legendPos=[0.20, 0.88, 0.70, 0.74]
        else:
            legendNCols = 1
            legendPos=[0.20, 0.88, 0.35, 0.74]
            
        DrawHistos(histos_met, legends, 0., 150., "Deep MET Corrected [GeV]", rmin, rmax, "Data / MC", "histo_zjets_ratio_deepmet_pt_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u_pt, legends, 0., 150., "u Corrected [GeV]", 0.7, 1.3, "Data / MC", "histo_zjets_ratio_u_pt_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u1, legends, -40.0, 100., "u_{#parallel} Corrected [GeV]", 0.7, 1.3, "Data / MC", "histo_zjets_ratio_u1_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u2, legends, -80.0, 80., "u_{#perp} Corrected [GeV]", 0.7, 1.3, "Data / MC", "histo_zjets_ratio_u2_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols) 

    def compMCRatios(hnames, colors, linestyles, legends, outputtag):
        # compare MC with MC
        histos_met =  [THStack2TH1(sampMan.hsmcs["histo_zjets_deepmet_pt_{}".format(hname)], "MCRatio") for hname in hnames]
        histos_u_pt = [THStack2TH1(sampMan.hsmcs["histo_zjets_u_pt_{}".format(hname)      ], "MCRatio") for hname in hnames]
        histos_u1   = [THStack2TH1(sampMan.hsmcs["histo_zjets_u1_{}".format(hname)        ], "MCRatio") for hname in hnames]
        histos_u2   = [THStack2TH1(sampMan.hsmcs["histo_zjets_u2_{}".format(hname)        ], "MCRatio") for hname in hnames] 

        histos_met_base = histos_met[0].Clone(histos_met[0].GetName()+"_Cloned")
        histos_u_pt_base = histos_u_pt[0].Clone(histos_u_pt[0].GetName()+"_Cloned")
        histos_u1_base = histos_u1[0].Clone(histos_u1[0].GetName()+"_Cloned")
        histos_u2_base = histos_u2[0].Clone(histos_u2[0].GetName()+"_Cloned")

        for hmet, hupt, hu1, hu2 in zip(histos_met, histos_u_pt, histos_u1, histos_u2):
            hmet.Divide(  histos_met_base )
            hmet.SetFillColor(0)
            hupt.Divide( histos_u_pt_base )
            hupt.SetFillColor(0)
            hu1.Divide(  histos_u1_base )
            hu1.SetFillColor(0)
            hu2.Divide(  histos_u2_base )
            hu2.SetFillColor(0)

        if len(hnames)>5:
            legendNCols = 2
            legendPos=[0.20, 0.88, 0.70, 0.74]
        else:
            legendNCols = 1
            legendPos=[0.20, 0.88, 0.35, 0.74]

        DrawHistos(histos_met, legends, 0., 150., "Deep MET Corrected [GeV]", rmin, rmax, "Ratio", "histo_zjets_MCratio_deepmet_pt_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u_pt, legends, 0., 150., "u Corrected [GeV]", 0.7, 1.3, "Ratio", "histo_zjets_MCratio_u_pt_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u1, legends, -40.0, 100., "u_{#parallel} Corrected [GeV]", 0.7, 1.3, "Ratio", "histo_zjets_MCratio_u1_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u2, legends, -80.0, 80., "u_{#perp} Corrected [GeV]", 0.7, 1.3, "Ratio", "histo_zjets_MCratio_u2_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)

    colors = [1,2]
    linestyles = [1,1]
    legends = ["amc@NLO", "MadGraph"]
    #compRatios(  ["corr_central", "corr_central_MG"], colors, linestyles, legends, "MCSample")
    #compMCRatios(["corr_central", "corr_central_MG"], colors, linestyles, legends, "MCSample")

    legends = ["2-Gaussian Fit", "Gaussian Smoothing"]
    #compRatios(  ["corr_central", "corr_central_GKS"], colors, linestyles, legends, "FitVSSmoothing")
    #compMCRatios(["corr_central", "corr_central_GKS"], colors, linestyles, legends, "FitVSSmoothing")

    legends = ["#mu_{r}=1, #mu_{f}=1", "#mu_{r}=1, #mu_{f}=0.5", "#mu_{r}=1, #mu_{f}=2", "#mu_{r}=0.5, #mu_{f}=1", "#mu_{r}=2, #mu_{f}=1", "#mu_{r}=0.5, #mu_{f}=0.5", "#mu_{r}=2, #mu_{f}=2"]
    postfixs = ["corr_central", "corr_mur_1_muf_0p5", "corr_mur_1_muf_2", "corr_mur_0p5_muf_1", "corr_mur_2_muf_1", "corr_mur_0p5_muf_0p5", "corr_mur_2_muf_2"]
    colors     = [1,2,3,4,5,6,7]
    linestyles = [1,1,1,1,1,1,1]
    #compRatios(  postfixs, colors, linestyles, legends, "QCDScale")
    #compMCRatios(postfixs, colors, linestyles, legends, "QCSScale") 


    doDump = False
    if doDump:
        print("Dump corrections into root file...")
        branchList = ROOT.vector('string')()
        for branchName in ["zpt", "deepmet_pt", "u1", "u2", "u1_corr_njets_3Gauss", "u2_corr_njets_3Gauss", "u_pt_corr_njets_3Gauss", "deepmet_pt_corr_njets_3Gauss","weight", "weight_WoVpt", "vtx_n", "jet_CSVLoose_n"]:
            branchList.push_back(branchName)
        #TTbarSamp.rdf.Snapshot("DeepMETCorr_TTbar", "results/deepmetcorr_TTbar.root", branchList)
        DYSamp.rdf.Snapshot("DeepMETCorr_DY",       "results/deepmetcorr_DY.root", branchList)

        branchList = ROOT.vector('string')()
        for branchName in ["zpt", "deepmet_pt", "u1", "u2", "weight", "weight_WoVpt", "vtx_n", "jet_CSVLoose_n"]:
            branchList.push_back(branchName)
        DataSamp.rdf.Snapshot("DeepMETCorr_Data",   "results/deepmetcorr_Data.root", branchList)


    print("Program end...")

    input()
    
    return 

if __name__ == "__main__":
   main()
