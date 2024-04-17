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

doDefaultCorrection = True
doGaussianSmooth = True
doBkgScaled = True

def main():
    print("Program start...")

    input_data    = "inputs/inputs_Z_UL_Post/input_data.txt"
    input_dy      = "inputs/inputs_Z_UL_Post/input_zjets.txt"
    input_ttbar   = "inputs/inputs_Z_UL_Post/input_ttbar.txt"
    input_WW2L    = "inputs/inputs_Z_UL_Post/input_WWTo2L2Nu.txt"
    input_WZ2L    = "inputs/inputs_Z_UL_Post/input_WZTo2Q2L.txt"
    input_ZZ2L    = "inputs/inputs_Z_UL_Post/input_ZZTo2L2Nu.txt"
    input_ZZ2L2Q  = "inputs/inputs_Z_UL_Post/input_ZZTo2Q2L.txt"
    input_dytau   = "inputs/inputs_Z_UL_Post/input_zjets_tautau.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", prepareVars=False, select=False)
    DYSamp    = Sample(input_dy,    xsec = 0,  color=5,  reweightzpt = False, legend="DY", name="DY", prepareVars=False, select=False)
    TTbarSamp = Sample(input_ttbar, xsec = 0,  color=46, reweightzpt = False, legend="t#bar{t}", name="ttbar", prepareVars=False, select=False)
    WW2LSamp  = Sample(input_WW2L,  xsec = 0,  color=38, reweightzpt = False, legend="WW2L",     name="WW2L",  prepareVars=False, select=False)
    WZ2LSamp  = Sample(input_WZ2L,  xsec = 0,  color=39, reweightzpt = False, legend="WZ2L",     name="WZ2L",  prepareVars=False, select=False)
    ZZ2LSamp  = Sample(input_ZZ2L,  xsec = 0,  color=37, reweightzpt = False, legend="ZZ2L",     name="ZZ2L",  prepareVars=False, select=False)
    ZZ2L2QSamp  = Sample(input_ZZ2L2Q,  xsec = 0, color=36, reweightzpt = False, legend="ZZ2L2Q", name="ZZ2L2Q",  prepareVars=False, select=False) 
    DYTauSamp   = Sample(input_dytau, xsec = 0, color=8,  reweightzpt = False, legend="DY#rightarrow#tau#tau", name="DYTauTau", prepareVars=False, select=False)

    sampMan = SampleManager(DataSamp, [DYSamp, DYTauSamp, WW2LSamp, WZ2LSamp, TTbarSamp, ZZ2LSamp, ZZ2L2QSamp])
    sampMan.groupMCs(["WW2L", "WZ2L", "ZZ2L", "ZZ2L2Q"], "Dibosons", 38, "Dibosons")
    
    DataSamp.Define("norm", "1")
    sampMan.DefineAll("weight_corr", "weight * norm")
    
    sampMan.SetDefaultWeightName("weight_corr")
    
    if doDefaultCorrection:
        # Default correction:
        ROOT.gROOT.ProcessLine('TFile* fitfunctions_Data = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_central.root")')
        ROOT.gROOT.ProcessLine('TH1F* h1_ptbins_Data_central   = (TH1F*)fitfunctions_Data->Get("h1_ptbins_Data_central")')
        ROOT.gROOT.ProcessLine('TH1F* h1_njetbins_Data_central = (TH1F*)fitfunctions_Data->Get("h1_njetbins_Data_central")')
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u1_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u1_njets_pt_central")')
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u2_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u2_njets_pt_central")')
        ROOT.gROOT.ProcessLine('TFile* fitfunctions_DY_central = TFile::Open("results/Fit/fitfunctions_DY_njets_pt_central.root")')
        ROOT.gROOT.ProcessLine('TList* tfs_MC_u1_njets_pt_central = (TList*)fitfunctions_DY_central->Get("tfs_DY_u1_njets_pt_central")')
        ROOT.gROOT.ProcessLine('TList* tfs_MC_u2_njets_pt_central = (TList*)fitfunctions_DY_central->Get("tfs_DY_u2_njets_pt_central")')

        DYSamp.Define("u1_corr_central",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_central, tfs_MC_u1_njets_pt_central, 0.00001)")
        DYSamp.Define("u2_corr_central",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_central, tfs_MC_u2_njets_pt_central, 0.00001)")
        DYSamp.Define("u_pt_corr_central", "TMath::Sqrt(u1_corr_central*u1_corr_central + u2_corr_central*u2_corr_central)")
        sampMan.DefineAll("u1_corr_central",     "u1"  , excludes=['DY'])
        sampMan.DefineAll("u2_corr_central",     "u2"  , excludes=['DY'])
        sampMan.DefineAll("u_pt_corr_central",   "u_pt", excludes=['DY'])
        sampMan.DefineAll("deepmet_corr_central", "METVec(Z_pt, Z_phi, u1_corr_central, u2_corr_central)")
        sampMan.DefineAll("deepmet_pt_corr_central", "deepmet_corr_central.Mod()") 
        sampMan.DefineAll("deepmet_phi_corr_central", "TVector2::Phi_mpi_pi(deepmet_corr_central.Phi())")
        
    if doGaussianSmooth:
        ## corrections from Gaussian Smearing
        ROOT.gROOT.ProcessLine('TFile* gaussSmoother_Data_njets_pt_central = TFile::Open("results/GaussSmoother/gaussSmoother_Data_njets_pt_central.root")')
        ROOT.gROOT.ProcessLine('TList* cdfs_Data_u1_njets_pt_central = (TList*)gaussSmoother_Data_njets_pt_central->Get("cdfs_Data_u1_njets_pt_central")')
        ROOT.gROOT.ProcessLine('TList* cdfs_Data_u2_njets_pt_central = (TList*)gaussSmoother_Data_njets_pt_central->Get("cdfs_Data_u2_njets_pt_central")')
        ROOT.gROOT.ProcessLine('TFile* gaussSmoother_DY_njets_pt_central = TFile::Open("results/GaussSmoother/gaussSmoother_DY_njets_pt_central.root")')
        ROOT.gROOT.ProcessLine('TList* cdfs_DY_u1_njets_pt_central = (TList*)gaussSmoother_DY_njets_pt_central->Get("cdfs_DY_u1_njets_pt_central")')
        ROOT.gROOT.ProcessLine('TList* cdfs_DY_u2_njets_pt_central = (TList*)gaussSmoother_DY_njets_pt_central->Get("cdfs_DY_u2_njets_pt_central")')

        DYSamp.Define("u1_corr_central_GKS",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, cdfs_Data_u1_njets_pt_central, cdfs_DY_u1_njets_pt_central, 0.00001)")
        DYSamp.Define("u2_corr_central_GKS",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, cdfs_Data_u2_njets_pt_central, cdfs_DY_u2_njets_pt_central, 0.00001)")
        DYSamp.Define("u_pt_corr_central_GKS", "TMath::Sqrt(u1_corr_central_GKS*u1_corr_central_GKS + u2_corr_central_GKS*u2_corr_central_GKS)")
        sampMan.DefineAll("u1_corr_central_GKS",     "u1"  , excludes=['DY'])
        sampMan.DefineAll("u2_corr_central_GKS",     "u2"  , excludes=['DY'])
        sampMan.DefineAll("u_pt_corr_central_GKS",   "u_pt", excludes=['DY'])
        sampMan.DefineAll("deepmet_corr_central_GKS", "METVec(Z_pt, Z_phi, u1_corr_central_GKS, u2_corr_central_GKS)")
        sampMan.DefineAll("deepmet_pt_corr_central_GKS", "deepmet_corr_central_GKS.Mod()") 
        sampMan.DefineAll("deepmet_phi_corr_central_GKS", "TVector2::Phi_mpi_pi(deepmet_corr_central_GKS.Phi())")
        
    if doBkgScaled and doDefaultCorrection:
        ROOT.gROOT.ProcessLine('TFile* fitfunctions_Data_bkgScale = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_central_bkgScale.root")')
        ROOT.gROOT.ProcessLine('TH1F* h1_ptbins_Data_central_bkgScale   = (TH1F*)fitfunctions_Data_bkgScale->Get("h1_ptbins_Data_central_bkgScale")')
        ROOT.gROOT.ProcessLine('TH1F* h1_njetbins_Data_central_bkgScale = (TH1F*)fitfunctions_Data_bkgScale->Get("h1_njetbins_Data_central_bkgScale")')
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u1_njets_pt_central_bkgScale  = (TList*)fitfunctions_Data_bkgScale->Get("tfs_Data_u1_njets_pt_central_bkgScale")')
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u2_njets_pt_central_bkgScale  = (TList*)fitfunctions_Data_bkgScale->Get("tfs_Data_u2_njets_pt_central_bkgScale")')
        
        # assume the central correction is the same as the default correction
        # and is already running
        DYSamp.Define("u1_corr_central_bkgScale",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central_bkgScale, h1_ptbins_Data_central_bkgScale, tfs_Data_u1_njets_pt_central_bkgScale, tfs_MC_u1_njets_pt_central, 0.00001)")
        DYSamp.Define("u2_corr_central_bkgScale",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central_bkgScale, h1_ptbins_Data_central_bkgScale, tfs_Data_u2_njets_pt_central_bkgScale, tfs_MC_u2_njets_pt_central, 0.00001)")
        DYSamp.Define("u_pt_corr_central_bkgScale", "TMath::Sqrt(u1_corr_central_bkgScale*u1_corr_central_bkgScale + u2_corr_central_bkgScale*u2_corr_central_bkgScale)")
        sampMan.DefineAll("u1_corr_central_bkgScale",     "u1"  , excludes=['DY'])
        sampMan.DefineAll("u2_corr_central_bkgScale",     "u2"  , excludes=['DY'])
        sampMan.DefineAll("u_pt_corr_central_bkgScale",   "u_pt", excludes=['DY'])
        sampMan.DefineAll("deepmet_corr_central_bkgScale", "METVec(Z_pt, Z_phi, u1_corr_central_bkgScale, u2_corr_central_bkgScale)")
        sampMan.DefineAll("deepmet_pt_corr_central_bkgScale", "deepmet_corr_central_bkgScale.Mod()") 
        sampMan.DefineAll("deepmet_phi_corr_central_bkgScale", "TVector2::Phi_mpi_pi(deepmet_corr_central_bkgScale.Phi())")
        

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    # z kinematics
    #sampMan.cacheDraw("Z_pt", "histo_zjets_zpt_WoZptWeight", 30, 0, 60, DrawConfig(xmin=0, xmax=60, xlabel='p^{ll}_{T} [GeV]'), weightname="weight_WoVpt")
    sampMan.cacheDraw("Z_pt", "histo_zjets_zpt", 30, 0, 60, DrawConfig(xmin=0, xmax=60, xlabel='p^{ll}_{T} [GeV]'))
    sampMan.cacheDraw("Z_eta", "histo_zjets_zeta", 30, -3, 3, DrawConfig(xmin=-3, xmax=3, xlabel='#eta^{ll}'))
    sampMan.cacheDraw("Z_phi", "histo_zjets_zphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='#phi^{ll}'))
    sampMan.cacheDraw("m_ll", "histo_zjets_mll", 30, 80, 100, DrawConfig(xmin=80, xmax=100, xlabel='m^{ll} [GeV]'))
    
    # muon kinematics
    sampMan.cacheDraw("leadMuon_pt", "histo_zjets_leadmuonpt", 30, 20, 80, DrawConfig(xmin=20, xmax=80, xlabel='p^{l}_{T} [GeV]'))
    sampMan.cacheDraw("leadMuon_eta", "histo_zjets_leadmuoneta", 30, -3, 3, DrawConfig(xmin=-3, xmax=3, xlabel='#eta^{l}'))
    sampMan.cacheDraw("leadMuon_phi", "histo_zjets_leadmuonphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='#phi^{l}'))
    sampMan.cacheDraw("subleadMuon_pt", "histo_zjets_subleadmuonpt", 30, 20, 80, DrawConfig(xmin=20, xmax=80, xlabel='p^{l}_{T} [GeV]'))
    sampMan.cacheDraw("subleadMuon_eta", "histo_zjets_subleadmuoneta", 30, -3, 3, DrawConfig(xmin=-3, xmax=3, xlabel='#eta^{l}'))
    sampMan.cacheDraw("subleadMuon_phi", "histo_zjets_subleadmuonphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='#phi^{l}'))
    
    # met
    sampMan.cacheDraw("MET_pt", "histo_zjets_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=150, xlabel='PF MET [GeV]', yrmin=0.4, yrmax=1.6, addOverflow=True, addUnderflow=True))
    sampMan.cacheDraw("MET_phi", "histo_zjets_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF #phi', ymax=1e10))
    sampMan.cacheDraw("PuppiMET_pt", "histo_zjets_puppimet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=150, xlabel='PUPPI MET [GeV]', yrmin=0.4, yrmax=1.6, addOverflow=True, addUnderflow=True))
    sampMan.cacheDraw("PuppiMET_phi", "histo_zjets_puppimet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PUPPI #phi', ymax=1e10))
    sampMan.cacheDraw("DeepMETResolutionTune_pt", "histo_zjets_deepmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET [GeV]', yrmin=0.4, yrmax=1.6, addOverflow=True, addUnderflow=True))
    sampMan.cacheDraw("DeepMETResolutionTune_phi", "histo_zjets_deepmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET #phi', ymax=1e10, yrmin=0.4, yrmax=1.6))
    
    # recoil
    sampMan.cacheDraw("u1", "histo_zjets_u1", u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} [GeV]', addOverflow=True, addUnderflow=True))
    sampMan.cacheDraw("u2", "histo_zjets_u2", u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp } [GeV]', addOverflow=True, addUnderflow=True))
    sampMan.cacheDraw("u_pt",    "histo_zjets_u_pt", u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u [GeV]', addOverflow=True, addUnderflow=True))
    
    # corrected deepmet
    def DrawCorrection(postfix):
        sampMan.cacheDraw("deepmet_pt_corr_"+postfix, "histo_zjets_deepmet_pt_corr_"+postfix, met_pt_bins, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET Corrected [GeV]', yrmin=0.4, yrmax=1.6, addOverflow=True, addUnderflow=True))
        sampMan.cacheDraw("deepmet_phi_corr_"+postfix, "histo_zjets_deepmet_phi_corr_"+postfix, 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET Corrected #phi'))
        sampMan.cacheDraw("u1_corr_"+postfix, "histo_zjets_u1_corr_"+postfix, u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]', addOverflow=True, addUnderflow=True))
        sampMan.cacheDraw("u2_corr_"+postfix, "histo_zjets_u2_corr_"+postfix, u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp } Corrected [GeV]', addOverflow=True, addUnderflow=True))
        sampMan.cacheDraw("u_pt_corr_"+postfix,    "histo_zjets_u_pt_corr_"+postfix, u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u Corrected [GeV]', addOverflow=True, addUnderflow=True))

    if doDefaultCorrection:
        DrawCorrection("central")
        DrawCorrection("central_GKS")
        DrawCorrection("central_bkgScale")
    
    sampMan.launchDraw()
    
    # compare ratios
    rmin = 0.7
    rmax = 1.3

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
        DrawHistos(histos_u2, legends, -80.0, 80., "u_{#perp } Corrected [GeV]", 0.7, 1.3, "Data / MC", "histo_zjets_ratio_u2_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols) 

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
        DrawHistos(histos_u2, legends, -80.0, 80., "u_{#perp } Corrected [GeV]", 0.7, 1.3, "Ratio", "histo_zjets_MCratio_u2_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)

    colors = [1,2, 3]
    linestyles = [1,1, 1]

    legends = ["2-Gaussian Fit", "Gaussian Smoothing", "With Background Scaling"]
    compRatios(  ["corr_central", "corr_central_GKS", "corr_central_bkgScale"], colors, linestyles, legends, "FitVSSmoothing")
    compMCRatios(["corr_central", "corr_central_GKS", "corr_central_bkgScale"], colors, linestyles, legends, "FitVSSmoothing")

    print("Program end...")

    input()
    
    return 

if __name__ == "__main__":
   main()
