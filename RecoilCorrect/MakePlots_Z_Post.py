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

def main():
    print("Program start...")

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

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

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150, 165, 180, 200])
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
    sampMan.cacheDraw("MET_pt", "histo_zjets_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='PF MET [GeV]'))
    sampMan.cacheDraw("MET_phi", "histo_zjets_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF #phi', ymax=1e10))
    sampMan.cacheDraw("PuppiMET_pt", "histo_zjets_puppimet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='PUPPI MET [GeV]'))
    sampMan.cacheDraw("PuppiMET_phi", "histo_zjets_puppimet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PUPPI #phi', ymax=1e10))
    sampMan.cacheDraw("DeepMETResolutionTune_pt", "histo_zjets_deepmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET [GeV]'))
    sampMan.cacheDraw("DeepMETResolutionTune_phi", "histo_zjets_deepmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET #phi', ymax=1e10))
    sampMan.cacheDraw("deepmet_pt_corr_central", "histo_zjets_deepmet_pt_corr_central", met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET Corrected [GeV]'))
    sampMan.cacheDraw("deepmet_phi_corr_central", "histo_zjets_deepmet_phi_corr_central", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET Corrected #phi', ymax=1e10))
    
    # recoil
    sampMan.cacheDraw("u1", "histo_zjets_u1", u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} [GeV]'))
    sampMan.cacheDraw("u2", "histo_zjets_u2", u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} [GeV]'))
    sampMan.cacheDraw("u_pt",    "histo_zjets_u_pt", u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u [GeV]'))
    sampMan.cacheDraw("u1_corr_central", "histo_zjets_u1_corr_central", u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]'))
    sampMan.cacheDraw("u2_corr_central", "histo_zjets_u2_corr_central", u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} Corrected [GeV]'))
    sampMan.cacheDraw("u_pt_corr_central",    "histo_zjets_u_pt_corr_central", u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u Corrected [GeV]'))

    # corrected deepmet
    def DrawCorrection(postfix):
        sampMan.cacheDraw("deepmet_pt_corr_"+postfix, "histo_zjets_deepmet_pt_corr_"+postfix, met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET Corrected [GeV]'))
        sampMan.cacheDraw("deepmet_phi_corr_"+postfix, "histo_zjets_deepmet_phi_corr_"+postfix, 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET Corrected #phi'))
        sampMan.cacheDraw("u1_corr_"+postfix, "histo_zjets_u1_corr_"+postfix, u1_bins, DrawConfig(xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]'))
        sampMan.cacheDraw("u2_corr_"+postfix, "histo_zjets_u2_corr_"+postfix, u2_bins, DrawConfig(xmin=-80., xmax=80., xlabel='u_{#perp} Corrected [GeV]'))
        sampMan.cacheDraw("u_pt_corr_"+postfix,    "histo_zjets_u_pt_corr_"+postfix, u_bins,  DrawConfig(xmin=0, xmax=150, xlabel='u Corrected [GeV]'))

    #DrawCorrection("central")
    
    sampMan.launchDraw()

    print("Program end...")

    input()
    
    return 

if __name__ == "__main__":
   main()
