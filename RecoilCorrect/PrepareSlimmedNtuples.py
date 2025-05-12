'''
Apply the event selections and save the slimmed ntuples with 
critical variables for downstream calibrations, uncertainties, and plotting
'''
from SampleManager import DrawConfig, Sample, SampleManager
import pickle
import CMS_lumi
from myFunction import DrawHistos, THStack2TH1
from tdrstyle import setTDRStyle
import ROOT
import numpy as np
from collections import OrderedDict
import sys
sys.path.append("../RecoilResol/CMSPLOTS")


ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)

dotest = False
VetoB = False
doDefaultCorrection = False


def main():
    print("Program start...")

    if dotest:
        input_data = "inputs/inputs_Z_UL/input_data_test.txt"
        input_dy = "inputs/inputs_Z_UL/input_zjets.txt"
        input_ttbar = "inputs/inputs_Z_UL/input_ttbar_test.txt"
    else:
        input_data = "inputs/inputs_Z_UL/input_data.txt"
        input_dy = "inputs/inputs_Z_UL/input_zjets_all.txt"
        input_ttbar = "inputs/inputs_Z_UL/input_ttbar.txt"
        input_WW2L = "inputs/inputs_Z_UL/input_WWTo2L2Nu.txt"
        input_WZ2L = "inputs/inputs_Z_UL/input_WZTo2Q2L.txt"
        input_ZZ2L = "inputs/inputs_Z_UL/input_ZZTo2L2Nu.txt"
        input_ZZ2L2Q = "inputs/inputs_Z_UL/input_ZZTo2Q2L.txt"
        input_dytau = "inputs/inputs_Z_UL/input_zjets_tautau.txt"

    DataSamp = Sample(input_data, isMC=False, legend="Data",
                      name="Data", bjetVeto=VetoB)
    DYSamp = Sample(input_dy,    xsec=2025.74*1e3,      color=5,
                    reweightzpt=False, legend="DY", name="DY", bjetVeto=VetoB)
    TTbarSamp = Sample(input_ttbar, xsec=831.76*0.105*1e3, color=46,
                       reweightzpt=False, legend="t#bar{t}", name="ttbar", bjetVeto=VetoB)
    if not dotest:
        WW2LSamp = Sample(input_WW2L,  xsec=12.178*1e3,       color=38,
                          reweightzpt=False, legend="WW2L",     name="WW2L",  bjetVeto=VetoB)
        WZ2LSamp = Sample(input_WZ2L,  xsec=5.26*1e3,         color=39,
                          reweightzpt=False, legend="WZ2L",     name="WZ2L",  bjetVeto=VetoB)
        ZZ2LSamp = Sample(input_ZZ2L,  xsec=0.564*1e3,        color=37,
                          reweightzpt=False, legend="ZZ2L",     name="ZZ2L",  bjetVeto=VetoB)
        ZZ2L2QSamp = Sample(input_ZZ2L2Q,  xsec=1.212*1e3,        color=36,
                            reweightzpt=False, legend="ZZ2L2Q",     name="ZZ2L2Q",  bjetVeto=VetoB)
        DYTauSamp = Sample(input_dytau, xsec=2025.74*1e3,      color=8,  reweightzpt=False,
                           legend="DY#rightarrow#tau#tau", name="DYTauTau", bjetVeto=VetoB)

    if not dotest:
        sampMan = SampleManager(
            DataSamp, [DYSamp, DYTauSamp, WW2LSamp, WZ2LSamp, TTbarSamp, ZZ2LSamp, ZZ2L2QSamp])
    else:
        sampMan = SampleManager(DataSamp, [DYSamp, TTbarSamp])
    sampMan.groupMCs(["WW2L", "WZ2L", "ZZ2L", "ZZ2L2Q"],
                     "Dibosons", 38, "Dibosons")

    sampMan.DefineAll("zpt", "Z_pt")
    sampMan.DefineAll("leadMuon_pt", "Muon_pt[0]")
    sampMan.DefineAll("leadMuon_eta", "Muon_eta[0]")
    sampMan.DefineAll("leadMuon_phi", "Muon_phi[0]")
    sampMan.DefineAll("subleadMuon_pt", "Muon_pt[1]")
    sampMan.DefineAll("subleadMuon_eta", "Muon_eta[1]")
    sampMan.DefineAll("subleadMuon_phi", "Muon_phi[1]")

    if doDefaultCorrection:
        # Default correction:
        ROOT.gROOT.ProcessLine(
            'TFile* fitfunctions_Data = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_central.root")')
        ROOT.gROOT.ProcessLine(
            'TH1F* h1_ptbins_Data_central   = (TH1F*)fitfunctions_Data->Get("h1_ptbins_Data_central")')
        ROOT.gROOT.ProcessLine(
            'TH1F* h1_njetbins_Data_central = (TH1F*)fitfunctions_Data->Get("h1_njetbins_Data_central")')
        ROOT.gROOT.ProcessLine(
            'TList* tfs_Data_u1_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u1_njets_pt_central")')
        ROOT.gROOT.ProcessLine(
            'TList* tfs_Data_u2_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u2_njets_pt_central")')
        ROOT.gROOT.ProcessLine(
            'TFile* fitfunctions_DY_central = TFile::Open("results/Fit/fitfunctions_DY_njets_pt_central.root")')
        ROOT.gROOT.ProcessLine(
            'TList* tfs_MC_u1_njets_pt_central = (TList*)fitfunctions_DY_central->Get("tfs_DY_u1_njets_pt_central")')
        ROOT.gROOT.ProcessLine(
            'TList* tfs_MC_u2_njets_pt_central = (TList*)fitfunctions_DY_central->Get("tfs_DY_u2_njets_pt_central")')

        DYSamp.Define("u1_corr_central",   "UCorrection_Quant(u1, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_central, tfs_MC_u1_njets_pt_central, 0.00001)")
        DYSamp.Define("u2_corr_central",   "UCorrection_Quant(u2, jet_n, Z_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_central, tfs_MC_u2_njets_pt_central, 0.00001)")
        DYSamp.Define("u_pt_corr_central",
                      "TMath::Sqrt(u1_corr_central*u1_corr_central + u2_corr_central*u2_corr_central)")
        sampMan.DefineAll("u1_corr_central",     "u1", excludes=['DY'])
        sampMan.DefineAll("u2_corr_central",     "u2", excludes=['DY'])
        sampMan.DefineAll("u_pt_corr_central",   "u_pt", excludes=['DY'])
        sampMan.DefineAll("deepmet_corr_central",
                          "METVec(Z_pt, Z_phi, u1_corr_central, u2_corr_central)")
        sampMan.DefineAll("deepmet_pt_corr_central",
                          "deepmet_corr_central.Mod()")
        sampMan.DefineAll("deepmet_phi_corr_central",
                          "TVector2::Phi_mpi_pi(deepmet_corr_central.Phi())")

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33,
                           36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150, 165, 180, 200])
    u1_bins = np.array([-40., -36., -32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10,
                       12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -
                       12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40,
                      43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    # z pt befor and after pt reweighting
    sampMan.cacheDraw("zpt", "histo_zjets_zpt_WoZptWeight", 30, 0, 60, DrawConfig(
        xmin=0, xmax=60, xlabel='p^{ll}_{T} [GeV]'), weightname="weight_WoVpt")
    sampMan.cacheDraw("zpt", "histo_zjets_zpt", 30, 0, 60, DrawConfig(
        xmin=0, xmax=60, xlabel='p^{ll}_{T} [GeV]'))

    sampMan.cacheDraw("MET_pt", "histo_zjets_pfmet_pt", met_pt_bins, DrawConfig(
        xmin=0, xmax=200, xlabel='PF MET [GeV]'))
    sampMan.cacheDraw("MET_phi", "histo_zjets_pfmet_phi", 30, phimin, phimax, DrawConfig(
        xmin=phimin, xmax=phimax, xlabel='PF #phi', ymax=1e10))
    sampMan.cacheDraw("PuppiMET_pt", "histo_zjets_puppimet_pt", met_pt_bins, DrawConfig(
        xmin=0, xmax=200, xlabel='PUPPI MET [GeV]'))
    sampMan.cacheDraw("PuppiMET_phi", "histo_zjets_puppimet_phi", 30, phimin, phimax, DrawConfig(
        xmin=phimin, xmax=phimax, xlabel='PUPPI #phi', ymax=1e10))

    # uncorrected deepmet
    sampMan.cacheDraw("DeepMETResolutionTune_pt", "histo_zjets_deepmet_pt",
                      met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET [GeV]'))
    sampMan.cacheDraw("DeepMETResolutionTune_phi", "histo_zjets_deepmet_phi", 30, phimin,
                      phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET #phi', ymax=1e10))
    sampMan.cacheDraw("u1", "histo_zjets_u1", u1_bins, DrawConfig(
        xmin=-40.0, xmax=100, xlabel='u_{#parallel} [GeV]'))
    sampMan.cacheDraw("u2", "histo_zjets_u2", u2_bins, DrawConfig(
        xmin=-80., xmax=80., xlabel='u_{#perp} [GeV]'))
    sampMan.cacheDraw("u_pt",    "histo_zjets_u_pt", u_bins,
                      DrawConfig(xmin=0, xmax=150, xlabel='u [GeV]'))
    sampMan.cacheDraw("jet_n", "histo_zjets_jet_n", 10, 0, 10,
                      DrawConfig(xmin=0, xmax=10, xlabel='Jet Multiplicity'))

    # corrected deepmet
    def DrawCorrection(postfix):
        sampMan.cacheDraw("deepmet_pt_corr_"+postfix, "histo_zjets_deepmet_pt_corr_"+postfix, met_pt_bins, DrawConfig(
            xmin=0, xmax=150, xlabel='Deep MET Corrected [GeV]', yrmin=0.4, yrmax=1.6, addOverflow=True, addUnderflow=True))
        sampMan.cacheDraw("deepmet_phi_corr_"+postfix, "histo_zjets_deepmet_phi_corr_"+postfix, 30,
                          phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET Corrected #phi'))
        sampMan.cacheDraw("u1_corr_"+postfix, "histo_zjets_u1_corr_"+postfix, u1_bins, DrawConfig(
            xmin=-40.0, xmax=100, xlabel='u_{#parallel} Corrected [GeV]', addOverflow=True, addUnderflow=True))
        sampMan.cacheDraw("u2_corr_"+postfix, "histo_zjets_u2_corr_"+postfix, u2_bins, DrawConfig(
            xmin=-80., xmax=80., xlabel='u_{#perp } Corrected [GeV]', addOverflow=True, addUnderflow=True))
        sampMan.cacheDraw("u_pt_corr_"+postfix,    "histo_zjets_u_pt_corr_"+postfix, u_bins,  DrawConfig(
            xmin=0, xmax=150, xlabel='u Corrected [GeV]', addOverflow=True, addUnderflow=True))

    if doDefaultCorrection:
        DrawCorrection("central")
        
    jets_variables_to_keep = None
    jets_variables_to_keep = sampMan.slimJets()

    sampMan.launchDraw()

    print("Save pre-selected and pre-processed variables to root file...")

    branches = ["Z_pt", "Z_eta", "Z_phi", "m_ll",
                "leadMuon_pt", "leadMuon_eta", "leadMuon_phi",
                "subleadMuon_pt", "subleadMuon_eta", "subleadMuon_phi",
                "DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi",
                "MET_pt", "MET_phi", "MET_significance", "MET_sumEt",
                "RawMET_pt", "RawMET_phi", "RawMET_sumEt",
                "RawPuppiMET_pt", "RawPuppiMET_phi", "RawPuppiMET_sumEt",
                "PuppiMET_pt", "PuppiMET_phi", "PuppiMET_phiJERDown", "PuppiMET_phiJERUp", "PuppiMET_phiJESDown", "PuppiMET_phiJESUp", "PuppiMET_phiUnclusteredDown", "PuppiMET_phiUnclusteredUp", "PuppiMET_ptJERDown", "PuppiMET_ptJERUp", "PuppiMET_ptJESDown", "PuppiMET_ptJESUp", "PuppiMET_ptUnclusteredDown", "PuppiMET_ptUnclusteredUp", "PuppiMET_sumEt",
                "u1", "u2", "u_pt", "u_phi",
                "weight", "weight_WoVpt", "PV_npvs", "PV_npvsGood",
                "jet_n", "jet_CSVLoose_n", "jet_CSVMedium_n", "jet_CSVTight_n",
                "metFilters"
                ]
    # more pileup related variables
    branches += [
        "fixedGridRhoFastjetAll", "fixedGridRhoFastjetCentral", "fixedGridRhoFastjetCentralNeutral", "fixedGridRhoFastjetCentralChargedPileUp", "fixedGridRhoFastjetCentralCalo"
    ]
    # more deepmet related variables
    branches += ["DeepMETPVRobust_pt", "DeepMETPVRobust_phi",
                 "DeepMETPVRobustNoPUPPI_pt", "DeepMETPVRobustNoPUPPI_phi",
                 "DeepMETResponseTune_pt", "DeepMETResponseTune_phi"]
    if doDefaultCorrection:
        branches += ["u1_corr_central", "u2_corr_central", "u_pt_corr_central",
                     "deepmet_pt_corr_central", "deepmet_phi_corr_central"]
    outdir = "/afs/cern.ch/work/y/yofeng/public/outputroot"
    outdir = "/eos/cms/store/user/yofeng/SlimmedOutput/"
    sampMan.snapShot(outdir, branches, jets_variables=jets_variables_to_keep)

    print("Program end...")

    input()

    return


if __name__ == "__main__":
    main()
