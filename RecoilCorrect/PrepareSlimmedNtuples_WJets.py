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
                      name="Data", bjetVeto=VetoB, isZSR=False, isWSR=True)
    sampMan = SampleManager(DataSamp, [])

    sampMan.DefineAll("leadMuon_pt", "Muon_pt[0]")
    sampMan.DefineAll("leadMuon_eta", "Muon_eta[0]")
    sampMan.DefineAll("leadMuon_phi", "Muon_phi[0]")


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

    sampMan.cacheDraw("MET_pt", "histo_wjets_pfmet_pt", met_pt_bins, DrawConfig(
        xmin=0, xmax=200, xlabel='PF MET [GeV]', showratio=False))
    sampMan.cacheDraw("MET_phi", "histo_wjets_pfmet_phi", 30, phimin, phimax, DrawConfig(
        xmin=phimin, xmax=phimax, xlabel='PF #phi', ymax=1e10, showratio=False))
    sampMan.cacheDraw("PuppiMET_pt", "histo_wjets_puppimet_pt", met_pt_bins, DrawConfig(
        xmin=0, xmax=200, xlabel='PUPPI MET [GeV]', showratio=False))
    sampMan.cacheDraw("PuppiMET_phi", "histo_wjets_puppimet_phi", 30, phimin, phimax, DrawConfig(
        xmin=phimin, xmax=phimax, xlabel='PUPPI #phi', ymax=1e10, showratio=False))

    # uncorrected deepmet
    sampMan.cacheDraw("DeepMETResolutionTune_pt", "histo_wjets_deepmet_pt",
                      met_pt_bins, DrawConfig(xmin=0, xmax=200, xlabel='Deep MET [GeV]', showratio=False))
    sampMan.cacheDraw("DeepMETResolutionTune_phi", "histo_wjets_deepmet_phi", 30, phimin,
                      phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Deep MET #phi', ymax=1e10, showratio=False))
    sampMan.cacheDraw("jet_n", "histo_wjets_jet_n", 10, 0, 10,
                      DrawConfig(xmin=0, xmax=10, xlabel='Jet Multiplicity', showratio=False))

    jets_variables_to_keep = None
    jets_variables_to_keep = sampMan.slimJets()

    sampMan.launchDraw()

    print("Save pre-selected and pre-processed variables to root file...")

    branches = [
                "leadMuon_pt", "leadMuon_eta", "leadMuon_phi",
                "DeepMETResolutionTune_pt", "DeepMETResolutionTune_phi",
                "MET_pt", "MET_phi", "MET_significance", "MET_sumEt",
                "RawMET_pt", "RawMET_phi", "RawMET_sumEt",
                "RawPuppiMET_pt", "RawPuppiMET_phi", "RawPuppiMET_sumEt",
                "PuppiMET_pt", "PuppiMET_phi", "PuppiMET_phiJERDown", "PuppiMET_phiJERUp", "PuppiMET_phiJESDown", "PuppiMET_phiJESUp", "PuppiMET_phiUnclusteredDown", "PuppiMET_phiUnclusteredUp", "PuppiMET_ptJERDown", "PuppiMET_ptJERUp", "PuppiMET_ptJESDown", "PuppiMET_ptJESUp", "PuppiMET_ptUnclusteredDown", "PuppiMET_ptUnclusteredUp", "PuppiMET_sumEt",
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
    outdir = "/afs/cern.ch/work/y/yofeng/public/outputroot"
    outdir = "/eos/cms/store/user/yofeng/SlimmedOutput_WJets/"
    sampMan.snapShot(outdir, branches, jets_variables=jets_variables_to_keep)

    print("Program end...")

    input()

    return


if __name__ == "__main__":
    main()
