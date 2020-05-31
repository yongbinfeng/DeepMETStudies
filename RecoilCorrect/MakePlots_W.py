import ROOT
import numpy as np
from collections import OrderedDict
import sys
sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
from tdrstyle import setTDRStyle
from myFunction import DrawHistos, THStack2TH1
import CMS_lumi
import pickle

from SampleManager import DrawConfig, Sample, SampleManager

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

dotest = 0

def main():
    print "Program start..."

    if dotest:
        input_data    = "inputs/inputs_W/input_data_test.txt"
        input_wjets   = "inputs/inputs_W/input_wjets_test.txt"
        # for the data-driven bkg estimation from the QCD dominated region
        input_data_qcdcr    = "inputs/inputs_QCDCR/input_data_test.txt"
        input_wjets_qcdcr   = "inputs/inputs_QCDCR/input_wjets_test.txt"
        input_zjets_qcdcr   = "inputs/inputs_QCDCR/input_zjets_test.txt"
        input_ttbar_qcdcr   = "inputs/inputs_QCDCR/input_ttbar_test.txt"
        input_lept_qcdcr    = "inputs/inputs_QCDCR/input_ttbar_singlelep_fromT_test.txt"
        input_leptbar_qcdcr = "inputs/inputs_QCDCR/input_ttbar_singlelep_fromTbar_test.txt"
    else:
        input_data    = "inputs/inputs_W/input_data.txt"
        input_wjets   = "inputs/inputs_W/input_wjets.txt"
        input_zjets   = "inputs/inputs_W/input_zjets.txt"
        input_ttbar   = "inputs/inputs_W/input_ttbar.txt"
        input_lept    = "inputs/inputs_W/input_ttbar_singlelep_fromT.txt"
        input_leptbar = "inputs/inputs_W/input_ttbar_singlelep_fromTbar.txt"
        #input_qcd     = "inputs/inputs_W/input_qcd.txt"
        # for the data-driven bkg estimation from the QCD dominated region
        input_data_qcdcr    = "inputs/inputs_QCDCR/input_data.txt"
        input_wjets_qcdcr   = "inputs/inputs_QCDCR/input_wjets.txt"
        input_zjets_qcdcr   = "inputs/inputs_QCDCR/input_zjets.txt"
        input_ttbar_qcdcr   = "inputs/inputs_QCDCR/input_ttbar.txt"
        input_lept_qcdcr    = "inputs/inputs_QCDCR/input_ttbar_singlelep_fromT.txt"
        input_leptbar_qcdcr = "inputs/inputs_QCDCR/input_ttbar_singlelep_fromTbar.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", isZSR=False, isWSR=True)
    WJetsSamp = Sample(input_wjets,    xsec = 61526.7*1e3,      color=38,  legend="WJets", name="WJets", isZSR=False, isWSR=True)
    if not dotest:
        ZJetsSamp = Sample(input_zjets, xsec = 6077.22*1e3,      color=5,  legend="DYJets", name="ZJets", isZSR=False, isWSR=True)
        TTbarSamp = Sample(input_ttbar, xsec = 831.76*0.105*1e3, color=46, legend="t#bar{t}", name="ttbar2lep", isZSR=False, isWSR=True)
        TTLepTSamp = Sample(input_lept, xsec = 831.76*0.219*1e3, color=47, legend="t#bar{t}", name="ttbarlept", isZSR=False, isWSR=True)
        TTLepTbarSamp = Sample(input_leptbar, xsec = 831.76*0.219*1e3, color=48, legend="t#bar{t}", name="ttbarleptbar", isZSR=False, isWSR=True)
        #QCDSamp   = Sample(input_qcd,   xsec = 720648000.0*0.00042*1e3, color=8, legend="QCD", name="QCD", isZSR=False, isWSR=True)

    DataQCDCRSamp  = Sample(input_data_qcdcr,  isMC=False, legend="Data_QCDCR", name="Data_QCDCR", isZSR=False, isWSR=True)
    WJetsQCDCRSamp = Sample(input_wjets_qcdcr, xsec = -61526.7*1e3, color=38, legend="WJets_QCDCR", name="WJets_QCDCR", isZSR=False, isWSR=True, applySF=False)
    ZJetsQCDCRSamp = Sample(input_zjets_qcdcr, xsec = -6077.22*1e3, color=5,  legend="ZJets_QCDCR", name="ZJets_QCDCR", isZSR=False, isWSR=True, applySF=False)
    TTbarQCDCRSamp = Sample(input_ttbar_qcdcr, xsec = -831.76*0.105*1e3, color=46, legend="t#bar{t} QCDCR", name="ttbar2lep_QCDCR", isZSR=False, isWSR=True, applySF=False)
    TTLepTQCDCRSamp = Sample(input_lept_qcdcr, xsec = -831.76*0.219*1e3, color=47, legend="t#bar{t} QCDCR", name="ttbarlept_QCDCR", isZSR=False, isWSR=True, applySF=False)
    TTLepTbarQCDCRSamp = Sample(input_leptbar_qcdcr, xsec = -831.76*0.219*1e3, color=48, legend="t#bar{t} QCDCR", name="ttbarleptb", isZSR=False, isWSR=True, applySF=False)

    if not dotest:
        sampMan = SampleManager(DataSamp, [WJetsSamp, ZJetsSamp, TTbarSamp, TTLepTSamp, TTLepTbarSamp, DataQCDCRSamp, WJetsQCDCRSamp, ZJetsQCDCRSamp, TTbarQCDCRSamp, TTLepTQCDCRSamp, TTLepTbarQCDCRSamp])
    else:
        sampMan = SampleManager(DataSamp, [WJetsSamp, DataQCDCRSamp, WJetsQCDCRSamp, ZJetsQCDCRSamp, TTbarQCDCRSamp, TTLepTQCDCRSamp, TTLepTbarQCDCRSamp])

    # create renormalized QCD scale variation weights for WJets
    def prepareWJetsQCDScaleVariations():
        print "prepare the renormalized QCD scale variations"
        QCDVariations = OrderedDict({"central": 0, "mur_1_muf_2":1, "mur_1_muf_0p5": 2, "mur_2_muf_1": 3, "mur_2_muf_2": 4, "mur_0p5_muf_1": 6, "mur_0p5_muf_0p5": 8})
        hcounts = OrderedDict()
        WJetsSamp.Define("const_1", "1.0")
        for var, index in QCDVariations.iteritems():
            WJetsSamp.Define("weight_"+var, "weight * EventWeights[%d]/EventWeights[0]"%index)
            hcounts[var] = WJetsSamp.rdf.Histo1D(("hcounts_"+var, "hcounts_"+var, 2, 0, 2), "const_1", "weight_"+var)

        nevents_norm = hcounts['central'].Integral()
        for var, hist in hcounts.iteritems():
            nevents = hcounts[var].Integral()
            print "For {}: Total: {}, weighted:{}. Norm factor {}".format(var, nevents_norm, nevents, nevents_norm/nevents)
            weightname = "weight_"+var
            WJetsSamp.Define(weightname+"_scaled", "{} * {}".format(weightname, nevents_norm/nevents))
        print "finished preparing renormalized QCD scale variations"
    prepareWJetsQCDScaleVariations()
            
    sampMan.groupMCs(["Data_QCDCR", "WJets_QCDCR", "ZJets_QCDCR", "ttbar2lep_QCDCR", "ttbarlept_QCDCR", "ttbarleptb"], "QCD", 17, "QCD", renormalize=False, renormalizefactor=3.29)
    sampMan.groupMCs(["ttbar2lep", "ttbarlept", "ttbarleptbar"], "ttbar", 46, "t#bar{t}")

    # QCD enriched isolation regions
    sampMan.ApplyCutSpecificMCs("(mu_pfIso[0]>0.30 && mu_pfIso[0]<0.60)", sampgroupnames=["QCD"])
    sampMan.DefineSpecificMCs("weight_Iso0", "weight * (mu_pfIso[0]>0.30 && mu_pfIso[0]<0.60)", sampgroupnames=["QCD"])
    sampMan.DefineAll("weight_Iso0", "weight", excludeGroups=['QCD'])
    sampMan.DefineSpecificMCs("weight_Iso1", "weight * (mu_pfIso[0]>0.30 && mu_pfIso[0]<0.60) * 3.32/3.29", sampgroupnames=["QCD"])
    sampMan.DefineAll("weight_Iso1", "weight", excludeGroups=['QCD'])
    sampMan.DefineSpecificMCs("weight_Iso2", "weight * (mu_pfIso[0]>0.30 && mu_pfIso[0]<0.60) * 3.20/3.29", sampgroupnames=["QCD"])
    sampMan.DefineAll("weight_Iso2", "weight", excludeGroups=['QCD'])

    sampMan.DefineAll("Mu_pt", "mu_pt[0]")
    sampMan.DefineAll("mT_PF",    "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], met_pt, met_phi)")
    sampMan.DefineAll("mT_Puppi", "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], puppimet_pt, puppimet_phi)")
    sampMan.DefineAll("mT_Deep",  "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt,  deepmet_phi)")
    sampMan.DefineAll("u_reco_pt_PF",    "calUT_PtPhi(mu_pt[0], mu_phi[0], met_pt, met_phi)")
    sampMan.DefineAll("u_reco_pt_Puppi", "calUT_PtPhi(mu_pt[0], mu_phi[0], puppimet_pt, puppimet_phi)")
    sampMan.DefineAll("u_reco_pt_Deep",  "calUT_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt,  deepmet_phi)")

    # correct DeepMET
    WJetsSamp.Define("trueW_index", "trueWIndex(trueW_n)")
    WJetsSamp.Define("W_pt", "(trueW_index >=0 && trueW_index<=1) ? trueW_pt[trueW_index] : 0")
    WJetsSamp.Define("W_phi", "(trueW_index>=0 && trueW_index<=1) ? trueW_phi[trueW_index] : 0")
    WJetsSamp.Define("Uvec", "UVec(mu_pt[0], mu_phi[0], deepmet_pt, deepmet_phi)")
    WJetsSamp.Define("u_pt",  "Uvec.Mod()")     
    WJetsSamp.Define("u_phi", "Uvec.Phi()")
    WJetsSamp.Define("u1",    "u_pt * TMath::Cos(u_phi + TMath::Pi() - W_phi)")
    WJetsSamp.Define("u2",    "u_pt * TMath::Sin(u_phi + TMath::Pi() - W_phi)")
 
    ROOT.gROOT.ProcessLine('TFile* fitfunctions_Data = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TH1F* h1_ptbins_Data_central   = (TH1F*)fitfunctions_Data->Get("h1_ptbins_Data_central")')
    ROOT.gROOT.ProcessLine('TH1F* h1_njetbins_Data_central = (TH1F*)fitfunctions_Data->Get("h1_njetbins_Data_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_Data_u1_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_Data_u2_njets_pt_central  = (TList*)fitfunctions_Data->Get("tfs_Data_u2_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TFile* fitfunctions_ZJets_NLO_central = TFile::Open("results/Fit/fitfunctions_ZJets_NLO_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TList* tfs_MC_u1_njets_pt_central = (TList*)fitfunctions_ZJets_NLO_central->Get("tfs_ZJets_NLO_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_MC_u2_njets_pt_central = (TList*)fitfunctions_ZJets_NLO_central->Get("tfs_ZJets_NLO_u2_njets_pt_central")') 

    WJetsSamp.Define("u1_corr_central",   "UCorrection_Quant(u1, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_central, tfs_MC_u1_njets_pt_central, 0.00001)")
    WJetsSamp.Define("u2_corr_central",   "UCorrection_Quant(u2, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_central, tfs_MC_u2_njets_pt_central, 0.00001)")
    WJetsSamp.Define("deepmet_corr_central", "METVec(mu_pt[0], mu_phi[0], u1_corr_central, u2_corr_central, W_phi)")
    WJetsSamp.Define("deepmet_pt_corr_central", "deepmet_corr_central.Mod()")
    WJetsSamp.Define("deepmet_phi_corr_central", "TVector2::Phi_mpi_pi(deepmet_corr_central.Phi())")
    WJetsSamp.Define("mT_Deep_corr_central", "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_central, deepmet_phi_corr_central)")
    sampMan.DefineAll("deepmet_pt_corr_central", "deepmet_pt", excludes=['WJets'])
    sampMan.DefineAll("deepmet_phi_corr_central", "deepmet_phi", excludes=['WJets'])
    sampMan.DefineAll("mT_Deep_corr_central", "mT_Deep", excludes=['WJets'])
    sampMan.DefineAll("u_reco_pt_Deep_corr_central",  "calUT_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_central,  deepmet_phi_corr_central)")

    # systematics variations
    # corrections from MadGraph
    ROOT.gROOT.ProcessLine('TFile* fitfunctions_ZJets_MG_central = TFile::Open("results/Fit/fitfunctions_ZJets_MG_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TList* tfs_MC_MG_u1_njets_pt_central = (TList*)fitfunctions_ZJets_MG_central->Get("tfs_ZJets_MG_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* tfs_MC_MG_u2_njets_pt_central = (TList*)fitfunctions_ZJets_MG_central->Get("tfs_ZJets_MG_u2_njets_pt_central")')
    
    WJetsSamp.Define("u1_corr_central_MG", "UCorrection_Quant(u1, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_central, tfs_MC_MG_u1_njets_pt_central, 0.00001)")
    WJetsSamp.Define("u2_corr_central_MG", "UCorrection_Quant(u2, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_central, tfs_MC_MG_u2_njets_pt_central, 0.00001)")
    WJetsSamp.Define("deepmet_corr_central_MG", "METVec(mu_pt[0], mu_phi[0], u1_corr_central_MG, u2_corr_central_MG, W_phi)")
    WJetsSamp.Define("deepmet_pt_corr_central_MG", "deepmet_corr_central_MG.Mod()")
    WJetsSamp.Define("deepmet_phi_corr_central_MG", "TVector2::Phi_mpi_pi(deepmet_corr_central_MG.Phi())")
    WJetsSamp.Define("mT_Deep_corr_central_MG", "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_central_MG, deepmet_phi_corr_central_MG)")
    sampMan.DefineAll("deepmet_pt_corr_central_MG", "deepmet_pt", excludes=['WJets'])
    sampMan.DefineAll("deepmet_phi_corr_central_MG", "deepmet_phi", excludes=['WJets'])
    sampMan.DefineAll("mT_Deep_corr_central_MG", "mT_Deep", excludes=['WJets'])
    sampMan.DefineAll("u_reco_pt_Deep_corr_central_MG",  "calUT_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_central_MG,  deepmet_phi_corr_central_MG)")

    # corrections from Gaussian Smearing
    ROOT.gROOT.ProcessLine('TFile* gaussSmoother_Data_njets_pt_central = TFile::Open("results/GaussSmoother/gaussSmoother_Data_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TList* cdfs_Data_u1_njets_pt_central = (TList*)gaussSmoother_Data_njets_pt_central->Get("cdfs_Data_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* cdfs_Data_u2_njets_pt_central = (TList*)gaussSmoother_Data_njets_pt_central->Get("cdfs_Data_u2_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TFile* gaussSmoother_ZJets_NLO_njets_pt_central = TFile::Open("results/GaussSmoother/gaussSmoother_ZJets_NLO_njets_pt_central.root")')
    ROOT.gROOT.ProcessLine('TList* cdfs_ZJets_NLO_u1_njets_pt_central = (TList*)gaussSmoother_ZJets_NLO_njets_pt_central->Get("cdfs_ZJets_NLO_u1_njets_pt_central")')
    ROOT.gROOT.ProcessLine('TList* cdfs_ZJets_NLO_u2_njets_pt_central = (TList*)gaussSmoother_ZJets_NLO_njets_pt_central->Get("cdfs_ZJets_NLO_u2_njets_pt_central")')
    
    WJetsSamp.Define("u1_corr_central_GKS",   "UCorrection_Quant(u1, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, cdfs_Data_u1_njets_pt_central, cdfs_ZJets_NLO_u1_njets_pt_central, 0.00001)")
    WJetsSamp.Define("u2_corr_central_GKS",   "UCorrection_Quant(u2, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, cdfs_Data_u2_njets_pt_central, cdfs_ZJets_NLO_u2_njets_pt_central, 0.00001)")
    WJetsSamp.Define("deepmet_corr_central_GKS", "METVec(mu_pt[0], mu_phi[0], u1_corr_central_GKS, u2_corr_central_GKS, W_phi)")
    WJetsSamp.Define("deepmet_pt_corr_central_GKS", "deepmet_corr_central_GKS.Mod()")
    WJetsSamp.Define("deepmet_phi_corr_central_GKS", "TVector2::Phi_mpi_pi(deepmet_corr_central_GKS.Phi())")
    WJetsSamp.Define("mT_Deep_corr_central_GKS", "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_central_GKS, deepmet_phi_corr_central_GKS)")
    sampMan.DefineAll("deepmet_pt_corr_central_GKS", "deepmet_pt", excludes=['WJets'])
    sampMan.DefineAll("deepmet_phi_corr_central_GKS", "deepmet_phi", excludes=['WJets'])
    sampMan.DefineAll("mT_Deep_corr_central_GKS", "mT_Deep", excludes=['WJets'])
    sampMan.DefineAll("u_reco_pt_Deep_corr_central_GKS",  "calUT_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_central_GKS,  deepmet_phi_corr_central_GKS)")

    # correction from ttbar qcd scale:
    def varyTTbarQCDScale(wstring):
        ROOT.gROOT.ProcessLine('TFile* fitfunctions_Data_{WS} = TFile::Open("results/Fit/fitfunctions_Data_njets_pt_{WS}.root")'.format(WS=wstring))
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u1_njets_pt_{WS}  = (TList*)fitfunctions_Data_{WS}->Get("tfs_Data_u1_njets_pt_{WS}")'.format(WS=wstring))
        ROOT.gROOT.ProcessLine('TList* tfs_Data_u2_njets_pt_{WS}  = (TList*)fitfunctions_Data_{WS}->Get("tfs_Data_u2_njets_pt_{WS}")'.format(WS=wstring))

        WJetsSamp.Define("u1_corr_"+wstring,   "UCorrection_Quant(u1, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u1_njets_pt_{WS}, tfs_MC_u1_njets_pt_central, 0.00001)".format(WS=wstring))
        WJetsSamp.Define("u2_corr_"+wstring,   "UCorrection_Quant(u2, jet_n, W_pt, h1_njetbins_Data_central, h1_ptbins_Data_central, tfs_Data_u2_njets_pt_{WS}, tfs_MC_u2_njets_pt_central, 0.00001)".format(WS=wstring))
        WJetsSamp.Define("deepmet_corr_"+wstring, "METVec(mu_pt[0], mu_phi[0], u1_corr_{WS}, u2_corr_{WS}, W_phi)".format(WS=wstring))
        WJetsSamp.Define("deepmet_pt_corr_"+wstring, "deepmet_corr_{WS}.Mod()".format(WS=wstring))
        WJetsSamp.Define("deepmet_phi_corr_"+wstring, "TVector2::Phi_mpi_pi(deepmet_corr_{WS}.Phi())".format(WS=wstring))
        WJetsSamp.Define("mT_Deep_corr_"+wstring, "calMT_fromMET_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_{WS}, deepmet_phi_corr_{WS})".format(WS=wstring))
        sampMan.DefineAll("deepmet_pt_corr_"+wstring, "deepmet_pt", excludes=['WJets'])
        sampMan.DefineAll("deepmet_phi_corr_"+wstring, "deepmet_phi", excludes=['WJets'])
        sampMan.DefineAll("mT_Deep_corr_"+wstring, "mT_Deep", excludes=['WJets'])
        sampMan.DefineAll("u_reco_pt_Deep_corr_"+wstring, "calUT_PtPhi(mu_pt[0], mu_phi[0], deepmet_pt_corr_{WS},  deepmet_phi_corr_{WS})".format(WS=wstring))

    varyTTbarQCDScale("mur_1_muf_0p5")
    varyTTbarQCDScale("mur_1_muf_2")
    varyTTbarQCDScale("mur_0p5_muf_1")
    varyTTbarQCDScale("mur_2_muf_1")
    varyTTbarQCDScale("mur_0p5_muf_0p5")
    varyTTbarQCDScale("mur_2_muf_2")

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -75., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 75, 80])
    u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])


    ymin = 1e3
    ymax = 5e9
    legendPos= [0.92, 0.88, 0.70, 0.62]
    sampMan.cacheDraw("Mu_pt", "h_wjets_mu_pt", 35, 30, 100, DrawConfig(xmin=30, xmax=100, xlabel='p_{T}^{#mu} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))

    sampMan.cacheDraw("met_pt", "h_wjets_met_pt", 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='PF MET [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
    sampMan.cacheDraw("puppimet_pt", "h_wjets_puppimet_pt", 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Puppi MET [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
    sampMan.cacheDraw("deepmet_pt", "h_wjets_deepmet_pt", 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))

    sampMan.cacheDraw("met_phi", "h_wjets_met_phi", 20, -ROOT.TMath.Pi(), ROOT.TMath.Pi(), DrawConfig(xmin=-ROOT.TMath.Pi(), xmax=ROOT.TMath.Pi(), xlabel='PF MET #phi', ymin=10., ymax=1e12))
    sampMan.cacheDraw("puppimet_phi", "h_wjets_puppimet_phi", 20, -ROOT.TMath.Pi(), ROOT.TMath.Pi(), DrawConfig(xmin=-ROOT.TMath.Pi(), xmax=ROOT.TMath.Pi(), xlabel='Puppi MET #phi', ymin=10., ymax=1e12))
    sampMan.cacheDraw("deepmet_phi", "h_wjets_deepmet_phi", 20, -ROOT.TMath.Pi(), ROOT.TMath.Pi(), DrawConfig(xmin=-ROOT.TMath.Pi(), xmax=ROOT.TMath.Pi(), xlabel='Deep MET #phi', ymin=10., ymax=1e12))

    sampMan.cacheDraw("mT_PF", "h_wjets_mT_PF", 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='PF m_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
    sampMan.cacheDraw("mT_Puppi", "h_wjets_mT_Puppi", 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Puppi m_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
    sampMan.cacheDraw("mT_Deep", "h_wjets_mT_Deep", 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep m_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))

    sampMan.cacheDraw("u_reco_pt_PF", "h_wjets_u_reco_pt_PF", 28, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='q_{T} PF [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
    sampMan.cacheDraw("u_reco_pt_Puppi", "h_wjets_u_reco_pt_Puppi", 28, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='q_{T} Puppi [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
    sampMan.cacheDraw("u_reco_pt_Deep", "h_wjets_u_reco_pt_Deep", 28, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep q_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))

    # corrected deepmet
    def DrawCorrection(postfix):
        sampMan.cacheDraw("deepmet_pt_corr_"+postfix, "h_wjets_deepmet_pt_corr_"+postfix, 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
        sampMan.cacheDraw("mT_Deep_corr_"+postfix, "h_wjets_mT_Deep_corr_"+postfix, 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep m_{T} Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))
        sampMan.cacheDraw("u_reco_pt_Deep_corr_"+postfix, "h_wjets_u_reco_pt_Deep_corr_"+postfix, 28, 0, 140,  DrawConfig(xmin=0, xmax=140, xlabel='Deep q_{T} Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos))

    DrawCorrection("central")
    DrawCorrection("central_MG")
    DrawCorrection("central_GKS")
    DrawCorrection("mur_1_muf_0p5")
    DrawCorrection("mur_1_muf_2")
    DrawCorrection("mur_0p5_muf_1")
    DrawCorrection("mur_2_muf_1")
    DrawCorrection("mur_0p5_muf_0p5")
    DrawCorrection("mur_2_muf_2")

    # vary the QCD Scale in WJets and compare the variations
    def CompareWJetsQCDScale(postfix, fixxsec = True):
        QCDVariations = {"central": 0, "mur_1_muf_2":1, "mur_1_muf_0p5": 2, "mur_2_muf_1": 3, "mur_2_muf_2": 4, "mur_0p5_muf_1": 6, "mur_0p5_muf_0p5": 8}
        weightname = "weight_wjetsqcdscale_"+postfix
        if fixxsec:
            # use the pre-calculated weights
            WJetsSamp.Define(weightname, "weight_"+postfix+"_scaled")
        else:
            WJetsSamp.Define(weightname, "EventWeights[%d]/EventWeights[0] * weight"%QCDVariations[postfix])
        sampMan.DefineAll(weightname, "weight", excludes=['WJets'])
        
        sampMan.cacheDraw("deepmet_pt", "h_wjets_deepmet_pt_wjetsqcdscale_"+postfix, 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("mT_Deep", "h_wjets_mT_Deep_wjetsqcdscale_"+postfix, 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep m_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("u_reco_pt_Deep", "h_wjets_u_reco_pt_Deep_wjetsqcdscale_"+postfix, 28, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep q_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("deepmet_pt_corr_central", "h_wjets_deepmet_pt_corr_central_wjetsqcdscale_"+postfix, 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("mT_Deep_corr_central", "h_wjets_mT_Deep_corr_central_wjetsqcdscale_"+postfix, 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep m_{T} Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("u_reco_pt_Deep_corr_central", "h_wjets_u_reco_pt_Deep_corr_central_wjetsqcdscale_"+postfix, 28, 0, 140,  DrawConfig(xmin=0, xmax=140, xlabel='Deep q_{T} Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)

    CompareWJetsQCDScale("central")
    CompareWJetsQCDScale("mur_1_muf_2")
    CompareWJetsQCDScale("mur_1_muf_0p5")
    CompareWJetsQCDScale("mur_2_muf_1")
    CompareWJetsQCDScale("mur_2_muf_2")
    CompareWJetsQCDScale("mur_0p5_muf_1")
    CompareWJetsQCDScale("mur_0p5_muf_0p5")

    def CompareWeights(weightname, postfix):
        sampMan.cacheDraw("deepmet_pt", "h_wjets_deepmet_pt_"+postfix, 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("mT_Deep", "h_wjets_mT_Deep_"+postfix, 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep m_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("u_reco_pt_Deep", "h_wjets_u_reco_pt_Deep_"+postfix, 28, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep q_{T} [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("deepmet_pt_corr_central", "h_wjets_deepmet_pt_corr_central_"+postfix, 50, 0, 150, DrawConfig(xmin=0, xmax=150, xlabel='Deep MET Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("mT_Deep_corr_central", "h_wjets_mT_Deep_corr_central_"+postfix, 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel='Deep m_{T} Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)
        sampMan.cacheDraw("u_reco_pt_Deep_corr_central", "h_wjets_u_reco_pt_Deep_corr_central_"+postfix, 28, 0, 140,  DrawConfig(xmin=0, xmax=140, xlabel='Deep q_{T} Corrected [GeV]', ymin=ymin, ymax=ymax, legendPos=legendPos), weightname=weightname)

    CompareWeights("weight_Iso0", "MuIso0")
    CompareWeights("weight_Iso1", "MuIso1")
    CompareWeights("weight_Iso2", "MuIso2")

    sampMan.launchDraw()

    def compMCRatios(hnames, colors, linestyles, legends, outputtag, ymin=0.8, ymax=1.2):
        # compare MC with MC
        histos_met =  [THStack2TH1(sampMan.hsmcs["h_wjets_deepmet_pt_{}".format(hname)    ], "MCRatio") for hname in hnames]
        histos_u_pt = [THStack2TH1(sampMan.hsmcs["h_wjets_u_reco_pt_Deep_{}".format(hname)], "MCRatio") for hname in hnames]
        histos_mt   = [THStack2TH1(sampMan.hsmcs["h_wjets_mT_Deep_{}".format(hname)       ], "MCRatio") for hname in hnames]

        histos_met_base  = histos_met[0].Clone( histos_met[0].GetName()+"_Cloned")
        histos_u_pt_base = histos_u_pt[0].Clone(histos_u_pt[0].GetName()+"_Cloned")
        histos_mt_base   = histos_mt[0].Clone(  histos_mt[0].GetName()+"_Cloned")

        for hmet, hupt, hmt in zip(histos_met, histos_u_pt, histos_mt):
            hmet.Divide(  histos_met_base )
            hmet.SetFillColor(0)
            hupt.Divide( histos_u_pt_base )
            hupt.SetFillColor(0)
            hmt.Divide(  histos_mt_base )
            hmt.SetFillColor(0)

        if len(hnames)>5:
            legendNCols = 2
            legendPos=[0.20, 0.88, 0.70, 0.74]
        elif len(hnames)>8:
            legendNCols = 3
            legendPos=[0.15, 0.99, 0.92, 0.72]
        else:
            legendNCols = 1
            legendPos=[0.20, 0.88, 0.35, 0.74]

        if "corr" in hnames[0]:
            metlabel = "Deep MET Corrected [GeV]"
            qtlabel  = "Deep q_{T} Corrected [GeV]"
            mtlabel  = "Deep m_{T} Corrected [GeV]"
        else:
            metlabel = "Deep MET [GeV]"
            qtlabel  = "Deep q_{T} [GeV]"
            mtlabel  = "Deep m_{T} [GeV]"
        DrawHistos(histos_met, legends, 0., 150., metlabel, ymin, ymax, "Ratio", "h_wjets_MCratio_deepmet_pt_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_u_pt, legends, 0., 140., qtlabel, ymin, ymax, "Ratio", "h_wjets_MCratio_u_reco_pt_Deep_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
        DrawHistos(histos_mt, legends, 0., 140., mtlabel, ymin, ymax, "Ratio", "h_wjets_MCratio_mT_Deep_{}".format(outputtag), dology = False, drawashist=True, mycolors=colors, linestyles=linestyles, legendPos=legendPos, legendNCols=legendNCols)
    

    colors = [1,2]
    linestyles = [1,1]
    legends = ["amc@NLO", "MadGraph"]
    compMCRatios(["corr_central", "corr_central_MG"], colors, linestyles, legends, "MCSample", ymin=0.9, ymax=1.1)

    legends = ["2-Gaussian Fit", "Gaussian Smoothing"]
    compMCRatios(["corr_central", "corr_central_GKS"], colors, linestyles, legends, "FitVSSmoothing", ymin=0.9, ymax=1.1)

    legends = ["#mu_{r}=1, #mu_{f}=1", "#mu_{r}=1, #mu_{f}=0.5", "#mu_{r}=1, #mu_{f}=2", "#mu_{r}=0.5, #mu_{f}=1", "#mu_{r}=2, #mu_{f}=1", "#mu_{r}=0.5, #mu_{f}=0.5", "#mu_{r}=2, #mu_{f}=2"]
    postfixs = ["corr_central", "corr_mur_1_muf_0p5", "corr_mur_1_muf_2", "corr_mur_0p5_muf_1", "corr_mur_2_muf_1", "corr_mur_0p5_muf_0p5", "corr_mur_2_muf_2"]
    colors     = [1,2,3,4,5,6,7]
    linestyles = [1,1,1,1,1,1,1]
    compMCRatios(postfixs, colors, linestyles, legends, "QCDScale", ymin=0.9, ymax=1.1)

    colors=      [1,2,3,4,5,6,7,8,9]
    linestyles = [1,1,1,1,1,1,1,1,1]
    legends = ["central", "MadGraph", "Gauss Smoothing", "#mu_{r}=1, #mu_{f}=0.5", "#mu_{r}=1, #mu_{f}=2", "#mu_{r}=0.5, #mu_{f}=1", "#mu_{r}=2, #mu_{f}=1", "#mu_{r}=0.5, #mu_{f}=0.5", "#mu_{r}=2, #mu_{f}=2"]
    postfixs = ["corr_central", "corr_central_MG", "corr_central_GKS", "corr_mur_1_muf_0p5", "corr_mur_1_muf_2", "corr_mur_0p5_muf_1", "corr_mur_2_muf_1", "corr_mur_0p5_muf_0p5", "corr_mur_2_muf_2"]
    compMCRatios(postfixs, colors, linestyles, legends, "Systematics", ymin=0.9, ymax=1.1)

    colors     = [1,2,3,4,5,6,7]
    linestyles = [1,1,1,1,1,1,1]
    legends = ["#mu_{r}=1, #mu_{f}=1", "#mu_{r}=1, #mu_{f}=0.5", "#mu_{r}=1, #mu_{f}=2", "#mu_{r}=0.5, #mu_{f}=1", "#mu_{r}=2, #mu_{f}=1", "#mu_{r}=0.5, #mu_{f}=0.5", "#mu_{r}=2, #mu_{f}=2"]
    postfixs = ["wjetsqcdscale_central", "wjetsqcdscale_mur_1_muf_0p5", "wjetsqcdscale_mur_1_muf_2", "wjetsqcdscale_mur_0p5_muf_1", "wjetsqcdscale_mur_2_muf_1", "wjetsqcdscale_mur_0p5_muf_0p5", "wjetsqcdscale_mur_2_muf_2"]
    compMCRatios(postfixs, colors, linestyles, legends, "WJetsQCDScale")

    colors     = [1,2,3,4,5,6,7]
    linestyles = [1,1,1,1,1,1,1]
    legends = ["#mu_{r}=1, #mu_{f}=1", "#mu_{r}=1, #mu_{f}=0.5", "#mu_{r}=1, #mu_{f}=2", "#mu_{r}=0.5, #mu_{f}=1", "#mu_{r}=2, #mu_{f}=1", "#mu_{r}=0.5, #mu_{f}=0.5", "#mu_{r}=2, #mu_{f}=2"]
    postfixs = ["corr_central_wjetsqcdscale_central", "corr_central_wjetsqcdscale_mur_1_muf_0p5", "corr_central_wjetsqcdscale_mur_1_muf_2", "corr_central_wjetsqcdscale_mur_0p5_muf_1", "corr_central_wjetsqcdscale_mur_2_muf_1", "corr_central_wjetsqcdscale_mur_0p5_muf_0p5", "corr_central_wjetsqcdscale_mur_2_muf_2"]
    compMCRatios(postfixs, colors, linestyles, legends, "corr_WJetsQCDScale")

    colors = [1,2,3]
    linestyles= [1,1,1]
    postfixs = ["corr_central_MuIso0", "corr_central_MuIso1", "corr_central_MuIso2"]
    legends= ["0.30<I<0.60", "0.30<I<0.45", "0.45<I<0.60"]
    legends= ["From m_{T} < 40", "From m_{T} < 30", "From m_{T} < 50"]
    compMCRatios(postfixs, colors, linestyles, legends, "corr_MuIso")

    return 

if __name__ == "__main__":
   main()
