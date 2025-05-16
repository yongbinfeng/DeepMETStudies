'''
plot the data-MC comparisons pre/post DeepMET corrections.
and the systematic uncs.
'''
from SampleManager import DrawConfig, Sample, SampleManager
import ROOT
import numpy as np
import sys
sys.path.append("../RecoilResol/CMSPLOTS")


ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(15)

deriveZptCorr = False
derivePileupCorr = False

applyZptCorr = not deriveZptCorr
applyPileupCorr = not derivePileupCorr

applyJER = True

doSnapshot = True


def main():
    print("Program start...")

    zptcorr_name = "corrections/Zpt.root"
    pileupcorr_name = "corrections/Pileup.root"
    h_zpt_corr = "h_zpt_ratio_DataVsMC"
    h_pileup_corr = "h_pileup_ratio_DataVsMC"

    suffix_zpt = "_WithZptWeight" if applyZptCorr else "_WoZptWeight"
    suffix_pileup = "_WithPileupWeight" if applyPileupCorr else "_WoPileupWeight"
    suffix = suffix_zpt + suffix_pileup

    input_data = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot/Data.root"
    input_dy = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/DY.root"
    input_ttbar = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/ttbar.root"
    input_WW2L = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/WW2L.root"
    input_WZ2L = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/WZ2L.root"
    input_ZZ2L = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/ZZ2L.root"
    input_ZZ2L2Q = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/ZZ2L2Q.root"
    input_dytau = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/DYTauTau.root"

    DataSamp = Sample(input_data, isMC=False, legend="Data",
                      name="Data", prepareVars=False, select=False)
    DYSamp = Sample(input_dy,    xsec=0,  color=5,  reweightzpt=False,
                    legend="DY", name="DY", prepareVars=False, select=False)
    TTbarSamp = Sample(input_ttbar, xsec=0,  color=46, reweightzpt=False,
                       legend="t#bar{t}", name="ttbar", prepareVars=False, select=False)
    WW2LSamp = Sample(input_WW2L,  xsec=0,  color=38, reweightzpt=False,
                      legend="WW2L",     name="WW2L",  prepareVars=False, select=False)
    WZ2LSamp = Sample(input_WZ2L,  xsec=0,  color=39, reweightzpt=False,
                      legend="WZ2L",     name="WZ2L",  prepareVars=False, select=False)
    ZZ2LSamp = Sample(input_ZZ2L,  xsec=0,  color=37, reweightzpt=False,
                      legend="ZZ2L",     name="ZZ2L",  prepareVars=False, select=False)
    ZZ2L2QSamp = Sample(input_ZZ2L2Q,  xsec=0, color=36, reweightzpt=False,
                        legend="ZZ2L2Q", name="ZZ2L2Q",  prepareVars=False, select=False)
    DYTauSamp = Sample(input_dytau, xsec=0, color=8,  reweightzpt=False,
                       legend="DY#rightarrow#tau#tau", name="DYTauTau", prepareVars=False, select=False)

    sampMan = SampleManager(
        DataSamp, [DYSamp, DYTauSamp, WW2LSamp, WZ2LSamp, TTbarSamp, ZZ2LSamp, ZZ2L2QSamp])
    sampMan.groupMCs(["WW2L", "WZ2L", "ZZ2L", "ZZ2L2Q"],
                     "Dibosons", 38, "Dibosons")

    DataSamp.Define("norm", "1")

    branches = DataSamp.rdf.GetColumnNames()
    branches = list(branches)
    print(branches)

    if applyZptCorr:
        ROOT.gROOT.ProcessLine(
            f'TFile* f_zpt = TFile::Open("{zptcorr_name}")')
        ROOT.gROOT.ProcessLine(
            f'TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("{h_zpt_corr}")')
        # make sure the normalization is not changed after reweighting
        DYSamp.Define("zptweight_bare", "ZptReWeight(Z_pt, h_zpt_ratio, 0)")
        DYSamp.Define("weight_zpt_bare", "weight * zptweight_bare")
        DYSamp.Define("dumbVal", "1")
        h_weighted = DYSamp.rdf.Histo1D(
            ("h_weighted", "h_weighted", 2, 0, 2.0), "dumbVal", "weight_zpt_bare")
        h_woweight = DYSamp.rdf.Histo1D(
            ("h_woweight", "h_woweight", 2, 0, 2.0), "dumbVal", "weight")
        zptnorm = h_weighted.Integral() / h_woweight.Integral()
        zptnorm = 1.0 / zptnorm
        print("without zpt reweighting norm: ", h_woweight.Integral())
        print("with zpt reweighting norm: ", h_weighted.Integral())
        print("Zpt reweighting norm: ", zptnorm)
        DYSamp.Define("zptweight", f"zptweight_bare * {zptnorm}")

        DYSamp.Define("weight_corr_zpt", "weight * zptweight * norm")
        TTbarSamp.Define("weight_corr_zpt", "weight * norm * 1.4")
        sampMan.DefineAll("weight_corr_zpt", "weight * norm",
                          excludes=['DY', 'ttbar'])
    else:
        sampMan.DefineAll("weight_corr_zpt", "weight * norm",
                          excludes=['ttbar'])
        TTbarSamp.Define("weight_corr_zpt", "weight * norm * 1.4")

    if applyPileupCorr:
        ROOT.gROOT.ProcessLine(
            f'TFile* f_pileup = TFile::Open("{pileupcorr_name}")')
        ROOT.gROOT.ProcessLine(
            f'TH1D* h_pileup_ratio  = (TH1D*)f_pileup->Get("{h_pileup_corr}")')
        # make sure the normalization is not changed after reweighting
        sampMan.DefineMC(
            "pileupweight_bare", "PileupReWeight(Pileup_nTrueInt, h_pileup_ratio)")
        sampMan.DefineMC("weight_pileup_bare", "weight * pileupweight_bare")
        hs_weighted = []
        hs_woweight = []
        for mc in sampMan.mcs:
            mc.Define("dumbVal_pileup", "1")
            h_weighted = mc.rdf.Histo1D(
                ("h_weighted", "h_weighted_" + mc.name, 2, 0, 2.0), "dumbVal_pileup", f"weight_pileup_bare")
            h_woweight = mc.rdf.Histo1D(
                ("h_woweight", "h_woweight_" + mc.name, 2, 0, 2.0), "dumbVal_pileup", "weight")
            hs_weighted.append(h_weighted)
            hs_woweight.append(h_woweight)

            pileupnorm = h_weighted.Integral() / h_woweight.Integral()
            pileupnorm = 1.0 / pileupnorm
            print(f"{mc.name} without pileup reweighting norm: ",
                  h_woweight.Integral())
            print(f"{mc.name} with pileup reweighting norm: ",
                  h_weighted.Integral())
            print(f"{mc.name} pileup reweighting norm: ", pileupnorm)
            mc.Define("pileupweight", f"pileupweight_bare * {pileupnorm}")
    else:
        sampMan.DefineMC("pileupweight", "1")
        
        
    branches_added = []
    if applyJER:
        # apply JER to all MCs
        # propagate the changes of jet pts to the MET
        branches_added = sampMan.ApplyJER()
        

    DataSamp.Define("pileupweight", "1")
    sampMan.DefineAll("weight_corr",
                      "weight_corr_zpt * pileupweight")

    sampMan.SetDefaultWeightName("weight_corr")

    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    # z kinematics
    zptbins = np.array([0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0,
                       48.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 200, 220, 240, 260, 280, 300])
    # sampMan.cacheDraw("Z_pt", "histo_zjets_zpt_WoZptWeight", zptbins, DrawConfig(xmin=0, xmax=180, xlabel='p^{ll}_{T} [GeV]'), weightname="weight_WoVpt")
    sampMan.cacheDraw("Z_pt", "histo_zjets_zpt" + suffix, zptbins, DrawConfig(
        xmin=0, xmax=300, xlabel='p^{ll}_{T} [GeV]', addOverflow=True, ymin=1.0, ymax=1e8))
    sampMan.cacheDraw("Z_eta", "histo_zjets_zeta" + suffix, 30, -3, 3,
                      DrawConfig(xmin=-3, xmax=3, xlabel='#eta^{ll}'))
    sampMan.cacheDraw("Z_phi", "histo_zjets_zphi" + suffix, 30, phimin, phimax, DrawConfig(
        xmin=phimin, xmax=phimax, xlabel='#phi^{ll}'))
    sampMan.cacheDraw("m_ll", "histo_zjets_mll" + suffix, 30, 80, 100,
                      DrawConfig(xmin=80, xmax=100, xlabel='m^{ll} [GeV]'))

    # muon kinematics
    sampMan.cacheDraw("leadMuon_pt", "histo_zjets_leadmuonpt" + suffix, 30,
                      20, 80, DrawConfig(xmin=20, xmax=80, xlabel='p^{l}_{T} [GeV]'))
    sampMan.cacheDraw("leadMuon_eta", "histo_zjets_leadmuoneta" + suffix,
                      30, -3, 3, DrawConfig(xmin=-3, xmax=3, xlabel='#eta^{l}'))
    sampMan.cacheDraw("leadMuon_phi", "histo_zjets_leadmuonphi" + suffix, 30, phimin,
                      phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='#phi^{l}'))
    sampMan.cacheDraw("subleadMuon_pt", "histo_zjets_subleadmuonpt" + suffix, 30,
                      20, 80, DrawConfig(xmin=20, xmax=80, xlabel='p^{l}_{T} [GeV]'))
    sampMan.cacheDraw("subleadMuon_eta", "histo_zjets_subleadmuoneta" + suffix,
                      30, -3, 3, DrawConfig(xmin=-3, xmax=3, xlabel='#eta^{l}'))
    sampMan.cacheDraw("subleadMuon_phi", "histo_zjets_subleadmuonphi" + suffix, 30,
                      phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='#phi^{l}'))

    # for pileup rweighting
    h_pileup = DYSamp.rdf.Histo1D(("h_pileup", "h_pileup", 99, 0, 99),
                                  "Pileup_nTrueInt", "weight_corr")

    # pileup histograms
    sampMan.cacheDraw("PV_npvs", "histo_zjets_npv" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='n_{PV}', yrmin=0, yrmax=2.0))
    sampMan.cacheDraw("PV_npvsGood", "histo_zjets_npvsgood" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='n_{PV}^{good}', yrmin=0, yrmax=2.0))
    sampMan.cacheDraw("fixedGridRhoFastjetAll", "histo_zjets_rho" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='#rho', yrmin=0, yrmax=2.0))
    sampMan.cacheDraw("fixedGridRhoFastjetCentral", "histo_zjets_rhocentral" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='#rho_{central}', yrmin=0, yrmax=2.0))
    sampMan.cacheDraw("fixedGridRhoFastjetCentralCalo", "histo_zjets_rhocentralcalo" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='#rho_{centralcalo}', yrmin=0, yrmax=2.0))
    sampMan.cacheDraw("fixedGridRhoFastjetCentralChargedPileUp", "histo_zjets_rhocentralchargedpileup" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='#rho_{centralchargedpileup}', yrmin=0, yrmax=2.0))
    sampMan.cacheDraw("fixedGridRhoFastjetCentralNeutral", "histo_zjets_rhocentralneutral" + suffix, 100, 0, 100,
                      DrawConfig(xmin=0, xmax=100, xlabel='#rho_{centralneutral}', yrmin=0, yrmax=2.0))

    sampMan.launchDraw()

    if deriveZptCorr:
        # DY always dominates the Zpt spectrum
        # use the Data/MC ratio to approximate the Zpt correction on DY MC
        hratio_zpt = list(
            sampMan.hratios["histo_zjets_zpt" + suffix].values())[0]
        hratio_zpt.SetName(h_zpt_corr)
        ofile = ROOT.TFile(zptcorr_name, "RECREATE")
        hratio_zpt.Write()
        ofile.Close()

    if derivePileupCorr:
        h_pileup = h_pileup.Clone("h_pileup_MC")
        h_pileup.Scale(1.0 / h_pileup.Integral())
        f_pileup_data_name = "corrections/PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root"
        f_pileup_data = ROOT.TFile(f_pileup_data_name)
        h_pileup_data = f_pileup_data.Get("pileup").Clone("h_pileup_data")
        h_pileup_data.Scale(1.0 / h_pileup_data.Integral())

        h_pileup_ratio = h_pileup_data.Clone(h_pileup_corr)
        h_pileup_ratio.Divide(h_pileup)

        ofile = ROOT.TFile(pileupcorr_name, "RECREATE")
        h_pileup_data.SetDirectory(ofile)
        h_pileup.SetDirectory(ofile)
        h_pileup_ratio.SetDirectory(ofile)
        h_pileup_data.Write()
        h_pileup.Write()
        h_pileup_ratio.Write()
        ofile.Close()

    if doSnapshot:
        # save output ntuples
        branches += [
            "weight_corr", "pileupweight", "weight_corr_zpt",
        ]
        sampMan.snapShot(
            "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withZptPileupCorrs", branches, addNorm=False, jer_variables=branches_added)

    print("Program end...")

    input()

    return


if __name__ == "__main__":
    main()
