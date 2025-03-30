'''
basic class to do the sample loading, filtering, defining new 
variables and plotting.
Currently customized to work best for the W/Z case. can be easily
extened to general usage.
'''
import ROOT
import numpy as np
import sys
import os
from collections import OrderedDict
sys.path.append("../RecoilResol/CMSPLOTS/")
from myFunction import DrawHistos  # noqa

MINMASS = 80
MAXMASS = 100
LEPPTMIN = 25.0
LEPETA = 2.4
LUMI = 16.7

ROOT.ROOT.EnableImplicitMT()
ROOT.gSystem.Load("Functions_cc.so")


class DrawConfig(object):
    """
    try to sync with the arguments in DrawHistos
    in /afs/cern.ch/work/y/yofeng/public/CMSPLOTS
    """

    def __init__(self, **kwargs):
        self.xmin = kwargs.get('xmin', 0.0)
        self.xmax = kwargs.get('xmax', 0.0)
        self.xlabel = kwargs.get('xlabel', 'p_{T} [GeV]')

        self.ymin = kwargs.get('ymin', 1e1)
        self.ymax = kwargs.get('ymax', 1e8)
        self.ylabel = kwargs.get('ylabel', 'Events / GeV')

        self.dologx = kwargs.get('dologx', False)
        self.dology = kwargs.get('dology', True)
        self.donormalize = kwargs.get('donormalize', False)
        self.donormalizebin = kwargs.get('donormalizebin', True)

        self.showratio = kwargs.get('showratio', True)
        self.ratiobase = kwargs.get('ratiobase', 1)
        self.hratiopanel = kwargs.get('hratiopanel', None)
        self.yrmin = kwargs.get('yrmin', 0.71)
        self.yrmax = kwargs.get('yrmax', 1.29)
        self.yrlabel = kwargs.get('yrlabel', 'Data / MC ')

        self.legends = kwargs.get('legends', [])
        self.legendPos = kwargs.get('legendPos', [0.92, 0.85, 0.70, 0.58])
        self.colors = kwargs.get('colors', [])

        self.redrawihist = kwargs.get('redrawihist', 0)

        self.addOverflow = kwargs.get('addOverflow', False)
        self.addUnderflow = kwargs.get('addUnderflow', False)

        self.outputname = kwargs.get('outputname', 'test')


class Sample(object):
    def __init__(self, inputfiles, isMC=True, xsec=1.0, color=1, reweightzpt=False, legend="", name="",
                 isZSR=True, isWSR=False, applySF=False, bjetVeto=False, prepareVars=True, select=True):
        if isMC:
            self.mcxsec = xsec
            if xsec != 0:
                self.nmcevt = self.getNMCEvt(inputfiles)
                self.normfactor = LUMI * self.mcxsec / self.nmcevt
            else:
                # use some pre-defined normalization values instead
                self.nmcevt = 0
                self.normfactor = 1.0
            self.color = color
            self.reweightzpt = reweightzpt
        else:
            self.nmcevt = 0
            self.mcxsec = 0
            self.color = 1
            self.normfactor = 1.0
            self.reweightzpt = False
        self.isMC = isMC
        self.legend = legend
        self.name = name
        self.isZSR = isZSR
        self.isWSR = isWSR
        assert not (
            self.isZSR and self.isWSR), "can not be ZSR and WSR at the same time !!"
        self.applySF = applySF
        # this variable is used in cases where the xsec is allowed to
        # free-float, or normalized to the rest of (data-MC_sum)
        self.renormalize = False
        self.renormalizefactor = 1.0
        self.bjetVeto = bjetVeto
        self.initRDF(inputfiles)
        if prepareVars:
            self.prepareVars()
        if select:
            self.select()
        else:
            self.rdf = self.rdf_org

        self._garbagerdfs = []

    def initRDF(self, inputfiles):
        # this seems to be a bug in ROOT.
        # the tree must be held by some vars, or else it would be
        # closed after this function and the rdf would crash
        self.tree = ROOT.TChain("Events")
        if inputfiles.endswith(".root"):
            self.tree.Add(inputfiles)
        else:
            print("Adding RDF from the list of files in ", inputfiles)
            for line in open(inputfiles, "r"):
                fname = line.rstrip()
                print(fname)
                self.tree.Add(fname)

        self.rdf_org = ROOT.ROOT.RDataFrame(self.tree)
        # print self.rdf_org.Count().GetValue()

    def getNMCEvtFromOneFile(self, inputfile, suffix=""):
        ifile = ROOT.TFile(inputfile)
        ttree = ifile.Get("Events")
        hist = ROOT.TH1F(f"h1_{suffix}", f"h1_{suffix}", 2, -2, 2)
        ttree.Draw(f"(genWeight >0 ? 1: -1) >> {hist.GetName()}")
        Npos = hist.GetBinContent(2)
        Nneg = hist.GetBinContent(1)
        ifile.Close()
        return Npos, Nneg

    def getNMCEvt(self, inputfiles):
        Npos = 0
        Nneg = 0
        print("count total number of MC events from: ", inputfiles)
        nfiles = 0
        if inputfiles.endswith(".root"):
            Npos, Nneg = self.getNMCEvtFromOneFile(inputfiles)
        else:
            for line in open(inputfiles, "r"):
                fname = line.rstrip()
                print(fname)
                Npos_temp, Nneg_temp = self.getNMCEvtFromOneFile(
                    fname, suffix=nfiles)
                Npos += Npos_temp
                Nneg += Nneg_temp
                nfiles += 1
        print("total number of events: {}, in which {} are positive, {} are negative".format(
            Npos+Nneg, Npos, Nneg))
        return Npos-Nneg

    def select(self):
        """ 
        Event selections.
        Not sure why but the rdf seems to have to 
        be different before and after Filter.
        So use rdf_org to hold the pre-selected one.
        """
        if self.isZSR:
            self.rdf_temp = self.rdf_org.Filter("HLT_IsoMu24 || HLT_IsoTkMu24") \
                                .Filter('nMuon > 1') \
                                .Filter("Muon_pt[0] > {} && Muon_pt[1] > {}".format(LEPPTMIN, LEPPTMIN))  \
                                .Filter("abs(Muon_eta[0]) < {} && abs(Muon_eta[1]) < {}".format(LEPETA, LEPETA))  \
                                .Filter("Muon_pfRelIso04_all[0] < 0.15 && Muon_pfRelIso04_all[1] < 0.15") \
                                .Filter("Muon_charge[0] * Muon_charge[1] < 0") \
                                .Filter("m_ll > {MINMASS} && m_ll < {MAXMASS}".format(MINMASS=MINMASS, MAXMASS=MAXMASS))

            # self.rdf_temp = self.rdf_org.Filter("m_ll > {MINMASS} && m_ll < {MAXMASS}".format(MINMASS=MINMASS, MAXMASS=MAXMASS)) \
            #                    .Filter("mu_pt[0] > {} && mu_pt[1] > {}".format(LEPPTMIN, LEPPTMIN))  \
            #                    .Filter("abs(mu_eta[0]) < {} && abs(mu_eta[1]) < {}".format(LEPETA, LEPETA))  \
            #                    .Define("lep_n", "mu_n+el_n").Filter("lep_n==2")
            #                    #.Filter("Z_pt < 40.0")
            if self.bjetVeto:
                self.rdf = self.rdf_temp.Filter("jet_CSVLoose_n<1")
            else:
                self.rdf = self.rdf_temp
        elif self.isWSR:
            self.rdf = self.rdf_org.Define(
                "lep_n", "mu_n+el_n").Filter("lep_n==1")
        else:
            self.rdf = self.rdf_org
        # print self.rdf.Count().GetValue()

    def ApplyCut(self, cutstring):
        rdf_postcut = self.rdf.Filter(cutstring)
        self._garbagerdfs.append(self.rdf)
        print("applied cut {} to Sample {}".format(cutstring, self.name))
        self.rdf = rdf_postcut

    def Define(self, varname, formula):
        self.rdf = self.rdf.Define(varname, formula)

    def prepareVars(self):
        self.rdf_org = self.rdf_org.Define("isData", "0" if self.isMC else "1")
        # define weight
        if self.isMC and self.applySF:
            self.rdf_org = self.rdf_org.Define(
                "weight_WoVpt", "( PUWeight * NLOWeight * mu_trigSF * mu_isoSF * mu_trkSF * mu_idSF )")
        elif self.isMC:
            self.rdf_org = self.rdf_org.Define(
                "weight_WoVpt", "( genWeight > 0 ? 1. : -1. )")
        else:
            self.rdf_org = self.rdf_org.Define("weight_WoVpt", "isData")

        self.rdf_org = self.rdf_org.Define(
            "jet_n", "Sum((Jet_jetId == 6) && Jet_pt > 20.0 && abs(Jet_eta)<3.0)") \
            .Define("jet_CSVLoose_n", "Sum((Jet_jetId == 6) && Jet_pt > 20.0 && abs(Jet_eta)<2.4 && Jet_btagDeepB > 0.1918)") \
            .Define("jet_CSVMedium_n", "Sum((Jet_jetId == 6) && Jet_pt > 20.0 && abs(Jet_eta)<2.4 && Jet_btagDeepB > 0.5847)") \
            .Define("jet_CSVTight_n", "Sum((Jet_jetId == 6) && Jet_pt > 20.0 && abs(Jet_eta)<2.4 && Jet_btagDeepB > 0.8767)")

        # met filters
        self.rdf_org = self.rdf_org.Define(
            "metFilters", "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_BadPFMuonDzFilter && Flag_eeBadScFilter")

        if self.isZSR:
            # make Z boson
            self.rdf_org = self.rdf_org.Define("VecZ", "VVecM(Muon_pt[0], Muon_eta[0], Muon_phi[0], Muon_mass[0], Muon_pt[1], Muon_eta[1], Muon_phi[1], Muon_mass[1])") \
                .Define("m_ll", "VecZ.M()") \
                .Define("Z_pt", "VecZ.Pt()") \
                .Define("Z_eta", "VecZ.Eta()") \
                .Define("Z_phi", "VecZ.Phi()")
            # make recoil vec
            self.rdf_org = self.rdf_org.Define("Uvec",  "UVec(Z_pt, Z_phi, DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)") \
                               .Define("u_pt",  "Uvec.Mod()")        \
                               .Define("u_phi", "Uvec.Phi()")        \
                               .Define("u1",    "u_pt * TMath::Cos(u_phi + TMath::Pi() - Z_phi)") \
                               .Define("u2",    "u_pt * TMath::Sin(u_phi + TMath::Pi() - Z_phi)") \
                # z pt reweight
            self.rdf_org = self.rdf_org.Define(
                "ZptWeight", "ZptReWeight(Z_pt, h_zpt_ratio, isData)" if self.reweightzpt else "1.")

            self.rdf_org = self.rdf_org.Define(
                "weight", "ZptWeight * weight_WoVpt")
        elif self.isWSR and self.isMC:
            # make W boson in the transverse plane
            # self.rdf_org = self.rdf_org.Define("Wvec", "VVec(mu_pt[0], mu_eta[0], mu_phi[0], mu_e[0], gen_met_pt, 0., gen_met_phi, gen_met_pt)") \
            #                           .Define("W_pt",  "Wvec.Pt()") \
            #                           .Define("W_phi", "Wvec.Phi()")
            # self.rdf_org = self.rdf_org.Define("W_pt", "trueW_pt[0]") \
            #                           .Define("W_phi", "trueW_phi[0]")
            # make recoil vec
            # self.rdf_org = self.rdf_org.Define("Uvec", "UVec(W_pt, W_phi, deepmet_pt, deepmet_phi)") \
            #                           .Define("u_pt",  "Uvec.Mod()")   \
            #                           .Define("u_phi", "Uvec.Phi()")   \
            # weights
            self.rdf_org = self.rdf_org.Define("weight", "weight_WoVpt")

        else:
            self.rdf_org = self.rdf_org.Define("weight", "weight_WoVpt")


class SampleManager(object):
    def __init__(self, data, mcs=[]):
        self.data = data
        self.mcs = mcs
        self.to_draw = OrderedDict()
        self.getNormFactors()
        self.initGroups()

        self.hdatas = {}
        self.hsmcs = {}
        self.hratios = {}

        self.weightname_default = "weight"
        self.outdir = "plots"

    def SetDefaultWeightName(self, weightname):
        self.weightname_default = weightname

    def SetOutDir(self, outdir):
        self.outdir = outdir

    def ApplyCutAll(self, cutstring):
        self.data.ApplyCut(cutstring)
        for mc in self.mcs:
            mc.ApplyCut(cutstring)

    def ApplyCutSpecificMCs(self, cutstring, sampnames=[], sampgroupnames=[]):
        for mc in self.mcs:
            if mc.name in sampnames:
                print("Apply Cut {} for sample {}".format(cutstring, mc.name))
                mc.ApplyCut(cutstring)
            elif mc.groupname in sampgroupnames:
                print("Apply Cut {} for sample {} in grouped sample {}".format(
                    cutstring, mc.name, mc.groupname))
                mc.ApplyCut(cutstring)

    def DefineAll(self, varname, formula, excludes=[], excludeGroups=[]):
        self.DefineData(varname, formula)
        self.DefineMC(varname, formula, excludes, excludeGroups)

    def DefineMC(self, varname, formula, excludes=[], excludeGroups=[]):
        for mc in self.mcs:
            if mc.name in excludes:
                print("Define {} with {}. Skip sample {}".format(
                    varname, formula, mc.name))
                continue
            if mc.groupname in excludeGroups:
                print("Define {} with {}. Skip sample {} in grouped sample {}".format(
                    varname, formula, mc.name, mc.groupname))
                continue
            mc.rdf = mc.rdf.Define(varname, formula)

    def DefineSpecificMCs(self, varname, formula, sampnames=[], sampgroupnames=[]):
        for mc in self.mcs:
            if mc.name in sampnames:
                print("Define {} with {} for sample {}".format(
                    varname, formula, mc.name))
                mc.rdf = mc.rdf.Define(varname, formula)
            elif mc.groupname in sampgroupnames:
                print("Define {} with {} for sample {} in grouped sample {}".format(
                    varname, formula, mc.name, mc.groupname))
                mc.rdf = mc.rdf.Define(varname, formula)

    def DefineData(self, varname, formula):
        self.data.rdf = self.data.rdf.Define(varname, formula)

    def getNormFactors(self):
        self.normfactors = []
        for mc in self.mcs:
            self.normfactors.append(mc.normfactor)
        print("normalization factors ", self.normfactors)

    def initGroups(self):
        # default group information is the same as the sample itself
        for mc in self.mcs:
            mc.groupname = mc.name
            mc.groupcolor = mc.color
            mc.grouplegend = mc.legend

    def groupMCs(self, mcnames_to_be_grouped, groupname, groupcolor, grouplegend, renormalize=False, renormalizefactor=1.0):
        for mc in self.mcs:
            if mc.name in mcnames_to_be_grouped:
                print("Group sample ", mc.name, " to ", groupname)
                mc.groupname = groupname
                mc.groupcolor = groupcolor
                mc.grouplegend = grouplegend
                if renormalize:
                    mc.renormalize = True
                mc.renormalizefactor = renormalizefactor
                mcnames_to_be_grouped.remove(mc.name)

        if mcnames_to_be_grouped:
            from termcolor import colored
            print(colored('some somples are not grouped yet, please check the names...: {}'.format(
                mcnames_to_be_grouped), 'red'))

    def cacheDraw(self, *args, **kwds):
        if len(args) == 6:
            self.cacheDraw_fb(*args, **kwds)
        elif len(args) == 4:
            self.cacheDraw_vb(*args, **kwds)
        else:
            raise ValueError(
                'The argument is problematic. It has to be either 4 or 6')

    def cacheDraw_fb(self, varname, hname, nbins, xmin, xmax, drawconfigs, weightname=None):
        """ 
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        if weightname is None:
            weightname = self.weightname_default
        h_data = self.data.rdf.Histo1D(
            (hname+"_data", hname, nbins, xmin, xmax), varname, weightname)
        h_mcs = []
        for imc in range(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append(mc.rdf.Histo1D(
                (hname+"_mc_"+str(imc), hname, nbins, xmin, xmax), varname, weightname))

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)

    def cacheDraw_vb(self, varname, hname, xbins, drawconfigs, weightname=None):
        """
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        if weightname is None:
            weightname = self.weightname_default
        assert isinstance(xbins, np.ndarray), "input must be np array"
        nbins = xbins.size - 1
        h_data = self.data.rdf.Histo1D(
            (hname+"_data", hname, nbins, xbins), varname, weightname)
        h_mcs = []
        for imc in range(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append(mc.rdf.Histo1D(
                (hname+"_mc_"+str(imc), hname, nbins, xbins), varname, weightname))

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)

    def _DrawPlot(self, h_data, h_mcs, drawconfigs, hname):
        legends = []

        h_data.Scale(1.0)
        h_data.SetLineColor(1)
        h_data.SetMarkerStyle(20)
        h_data.SetMarkerSize(1)
        legends.append("Data")

        hgroupedmcs = OrderedDict()
        group_to_renormalize = ""
        for imc in range(len(h_mcs)):
            # scale the MC to the xsec
            h_mcs[imc].Scale(self.normfactors[imc])
            if self.mcs[imc].renormalizefactor != 1.0:
                print("renormalize MC {} with a factor or {}".format(
                    self.mcs[imc].name, self.mcs[imc].renormalizefactor))
                h_mcs[imc].Scale(self.mcs[imc].renormalizefactor)

            groupname = self.mcs[imc].groupname
            if groupname not in hgroupedmcs:
                hgroupedmcs[groupname] = h_mcs[imc].Clone(
                    h_mcs[imc].GetName().replace("_mc_", "_grouped_"+groupname))
                hgroupedmcs[groupname].SetLineColor(self.mcs[imc].groupcolor)
                hgroupedmcs[groupname].SetFillColor(self.mcs[imc].groupcolor)

                legends.append(self.mcs[imc].grouplegend)

                if self.mcs[imc].renormalize:
                    group_to_renormalize = groupname
            else:
                # group mc already exist. Add to the histogram
                hgroupedmcs[groupname].Add(h_mcs[imc].GetValue())

        if group_to_renormalize:
            ndata = h_data.Integral(0, h_data.GetNbinsX()+1)
            nmc = 0
            for gname, ghisto in hgroupedmcs.items():
                if gname != group_to_renormalize:
                    nmc += ghisto.Integral(0, ghisto.GetNbinsX()+1)
            weight = float(ndata-nmc) / hgroupedmcs[group_to_renormalize].Integral(
                0, hgroupedmcs[group_to_renormalize].GetNbinsX()+1)
            print("Renormalize group for {}, with {} data, {} MC and weight {}".format(
                group_to_renormalize, ndata, nmc, weight))
            hgroupedmcs[group_to_renormalize].Scale(weight)

        hsname = "hs_" + drawconfigs.outputname
        hs_gmc = ROOT.THStack(hsname, hsname)
        for h_gmc in reversed(list(hgroupedmcs.values())):
            if drawconfigs.donormalizebin:
                # scale to bin width
                h_gmc.Scale(1.0, "width")
            hs_gmc.Add(h_gmc)

        if not drawconfigs.legends:
            drawconfigs.legends = legends

        if drawconfigs.outputname == "test":
            # change default outputname to the histo name
            drawconfigs.outputname = hname

        if drawconfigs.donormalizebin:
            h_data.Scale(1.0, "width")

        self.hdatas[drawconfigs.outputname] = h_data
        self.hsmcs[drawconfigs.outputname] = hs_gmc
        self.hratios[drawconfigs.outputname] = DrawHistos([h_data, hs_gmc], drawconfigs.legends, drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax=drawconfigs.yrmax,
                                                          yrmin=drawconfigs.yrmin, yrlabel=drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos=drawconfigs.legendPos, redrawihist=drawconfigs.redrawihist, outdir=self.outdir, addOverflow=drawconfigs.addOverflow, addUnderflow=drawconfigs.addUnderflow, hratiopanel=drawconfigs.hratiopanel)

    def launchDraw(self):
        """
        launch all the draw options in to_draw
        """
        for hname, plotinfo in self.to_draw.items():
            self._DrawPlot(plotinfo[0], plotinfo[1], plotinfo[2], hname)
        print("finished drawing..")
        # clear the to_draw list
        self.to_draw = OrderedDict()

    def snapShot(self, outdir, branches, addNorm=True):
        """
        write the ntuples to a root file
        """
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        branches_mc = branches.copy()
        branches_mc += ["Pileup_nTrueInt",
                        "Pileup_pudensity", "Pileup_gpudensity", "Pileup_nPU", "Pileup_sumEOOT", "Pileup_sumLOOT"]
        if addNorm:
            branches_mc += ["norm"]

        for mc in self.mcs:
            if addNorm:
                # include the normalizing factor
                mc.rdf = mc.rdf.Define(
                    "norm", "1.0 * {}".format(mc.normfactor))
            print("snapshot for ", mc.name)
            mc.rdf.Snapshot("Events", os.path.join(
                outdir, mc.name+".root"), branches_mc)

        print("snapshot for ", self.data.name)
        self.data.rdf.Snapshot("Events", os.path.join(
            outdir, self.data.name+".root"), branches)

        print("finished snapshot..")
