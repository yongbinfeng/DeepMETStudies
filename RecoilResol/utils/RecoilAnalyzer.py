from collections import OrderedDict
import sys
import os
from CMSPLOTS.myFunction import getResolution, SubtractProfiles
from utils.utils import prepVars
import ROOT


class RecoilAnalyzer(object):
    """
    produce the response and resolution of different recoil estimators.
    To make use of the 'lazy' actions in RDataFrame, for all the responses and 
    resolutions, we try to 'prepare' them all together first, and then 'get' all 
    of them. In this case the event loop only need to be run once!
    """

    def __init__(self, rdf, recoils, rdfMC=None, rdfBkg = None, name="RecoilAnalyzer", useRMS=False, weightname="weight", verbose = False):
        self.rdf = rdf
        self.rdfMC = rdfMC
        self.rdfBkg = rdfBkg
        self.recoils = recoils
        self.name = name
        self.useRMS = useRMS
        self.weightname = weightname
        self.verbose = verbose
        self.profs = OrderedDict()
        self.histo2ds_paral_diff = OrderedDict()
        self.histo2ds_perp = OrderedDict()
        self.histos1d_paral_diff = OrderedDict()
        self.histos1d_perp = OrderedDict()

    def prepareVars(self):
        # prepare the paral, perp, paral_diff, response, etc
        # vars with respect to u_GEN in RDataFrame
        # for the response and resolution calculations
        for itype in self.recoils:
            self.rdf = prepVars(
                self.rdf, "u_{RECOIL}".format(RECOIL=itype), "u_GEN")
            if self.rdfMC:
                self.rdfMC = prepVars(
                    self.rdfMC, "u_{RECOIL}".format(RECOIL=itype), "u_GEN")
            if self.rdfBkg:
                self.rdfBkg = prepVars(
                    self.rdfBkg, "u_{RECOIL}".format(RECOIL=itype), "u_GEN")

    def prepareResponses(self, xvar, xbins, option="i"):
        self.profs[xvar] = OrderedDict()
        # xbins should be numpy arrays
        nbins_x = xbins.size-1
        for itype in self.recoils:
            hparal = "h_{RECOIL}_paral_VS_{XVAR}_{NAME}".format(
                RECOIL=itype, XVAR=xvar, NAME=self.name)
            if self.verbose:
                print(hparal, hparal, nbins_x, xbins, option, xvar, "u_{RECOIL}_paral".format(RECOIL=itype))
                
            self.profs[xvar][itype] = self.rdf.Profile1D(
                (hparal, hparal, nbins_x, xbins, option), xvar, "u_{RECOIL}_paral".format(RECOIL=itype), self.weightname)
        # add GEN
        hparal_Gen = "h_GEN_VS_{XVAR}_{NAME}".format(XVAR=xvar, NAME=self.name)
        self.profs[xvar]['GEN'] = self.rdf.Profile1D(
            (hparal_Gen, hparal_Gen, nbins_x, xbins, option), xvar, "u_GEN_pt", self.weightname)

        if self.rdfMC:
            for itype in self.recoils:
                hparal = "h_{RECOIL}_paral_VS_{XVAR}_MC_{NAME}".format(
                    RECOIL=itype, XVAR=xvar, NAME=self.name)
                self.profs[xvar][itype + "_MC"] = self.rdfMC.Profile1D(
                    (hparal, hparal, nbins_x, xbins, option), xvar, "u_{RECOIL}_paral".format(RECOIL=itype), self.weightname)
            # add GEN
            hparal_Gen = "h_GEN_VS_{XVAR}_MC_{NAME}".format(
                XVAR=xvar, NAME=self.name)
            self.profs[xvar]['GEN_MC'] = self.rdfMC.Profile1D(
                (hparal_Gen, hparal_Gen, nbins_x, xbins, option), xvar, "u_GEN_pt", self.weightname)
            
        if self.rdfBkg:
            for itype in self.recoils:
                hparal = "h_{RECOIL}_paral_VS_{XVAR}_Bkg_{NAME}".format(
                    RECOIL=itype, XVAR=xvar, NAME=self.name)
                self.profs[xvar][itype + "_Bkg"] = self.rdfBkg.Profile1D(
                    (hparal, hparal, nbins_x, xbins, option), xvar, "u_{RECOIL}_paral".format(RECOIL=itype), self.weightname)
            # add GEN
            hparal_Gen = "h_GEN_VS_{XVAR}_Bkg_{NAME}".format(
                XVAR=xvar, NAME=self.name)
            self.profs[xvar]['GEN_Bkg'] = self.rdfBkg.Profile1D(
                (hparal_Gen, hparal_Gen, nbins_x, xbins, option), xvar, "u_GEN_pt", self.weightname)

    def prepareResolutions(self, xvar, xbins, nbins_y, ymin, ymax):
        # for the resolutions
        # do both paral_diff and perp
        self.histo2ds_paral_diff[xvar] = OrderedDict()
        self.histo2ds_perp[xvar] = OrderedDict()
        nbins_x = xbins.size-1
        for itype in self.recoils:
            h_paral_diff = "h_{RECOIL}_paral_diff_VS_{XVAR}_{NAME}".format(
                RECOIL=itype, XVAR=xvar, NAME=self.name)
            h_perp = "h_{RECOIL}_perp_VS_{XVAR}_{NAME}".format(
                RECOIL=itype, XVAR=xvar, NAME=self.name)

            self.histo2ds_paral_diff[xvar][itype] = self.rdf.Histo2D(
                (h_paral_diff, h_paral_diff, nbins_x, xbins, nbins_y, ymin, ymax), xvar, "u_{RECOIL}_paral_diff".format(RECOIL=itype), self.weightname)
            self.histo2ds_perp[xvar][itype] = self.rdf.Histo2D(
                (h_perp,       h_perp,       nbins_x, xbins, nbins_y, ymin, ymax), xvar, "u_{RECOIL}_perp".format(RECOIL=itype), self.weightname)

        if self.rdfMC:
            for itype in self.recoils:
                h_paral_diff = "h_{RECOIL}_paral_diff_VS_{XVAR}_MC_{NAME}".format(
                    RECOIL=itype, XVAR=xvar, NAME=self.name)
                h_perp = "h_{RECOIL}_perp_VS_{XVAR}_MC_{NAME}".format(
                    RECOIL=itype, XVAR=xvar, NAME=self.name)

                self.histo2ds_paral_diff[xvar][itype + "_MC"] = self.rdfMC.Histo2D(
                    (h_paral_diff, h_paral_diff, nbins_x, xbins, nbins_y, ymin, ymax), xvar, "u_{RECOIL}_paral_diff".format(RECOIL=itype), self.weightname)
                self.histo2ds_perp[xvar][itype + "_MC"] = self.rdfMC.Histo2D(
                    (h_perp,       h_perp,       nbins_x, xbins, nbins_y, ymin, ymax), xvar, "u_{RECOIL}_perp".format(RECOIL=itype), self.weightname)
                
        if self.rdfBkg:
            for itype in self.recoils:
                h_paral_diff = "h_{RECOIL}_paral_diff_VS_{XVAR}_Bkg_{NAME}".format(
                    RECOIL=itype, XVAR=xvar, NAME=self.name)
                h_perp = "h_{RECOIL}_perp_VS_{XVAR}_Bkg_{NAME}".format(
                    RECOIL=itype, XVAR=xvar, NAME=self.name)

                self.histo2ds_paral_diff[xvar][itype + "_Bkg"] = self.rdfBkg.Histo2D(
                    (h_paral_diff, h_paral_diff, nbins_x, xbins, nbins_y, ymin, ymax), xvar, "u_{RECOIL}_paral_diff".format(RECOIL=itype), self.weightname)
                self.histo2ds_perp[xvar][itype + "_Bkg"] = self.rdfBkg.Histo2D(
                    (h_perp,       h_perp,       nbins_x, xbins, nbins_y, ymin, ymax), xvar, "u_{RECOIL}_perp".format(RECOIL=itype), self.weightname)
                

    def getResponses(self, xvar):
        hresponses = OrderedDict()
        for itype in self.recoils:
            if self.verbose:
                print("response for ", xvar, itype)
            hnum = self.profs[xvar][itype].Clone(
                self.profs[xvar][itype].GetName() + "_numerator")
            if self.rdfBkg:
                hnum = SubtractProfiles(
                    hnum, self.profs[xvar][itype + "_Bkg"])
            hden = self.profs[xvar]['GEN'].Clone(
                self.profs[xvar]['GEN'].GetName() + "_denominator_for_" + self.profs[xvar][itype].GetName())
            if self.rdfBkg:
                hden = SubtractProfiles(
                    hden, self.profs[xvar]['GEN_Bkg'])
            
            # profile divide uncertainties seem to be different from the TH1 divide uncertainties
            # https://root.cern.ch/doc/master/classTProfile.html#a8ce37db032734f54751347ddbbbebc09
            # use the TH1 divide uncertainties for now
            # projectionX() is needed to get the TH1
            # so save the TH1 instead
            
            hnum = hnum.ProjectionX()
            hden = hden.ProjectionX()
            hresponses[itype] = hnum.Clone(self.profs[xvar][itype].GetName() + "_response")
            hresponses[itype].Divide(hden)

        if self.rdfMC:
            # mc: no need to subtract bkg
            for itype in self.recoils:
                hnum = self.profs[xvar][itype + "_MC"].Clone(
                    self.profs[xvar][itype + "_MC"].GetName() + "_numerator_MC")
                hden = self.profs[xvar]['GEN_MC'].Clone(
                    self.profs[xvar]['GEN_MC'].GetName() + "_denominator_for_" + self.profs[xvar][itype + "_MC"].GetName())
                # projectionX() is needed to get the TH1
                # so save the TH1 instead
                hnum = hnum.ProjectionX()
                hden = hden.ProjectionX()
                hresponses[itype + "_MC"] = hnum.Clone(self.profs[xvar][itype + "_MC"].GetName() + "_response_MC")
                hresponses[itype + "_MC"].Divide(hden)
                #hresponses[itype + "_MC"] = self.profs[xvar][itype + "_MC"].Clone(
                #    self.profs[xvar][itype + "_MC"].GetName() + "_response").ProjectionX()
                #hresponses[itype + "_MC"].Divide(self.profs[xvar]['GEN_MC'].ProjectionX())
        return hresponses

    def getResolutions(self, xvar):
        hresols_paral = OrderedDict()
        hresols_perp = OrderedDict()
        for itype in self.recoils:
            hparal_diff = self.histo2ds_paral_diff[xvar][itype].Clone(self.histo2ds_paral_diff[xvar][itype].GetName() + "_resol")
            if self.rdfBkg:
                print("subtracting bkg for ", xvar, itype)
                hparal_diff.Add(self.histo2ds_paral_diff[xvar][itype + "_Bkg"].GetValue(), -1.0)
            hresols_paral[itype] = getResolution(
                hparal_diff, useRMS=self.useRMS)
            
            hperp_diff = self.histo2ds_perp[xvar][itype].Clone(self.histo2ds_perp[xvar][itype].GetName() + "_resol")
            if self.rdfBkg:
                hperp_diff.Add(self.histo2ds_perp[xvar][itype + "_Bkg"].GetValue(), -1.0)
            hresols_perp[itype] = getResolution(
                hperp_diff, useRMS=self.useRMS)

        if self.rdfMC:
            for itype in self.recoils:
                hresols_paral[itype + "_MC"] = getResolution(
                    self.histo2ds_paral_diff[xvar][itype + "_MC"], useRMS=self.useRMS)
                hresols_perp[itype + "_MC"] = getResolution(
                    self.histo2ds_perp[xvar][itype + "_MC"], useRMS=self.useRMS)
        return hresols_paral, hresols_perp

    def prepareResolutions1D(self, nbins_x, xmin, xmax):
        for itype in self.recoils:
            h_paral_diff_1D = "h_{RECOIL}_paral_diff_1D_{NAME}".format(
                RECOIL=itype, NAME=self.name)
            h_perp_1D = "h_{RECOIL}_perp_1D_{NAME}".format(
                RECOIL=itype, NAME=self.name)

            self.histos1d_paral_diff[itype] = self.rdf.Histo1D(
                (h_paral_diff_1D, h_paral_diff_1D, nbins_x, xmin, xmax), "u_{RECOIL}_paral_diff".format(RECOIL=itype), self.weightname)
            self.histos1d_perp[itype] = self.rdf.Histo1D(
                (h_perp_1D,       h_perp_1D,       nbins_x, xmin, xmax), "u_{RECOIL}_perp".format(RECOIL=itype), self.weightname)

        if self.rdfMC:
            for itype in self.recoils:
                h_paral_diff_1D = "h_{RECOIL}_paral_diff_1D_MC_{NAME}".format(
                    RECOIL=itype, NAME=self.name)
                h_perp_1D = "h_{RECOIL}_perp_1D_MC_{NAME}".format(
                    RECOIL=itype, NAME=self.name)

                self.histos1d_paral_diff[itype + "_MC"] = self.rdfMC.Histo1D(
                    (h_paral_diff_1D, h_paral_diff_1D, nbins_x, xmin, xmax), "u_{RECOIL}_paral_diff".format(RECOIL=itype), self.weightname)
                self.histos1d_perp[itype + "_MC"] = self.rdfMC.Histo1D(
                    (h_perp_1D,       h_perp_1D,       nbins_x, xmin, xmax), "u_{RECOIL}_perp".format(RECOIL=itype), self.weightname)

    def getResolutions1D(self):
        vresols_paral = OrderedDict()
        vresols_perp = OrderedDict()
        import numpy as np
        probs = np.array([0.50, 0.16, 0.84])
        quants = np.array([0., 0., 0.])

        for itype in self.recoils:
            self.histos1d_paral_diff[itype].GetQuantiles(3, quants, probs)
            vresols_paral[itype] = (
                quants[0], (quants[2]-quants[1])/2.0, quants[0]-quants[1], quants[2]-quants[0])

            self.histos1d_perp[itype].GetQuantiles(3, quants, probs)
            vresols_perp[itype] = (
                quants[0], (quants[2]-quants[1])/2.0, quants[0]-quants[1], quants[2]-quants[0])

        return vresols_paral, vresols_perp

    def saveHistos(self, outname):
        if "/" in outname:
            dirname = outname[:outname.rfind("/")]
            if not os.path.exists(dirname):
                print("Creating directory: ", dirname)
                os.makedirs(dirname)
        fout = ROOT.TFile(outname, "RECREATE")
        for xvar in self.profs:
            for itype in self.profs[xvar]:
                self.profs[xvar][itype].Write()
        for xvar in self.histo2ds_paral_diff:
            for itype in self.histo2ds_paral_diff[xvar]:
                self.histo2ds_paral_diff[xvar][itype].Write()
        for xvar in self.histo2ds_perp:
            for itype in self.histo2ds_perp[xvar]:
                self.histo2ds_perp[xvar][itype].Write()
        for itype in self.histos1d_paral_diff:
            self.histos1d_paral_diff[itype].Write()
        for itype in self.histos1d_perp:
            self.histos1d_perp[itype].Write()
        fout.Close()
