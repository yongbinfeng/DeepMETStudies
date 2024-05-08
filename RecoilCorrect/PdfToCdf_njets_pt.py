'''
conver the pdfs from the fits to cdfs, which will be used 
in the quantile corrections.
it also allows Gaussian kernel smoothings to smooth the u1(uparal)
and u2(uperp) distributions, intead of doing fits.
'''

import ROOT
import numpy as np
import sys, os
sys.path.append("../RecoilResol/CMSPLOTS")
from myFunction import DrawHistos
from Utils.utils import getPtBins, getJetBins
from Smoother import Smoother

ROOT.gROOT.SetBatch(True)

def PdfToCdfTF(pdf, var, postfix=""):
    cdf = pdf.createCdf( ROOT.RooArgSet(var) )
    cdf_tf = cdf.asTF( ROOT.RooArgList(var) )
    cdf_tf.SetName( pdf.GetName()  + "_cdf_tf_" + postfix)
    cdf_tf.SetTitle( pdf.GetName() + "_cdf_tf_" + postfix)
    return cdf_tf

def convertPdfToCdfTF(ws, nbins, jetbin, uname, postfix=""):
    tfs = ROOT.TList()
    assert uname in ["u1", "u2"], "uname must be either u1 or u2"
    for ibin in range(nbins):
        varname = "u_{}_{}_pt{}".format(uname, jetbin, ibin)
        var = ws.var(varname)
        pdfname = "sig_{}_{}_pt{}".format(uname, jetbin, ibin)
        pdf = ws.pdf(pdfname)

        tf = PdfToCdfTF(pdf, var, postfix)
        tfs.Add( tf )
        #tfs.append( tf )
    return tfs

def HistToCdfTF(hist, var, postfix=""):
    """
    input is a root histogram,
    construct it into RooHistPdf first,
    get the cdf from RooHistPdf and
    return as the TF1 format.
    To be consistent with the roofit method.
    """
    rdh = ROOT.RooDataHist("datahist_gs_"+postfix, "datahist_gs_"+postfix, ROOT.RooArgList(var, "argdatahist_gs_"+postfix), hist)
    pdf = ROOT.RooHistPdf( "histpdf_gs_" +postfix, "histpdf_gs_" +postfix, ROOT.RooArgSet(var),                             rdh )
    cdf_tf = PdfToCdfTF(pdf, var, postfix)
    return rdh, pdf, cdf_tf

def rdhToGKS(ws, nbins, jetbin, uname, postfix=""):
    """
    extract RooDataHist from ws and
    create the Gaussian Kernel Smoother
    """
    dhs = ROOT.TList()
    gkss = ROOT.TList()
    rdhs_gks = ROOT.TList()
    pdfs_gks = ROOT.TList()
    cdfs_gks = ROOT.TList()
    assert uname in ["u1", "u2"], "uname must be either u1 or u2"
    for ibin in range(nbins):
        varname = "u_{}_{}_pt{}".format(uname, jetbin, ibin)
        var = ws.var(varname)
        rdhname = "datahist_{}_{}_pt{}".format(uname, jetbin, ibin)
        rdh = ws.data(rdhname)
        dh = rdh.createHistogram("dh_{}_{}_pt{}_{}".format(uname, jetbin, ibin, postfix), var)

        # apply Gaussian Kernel Smoother
        smo = Smoother()
        smo.histo = dh
        smo.computeSmoothHisto()
        gks = smo.smoothHisto
        #cdf = gks.GetCumulative()
        rdh_gks, pdf_gks, cdf = HistToCdfTF( gks, var, postfix)
        print(dh.GetName(), gks.GetName(), cdf.GetName())

        dhs.Add(  dh )
        gkss.Add( gks )
        rdhs_gks.Add( rdh_gks )
        pdfs_gks.Add( pdf_gks )
        cdfs_gks.Add( cdf )

        doPlot = True
        if doPlot:
            # plot the comparison between orignal and smoothed
            xtitle = "u_{#parallel} [GeV]" if uname=="u1" else "u_{#perp} [GeV]"
            DrawHistos([dh, gks], ["Original", "Smoothed"], dh.GetXaxis().GetXmin(), dh.GetXaxis().GetXmax(), xtitle, 0., 1.2*dh.GetMaximum(), "Events / GeV", "hcomp_{}_{}_pt{}_{}".format(uname, jetbin, ibin, postfix), dology=False, showpull=True, mycolors=[1, 2], drawoptions=["", "HIST"], legendoptions=["PE", "L"], ratiobase=1, doNewman=True, legendPos=[0.88, 0.84, 0.72, 0.74])
    return dhs, gkss, rdhs_gks, pdfs_gks, cdfs_gks


def main(cat = 0, opt = 0):
    
    isData = False
    isAMCNLO = False
    scaleBkg = False
    
    doPdfToCdf = False
    doGKS = False
    
    assert cat in [0, 1, 2], "cat must be either 0, 1, or 2"
    assert opt in [0, 1], "opt must be either 0 or 1"
    
    if cat == 0:
        isData = True
    elif cat == 1:
        isAMCNLO = True
    elif cat == 2:
        isData = True
        scaleBkg = True
        
    if opt == 0:
        doPdfToCdf = True
    elif opt == 1:
        doGKS = True
    
    wstrings = ["central", "mur_1_muf_0p5", "mur_1_muf_2", "mur_0p5_muf_1", "mur_2_muf_1", "mur_0p5_muf_0p5", "mur_2_muf_2"]
    WSTRING = wstrings[0]

    if scaleBkg and isData:
        WSTRING += "_bkgScale"

    if isData:
        sampname = "Data"
    elif isAMCNLO:
        sampname = "DY"

    inputname = "results/Fit/results_fit_{}_njets_pt_{}.root".format(sampname, WSTRING)

    ifile = ROOT.TFile(inputname)
    ws = ifile.Get("fit")

    ptbins = getPtBins()
    njetbins = getJetBins()

    if doPdfToCdf:
        print("start converting pdfs from fits to cdfs from ", inputname)
        ofilename = "results/Fit/fitfunctions_{}_njets_pt_{}.root".format(sampname, WSTRING)
        ofile = ROOT.TFile(ofilename, "RECREATE")
        tfs_u1 = ROOT.TList()
        tfs_u2 = ROOT.TList()
        for ijet in range(njetbins.size-1):
            nbins_pt = ptbins.size-1
            jetbin = 'njet%d'%ijet
            tfs_u1_njet = convertPdfToCdfTF(ws, nbins_pt, jetbin, "u1", sampname)
            tfs_u2_njet = convertPdfToCdfTF(ws, nbins_pt, jetbin, "u2", sampname)

            tfs_u1.Add( tfs_u1_njet )
            tfs_u2.Add( tfs_u2_njet )

        tfs_u1.Write("tfs_{}_u1_njets_pt_{}".format(sampname, WSTRING), 1)
        tfs_u2.Write("tfs_{}_u2_njets_pt_{}".format(sampname, WSTRING), 1)

        h1_ptbins_name = "h1_ptbins_{}_{}".format(sampname, WSTRING)
        h1_ptbins = ROOT.TH1F(h1_ptbins_name, h1_ptbins_name, ptbins.size-1, ptbins)
        for ipt in range(ptbins.size-1):
            h1_ptbins.SetBinContent(ipt+1, ipt)

        h1_njetbins_name = "h1_njetbins_{}_{}".format(sampname, WSTRING)
        h1_njetbins = ROOT.TH1F(h1_njetbins_name, h1_njetbins_name, njetbins.size-1, njetbins)
        for ijet in range(njetbins.size-1):
            h1_njetbins.SetBinContent(ijet+1, ijet)

        h1_ptbins.SetDirectory( ofile )
        h1_ptbins.Write()

        h1_njetbins.SetDirectory( ofile )
        h1_njetbins.Write()

        ofile.Close()

        print("write to ", ofilename)

    if doGKS:
        print("starting creating Gaussian smeared smoother and create cdfs from ", inputname)
        odir = "results/GaussSmoother"
        if not os.path.exists(odir):
            os.makedirs(odir)
        ofilename = "{}/gaussSmoother_{}_njets_pt_{}.root".format(odir, sampname, "central")
        ofile = ROOT.TFile(ofilename, "RECREATE")
        cdfs_gks_u1 = ROOT.TList()
        cdfs_gks_u2 = ROOT.TList()
        for ijet in range(njetbins.size-1):
            nbins_pt = ptbins.size-1
            jetbin = "njet{}".format(ijet)
            _, _, _, _, cdfs_gks_u1_njet = rdhToGKS(ws, nbins_pt, jetbin, "u1", sampname)
            _, _, _, _, cdfs_gks_u2_njet = rdhToGKS(ws, nbins_pt, jetbin, "u2", sampname)

            cdfs_gks_u1.Add( cdfs_gks_u1_njet )
            cdfs_gks_u2.Add( cdfs_gks_u2_njet )

        cdfs_gks_u1.Write("cdfs_{}_u1_njets_pt_{}".format(sampname, "central"), 1)
        cdfs_gks_u2.Write("cdfs_{}_u2_njets_pt_{}".format(sampname, "central"), 1)

        h1_ptbins_name = "h1_ptbins_{}_{}".format(sampname, "central")
        h1_ptbins = ROOT.TH1F(h1_ptbins_name, h1_ptbins_name, ptbins.size-1, ptbins)
        for ipt in range(ptbins.size-1):
            h1_ptbins.SetBinContent(ipt+1, ipt)

        h1_njetbins_name = "h1_njetbins_{}_{}".format(sampname, "central")
        h1_njetbins = ROOT.TH1F(h1_njetbins_name, h1_njetbins_name, njetbins.size-1, njetbins)
        for ijet in range(njetbins.size-1):
            h1_njetbins.SetBinContent(ijet+1, ijet)

        h1_ptbins.SetDirectory( ofile )
        h1_ptbins.Write()

        h1_njetbins.SetDirectory( ofile )
        h1_njetbins.Write()

        ofile.Close()



    print("finished")
    #input()

if __name__ == "__main__":
    #main(0, 0)
    #main(0, 1)
    #main(1, 0)
    main(1, 1)
    #main(2, 0)