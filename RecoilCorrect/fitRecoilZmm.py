'''
script to run the multi-Gaussian fits on Z data/MC samples 
of u1(uparal) and u2(uperp), in bins of Z pt and jet multiplicity
'''

import ROOT
import numpy as np
from collections import OrderedDict
import sys, os
sys.path.append("../RecoilResol/CMSPLOTS/")
from tdrstyle import setTDRStyle
import CMS_lumi
from Utils.utils import getPtBins, getJetBins
import pickle

isData = True
isAMCNLO = False
isMADGRAPH = False

subtractBkg = True

assert isData+isAMCNLO+isMADGRAPH==1, "must pick one sample from data, or amc@nlo, or madgraph"

NGAUSSIAN = 2

ROOT.gROOT.SetBatch(True)

def create_standard_ratio_canvas(name="",xsize=600,ysize=600) :
    """ ported from SampleManager in WGamma Analysis. setup plot canvas with ratio inset"""
    curr_canvases = {}
    curr_canvases['base'] = ROOT.TCanvas('basecan'+name, 'basecan', xsize, ysize)

    padsize1 = 0.67
    padsize2 = 0.33
    curr_canvases['bottom'] = ROOT.TPad('bottompad'+name, 'bottompad', 0.0, 0.0, 1, padsize2)
    curr_canvases['top'] = ROOT.TPad('toppad'+name, 'toppad', 0.0, padsize2, 1, 1)
    curr_canvases['top'].SetTopMargin(0.05*0.667/(padsize1*padsize1))
    curr_canvases['top'].SetBottomMargin(0.012/padsize1)
    #curr_canvases['top'].SetLeftMargin(0.1)
    #curr_canvases['top'].SetRightMargin(0.05)
    curr_canvases['bottom'].SetTopMargin(0.010/padsize2)
    curr_canvases['bottom'].SetBottomMargin(0.13/padsize2)
    curr_canvases['bottom'].SetGridy(1)
    curr_canvases['base'].cd()
    curr_canvases['bottom'].Draw()
    curr_canvases['top'].Draw()

    return curr_canvases

def ratio_formatting(subframe, padsize2=0.33):
    """ resize ratio label sizes as they are plotted in smaller pad """
    xAxs = subframe.GetXaxis()
    yAxs = subframe.GetYaxis()
    xAxs.SetTitleSize(0.050/padsize2)
    xAxs.SetLabelSize(0.045/padsize2)
    xAxs.SetTitleOffset(1.1)
    yAxs.SetTitleSize(0.050/padsize2)
    yAxs.SetLabelSize(0.045/padsize2)
    yAxs.SetTitleOffset(1.3*padsize2)
    yAxs.CenterTitle()
    yAxs.SetNdivisions(8)
    return

def plotRecoilFit(var, datahist, pdf, nGaussians, nParams, fitresult, njetmin, njetmax, ptmin, ptmax, xmean, xrms, result, postfix):
    """ 
    plot fit result.
    the function names pre-defined in doFit.
    This function needs to be modified if used in other cases.
    """
    frame = var.frame()
    datahist.plotOn(frame,
                    ROOT.RooFit.MarkerStyle(ROOT.kFullCircle),
                    ROOT.RooFit.MarkerSize(0.8),
                    ROOT.RooFit.DrawOption("ZP")
                   )
    pdf.plotOn(     frame, 
                    ROOT.RooFit.FillColor(7), 
                    ROOT.RooFit.VisualizeError(fitresult, 1), 
                    ROOT.RooFit.DrawOption("F")
              )
    # redraw datahist
    # funny 'feature' in roofit: the datahist has to be drawn before pdf,
    # otherwise the normalization would be an issue in the rooplot.
    # so draw datathist first, then pdf, then redraw datahist
    datahist.plotOn(frame,
                    ROOT.RooFit.MarkerStyle(ROOT.kFullCircle),
                    ROOT.RooFit.MarkerSize(0.8),
                    ROOT.RooFit.DrawOption("ZP")
                   )
    # plot components
    colors = [ ROOT.kRed, ROOT.kGreen, ROOT.kCyan, ROOT.kOrange, ROOT.kPink]
    for iG in range(nGaussians):
        post = "{}_{}".format(iG, postfix)
        pdf.plotOn(     frame, 
                        ROOT.RooFit.Components("gauss"+post), 
                        ROOT.RooFit.LineStyle( ROOT.kDashed ),
                        ROOT.RooFit.LineColor( colors[iG] )
                  )
    # again, redraw the sig pdf such that the pulls in the following
    # is correct. This is also due to the 'feature' in RooFit...
    pdf.plotOn(     frame,
                    ROOT.RooFit.LineColor(ROOT.kBlue)
              )

    setTDRStyle()
    CMS_lumi.lumi_sqrtS = "(13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Internal"

    curr_canvases = create_standard_ratio_canvas(name=postfix)
    curr_canvases['top'].cd()
    padsize1 = 0.67

    label = "u1" if "u1" in postfix else "u2"
    xtitle = "u_{#parallel} [GeV]" if label=='u1' else 'u_{#perp } [GeV]'

    frame.GetXaxis().SetTitle( xtitle )
    frame.GetYaxis().SetTitle("A.U.")
    frame.GetYaxis().SetTitleSize(0.050/padsize1)
    frame.GetYaxis().SetLabelSize(0.045/padsize1)
    frame.GetXaxis().SetTitleSize(0.050/padsize1)
    frame.GetXaxis().SetLabelSize(0.045/padsize1)
    frame.GetYaxis().SetTitleOffset(1.3*padsize1)
    frame.GetXaxis().SetTitleOffset(1.1*padsize1)

    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleSize(0)
    frame.Draw()

    chi2 = frame.chiSquare( nParams )
    latex = ROOT.TLatex()
    latex.SetTextSize(0.027/padsize1)
    latex.SetTextFont(42)
    latex.SetTextColor(1)
    latex.DrawLatexNDC(0.20, 0.85-0.05*0, "%.1f< Njets <%.1f"%(njetmin, njetmax))
    latex.DrawLatexNDC(0.20, 0.85-0.05*1, "%.1f< p^{Z}_{T} <%.1f"%(ptmin, ptmax))
    #latex.DrawLatexNDC(0.20, 0.85-0.05*2, "#chi^{2}/ndof = %.2f"%chi2)
    latex.DrawLatexNDC(0.20, 0.85-0.05*2, "Mean = %.1f"%xmean)
    latex.DrawLatexNDC(0.20, 0.85-0.05*3, "RMS = %.2f"%xrms)
    icol = 0
    for iG in range(nGaussians):
        latex.DrawLatexNDC(0.72, 0.85-0.05*icol, "#mu_{%d} = %.2f #pm %.2f"%(   iG, result['mean'+str(iG) ][0],    result['mean'+str(iG) ][1] ) )
        icol += 1
        latex.DrawLatexNDC(0.72, 0.85-0.05*icol, "#sigma_{%d} = %.2f #pm %.2f"%(iG, result['sigma'+str(iG)][0],    result['sigma'+str(iG)][1] ) )
        icol += 1
        if iG>0:
            latex.DrawLatexNDC(0.72, 0.85-0.05*icol, "f_{%d} = %.2f #pm %.2f"%( iG, result['frac'+str(iG) ][0],    result['frac'+str(iG) ][1] ) )
            icol += 1
    frame.addObject( latex )

    CMS_lumi.lumi_13TeV = "16.8 fb^{-1}" if isData else ""
    CMS_lumi.CMS_lumi( curr_canvases['top'], 4, 0)

    hpull = frame.pullHist()
    subframe = var.frame()
    subframe.addPlotable(hpull, "P")
    subframe.SetMaximum(5.0)
    subframe.SetMinimum(-5.0)
    subframe.GetYaxis().SetTitle("Pull")
    subframe.GetXaxis().SetTitle(xtitle)
    curr_canvases['bottom'].cd()
    subframe.Draw()
    ratio_formatting( subframe )

    if isData:
        tag = "Data"
    elif isAMCNLO:
        tag = "MC"
    elif isMADGRAPH:
        tag = "MG_MC"
    curr_canvases['base'].Print("plots/fit_{}_{}_{}Gauss.pdf".format(tag, postfix, nGaussians))
    curr_canvases['base'].Print("plots/fit_{}_{}_{}Gauss.png".format(tag, postfix, nGaussians))

    frame.SetMaximum(frame.GetMaximum()*10.0)
    frame.SetMinimum(1e-1)
    curr_canvases['top'].SetLogy()
    curr_canvases['base'].Print("plots/fit_{}_{}_{}Gauss_logy.pdf".format(tag, postfix, nGaussians))

    curr_canvases['base'].Close()

    return chi2


def doFit(hist, njetmin, njetmax, ptmin, ptmax, ws, postfix="u1_pt0", nGaussians = 3):
    """ do the fit on u1(uparal) or u2(uperp) """ 

    assert "u1" in postfix or "u2" in postfix, "the fit has to be on either u1 or u2 !"

    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()
    xmean = hist.GetMean()
    xrms  = hist.GetRMS()
    xbins = hist.GetNbinsX()

    var = ROOT.RooRealVar( "u_"+postfix, "u_"+postfix, xmin, xmax, "GeV")
    var.setBins(xbins)
    datahist = ROOT.RooDataHist("datahist_"+postfix, "datahist_"+postfix, ROOT.RooArgList(var, "argdatahist_"+postfix), hist)

    # parameters to be fit
    means = []
    sigmas = []
    fracs = []
    # gaussians
    gaussians = []
    constGaussians = []

    # intial value and range of parameters
    # these are from tests
    sigma_initials = []
    sigma_initials.append( (0.7*xrms,  0.,       1.0*xrms)  ) # sigma0
    sigma_initials.append( (1.4*xrms,  1.0*xrms, 4.5*xrms)  ) # sigma1
    sigma_initials.append( (5.0*xrms,  4.0*xrms, 8.0*xrms)  ) # sigma2
    sigma_initials.append( (8.0*xrms,  6.0*xrms, 10.0*xrms) ) # sigam3
    sigma_initials.append( (12.0*xrms, 9.0*xrms, 20.0*xrms) ) # sigma4

    frac_initials = []
    frac_initials.append( (0.26,  0.10, 0.50) ) # frac1
    frac_initials.append( (0.10,  0.0,  0.20) ) # frac2
    frac_initials.append( (0.05,  0.0,  0.20) ) # frac3
    frac_initials.append( (0.05,  0.0,  0.20) ) # frac4
    for iG in range(nGaussians):
        post = "{}_{}".format(iG, postfix)
        means.append(  ROOT.RooRealVar("mean" +post,  "mean" +post,   xmean,      xmin-50.,   xmax+50.,   "GeV") )
        sigmas.append( ROOT.RooRealVar("sigma"+post,  "sigma"+post,   sigma_initials[iG][0], sigma_initials[iG][1],  sigma_initials[iG][2],  "GeV") )

        gaussians.append( ROOT.RooGaussian("gauss"+post, "gauss"+post, var, means[iG], sigmas[iG]) )
        constGaussians.append(ROOT.RooGaussian("constGauss"+post, "constGauss"+post, means[iG], ROOT.RooFit.RooConst(xmean), ROOT.RooFit.RooConst(0.15*xrms)))

        if iG>0:
            # number of fractions is 1 less than the others: 
            # frac_n * Gaussian_n + frac_n-1 * Gaussian_n-1 + ... + frac_2*Gaussian_2 + frac_1*Gaussian_1 + (1-frac_n-...-frac_1)*Gaussian_0
            fracs.append( ROOT.RooRealVar("frac"+post, "frac"+post,    frac_initials[iG-1][0], frac_initials[iG-1][1],    frac_initials[iG-1][2],   "GeV") )

    # nGaussians mean + nGaussians sigma + nGaussians-1 frac
    nParams = nGaussians + nGaussians + nGaussians-1
    #if "u2" in postfix:
    #    # set the mean always to 0 for u2
    #    for mean in means:
    #        mean.setVal(0)
    #        mean.setConstant(1)
    #        nParams -= 1

    # construct fit model
    rshapes = ROOT.RooArgList()
    for igaussian in reversed(gaussians):
        rshapes.add(igaussian)

    rfracs = ROOT.RooArgList()
    for ifrac in reversed(fracs):
        rfracs.add(ifrac)

    sig = ROOT.RooAddPdf("sig_"+postfix, "sig_"+postfix, rshapes, rfracs)

    # run fit
    fitresult = sig.fitTo(
                            datahist, 
                            #ROOT.RooFit.SumW2Error(True),
                            ROOT.RooFit.NumCPU(4), 
                            ROOT.RooFit.Minimizer("Minuit2","minimize"), 
                            #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss0) ), 
                            #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss2) ),
                            ROOT.RooFit.Minos(),
                            ROOT.RooFit.Strategy(2),
                            ROOT.RooFit.Save()
                          )

    nTries = 0
    #while ( ( fitresult.status()>0 or fitresult.covQual()<3 ) and nTries<10 ):
    while ( fitresult.status()>0 and nTries<1 ):
        fitresult = sig.fitTo(
                                datahist,
                                #ROOT.RooFit.SumW2Error(True),
                                ROOT.RooFit.NumCPU(4),
                                ROOT.RooFit.Minimizer("Minuit2","scan"),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss0) ),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss2) ),
                                ROOT.RooFit.Minos(),
                                ROOT.RooFit.Strategy(2),
                                ROOT.RooFit.Save()
                            )
        fitresult = sig.fitTo(
                                datahist,
                                #ROOT.RooFit.SumW2Error(True),
                                ROOT.RooFit.NumCPU(4),
                                ROOT.RooFit.Minimizer("Minuit2","migrad"),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss0) ),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss2) ),
                                ROOT.RooFit.Hesse(),
                                ROOT.RooFit.Strategy(2),
                                ROOT.RooFit.Save()
                            )
        fitresult = sig.fitTo(
                                datahist,
                                #ROOT.RooFit.SumW2Error(True),
                                ROOT.RooFit.NumCPU(4),
                                ROOT.RooFit.Minimizer("Minuit2","improve"),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss0) ),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss2) ),
                                ROOT.RooFit.Minos(),
                                ROOT.RooFit.Strategy(2),
                                ROOT.RooFit.Save()
                            )
        fitresult = sig.fitTo(
                                datahist,
                                #ROOT.RooFit.SumW2Error(True),
                                ROOT.RooFit.NumCPU(4),
                                ROOT.RooFit.Minimizer("Minuit2","minimize"),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss0) ),
                                #ROOT.RooFit.ExternalConstraints( ROOT.RooArgSet(constGauss2) ),
                                ROOT.RooFit.Minos(),
                                ROOT.RooFit.Strategy(2),
                                ROOT.RooFit.Save()
                            )
        nTries += 1

    result = OrderedDict()
    for iG in range(nGaussians):
        result['mean'+str(iG) ]  = (means[iG].getVal(),  means[iG].getError()  )
        result['sigma'+str(iG)]  = (sigmas[iG].getVal(), sigmas[iG].getError() )
        if iG>0:
            result['frac'+str(iG)] = (fracs[iG-1].getVal(), fracs[iG-1].getError() )

    # plot the fit
    chi2 = plotRecoilFit(var, datahist, sig, nGaussians, nParams, fitresult, njetmin, njetmax, ptmin, ptmax, xmean, xrms, result, postfix)

    result['chi2']   = (chi2,            0.               )

    getattr( ws, "import" ) ( sig )
    getattr( ws, "import" ) ( datahist )

    return result


def main():
    
    if isData:
        sampname = "Data"
    elif isAMCNLO:
        sampname = "ZJets_NLO"
    elif isMADGRAPH:
        sampname = "ZJets_MG"

    haveScaleVariations = {
                            "Data"     : 0,
                            "ZJets_NLO": 0,
                            "ZJets_MG" : 0,
                            "TTbar"    : 0,
                            "WW2L"     : 0,
                            "WZ2L"     : 0,
                            "ZZ2L"     : 0,
                            "ZZ2L2Q"   : 0,
                          }

    wnames = {
                1: "weight_mur_1_muf_1",
                0: "central",
              }

    ptbins = getPtBins()
    njetbins = getJetBins()
            
    print(ptbins)
    print(njetbins)

    #ptbins = np.array([0,0.5,1.0])
    #njetbins = np.array([-0.5,0.5])

    # read u1 (uparal) and u2 (uperp) distributions 
    # in differen Z pt and jet multiplicity bins
    def readU1U2(sampname, wname):
        postfix="{}_{}".format(sampname, wname)
        filename = "results/U1U2/histos_u1u2_{}_njets_pt_{}.root".format(sampname, wname)

        print("read u1, u2 histograms from file ", filename)

        hfile = ROOT.TFile(filename)
        histos_u1 = OrderedDict() 
        histos_u2 = OrderedDict()
        for ijet in range(njetbins.size-1):
            njetmin = njetbins[ijet]
            njetmax = njetbins[ijet+1]
            histos_u1[(njetmin, njetmax)] = OrderedDict()
            histos_u2[(njetmin, njetmax)] = OrderedDict()
            for ipt in range(ptbins.size-1):
                ptmin = ptbins[ipt]
                ptmax = ptbins[ipt+1]

                wstring = "njetbin_{}_ptbin_{}_{}".format(ijet, ipt, postfix)
                print(wstring)
                hname_u1 = "hist_uparal_{}".format(wstring)
                hname_u2 = "hist_uperp_{}".format(wstring)
                histos_u1[(njetmin, njetmax)][(ptmin, ptmax)] = hfile.Get( hname_u1 ).Clone(hname_u1.replace("hist_", "h_"))
                histos_u2[(njetmin, njetmax)][(ptmin, ptmax)] = hfile.Get( hname_u2 ).Clone(hname_u2.replace("hist_", "h_"))
                histos_u1[(njetmin, njetmax)][(ptmin, ptmax)].SetDirectory(0)
                histos_u2[(njetmin, njetmax)][(ptmin, ptmax)].SetDirectory(0)

        return histos_u1, histos_u2

    histos_u1, histos_u2 = readU1U2(sampname, wnames[haveScaleVariations[sampname]])

    if isData and subtractBkg:
        # read ttbar bkg template and subtract from data
        sampnames = ["TTbar", "WW2L", "WZ2L", "ZZ2L", "ZZ2L2Q"]

        for sampBkgname in sampnames:
            histosBkg_u1, histosBkg_u2 = readU1U2(sampBkgname, wnames[haveScaleVariations[sampBkgname]])

            for ijetbin in list(histos_u1.keys()):
                for iptbin in list(histos_u1[ijetbin].keys()):
                    h_u1 = histos_u1[ijetbin][iptbin]
                    h_u2 = histos_u2[ijetbin][iptbin]
                    hBkg_u1 = histosBkg_u1[ijetbin][iptbin]
                    hBkg_u2 = histosBkg_u2[ijetbin][iptbin]

                    # subtract ttbar bkg contribution from data
                    h_u1.Add(hBkg_u1, -1.0)
                    h_u2.Add(hBkg_u2, -1.0)
        
        print("finished subtracting ttbar background")

    ws = ROOT.RooWorkspace('fit')

    results_u1 = OrderedDict()
    results_u2 = OrderedDict()
    ijet = 0
    for ijetbin in list(histos_u1.keys()):
        ipt = 0
        results_u1[ijetbin] = OrderedDict()
        results_u2[ijetbin] = OrderedDict()
        for iptbin in list(histos_u1[ijetbin].keys()):
            results_u1[ijetbin][iptbin] = doFit( histos_u1[ijetbin][iptbin], ijetbin[0], ijetbin[1], iptbin[0], iptbin[1], ws, postfix="u1_njet{}_pt{}".format(ijet, ipt), nGaussians=NGAUSSIAN )
            results_u2[ijetbin][iptbin] = doFit( histos_u2[ijetbin][iptbin], ijetbin[0], ijetbin[1], iptbin[0], iptbin[1], ws, postfix="u2_njet{}_pt{}".format(ijet, ipt), nGaussians=NGAUSSIAN )
            ipt += 1
        ijet += 1

    print(results_u1)
    print(results_u2)


    doDump = True
    if doDump:
        tag = "{}Gauss".format(NGAUSSIAN)
        outdir = "results/Fit/"
        if os.path.exists(outdir)==False:
            os.makedirs(outdir)
        ofilename = "results/Fit/results_{}_njets_pt_{}.pickle".format(sampname, "central")

        print("write to pickle")
        ofile = open( ofilename, 'wb' )
        pickle.dump( [results_u1, results_u2], ofile )
        ofile.close()

        print("write to workspace")
        ws.writeToFile( ofilename.replace("results_", "results_fit_").replace("pickle", "root") )

    print("Finished..")

    return 

if __name__ == "__main__":
   main()
