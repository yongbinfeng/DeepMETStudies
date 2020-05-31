import ROOT
import math
import array

class Smoother:
    def __init__(self):
        self.histo = None
        self.weights = None
        self.smoothHisto = None
        self.logScale = False
        self.gausWidth = 1.

    def getSmoothedValue(self, x):
        if self.logScale and x<=0.:
            raise StandardError("ERROR: use log scale and x<=0.")
        sumw = 0.
        sumwy = 0.
        nbins = self.histo.GetNbinsX()
        for b in range(1,nbins+1):
            xi = self.histo.GetXaxis().GetBinCenter(b)
            yi = self.histo.GetBinContent(b)
            if self.logScale and xi<=0.:
                raise StandardError("ERROR: use log scale and xi<=0.")
            dx = 0.
            if self.logScale:
                dx = (math.log(x) - math.log(xi))/self.gausWidth
            else:
                dx = (x-xi)/self.gausWidth
            wi = ROOT.TMath.Gaus(dx)
            if self.weights:
                wi *= self.weights.GetBinContent(b)
            sumw += wi
            sumwy += wi*yi
        value = 0.
        if sumw>0.:
            value = sumwy/sumw
        return value

    def computeSmoothHisto(self):
        if not self.histo:
            raise StandardError("ERROR: non existing input histo")
        self.smoothHisto = self.histo.Clone(self.histo.GetName()+"_smooth")
        self.smoothHisto.__class__ = ROOT.TH1D
        self.smoothHisto.SetDirectory(0)
        nbins = self.smoothHisto.GetNbinsX()
        for b in range(1,nbins+1):
            x = self.smoothHisto.GetBinCenter(b)
            smoothedValue = self.getSmoothedValue(x)
            self.smoothHisto.SetBinContent(b,smoothedValue)

    def getContinuousSmoothHisto(self):
        if not self.histo:
            raise StandardError("ERROR: non existing input histo")
        mini = self.histo.GetXaxis().GetBinLowEdge(1)
        maxi = self.histo.GetXaxis().GetBinUpEdge(self.histo.GetNbinsX())
        if self.logScale and mini<0.:
            raise StandardError("ERROR: use log scale and min value<0")
        if mini==0.:
            mini = self.histo.GetXaxis().GetBinUpEdge(1)/10.
        nbins = 1000
        bins = []
        if self.logScale:
            dx = (math.log(maxi) - math.log(mini))/nbins
            bins = [math.exp(math.log(mini)+i*dx) for i in range(0,nbins+1)]
        else:
            dx = (maxi - mini)/nbins
            bins = [mini+i*dx for i in range(0,nbins+1)]
        smoothHisto = ROOT.TH1D(self.histo.GetName()+"_cont",self.histo.GetTitle(), nbins, array.array('f',bins))
        for b in range(1,nbins+1):
            x = smoothHisto.GetBinCenter(b)
            smoothedValue = self.getSmoothedValue(x)
            smoothHisto.SetBinContent(b,smoothedValue)
        return smoothHisto
