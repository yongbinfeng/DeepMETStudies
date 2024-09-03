import ROOT
import sys, os
sys.path.append("../RecoilResol/CMSPLOTS")
sys.path.append("../RecoilResol/")
from CMSPLOTS.myFunction import DrawHistos
from collections import OrderedDict
from utils.utils import getnVtxBins
import argparse

noLumi = False
writeOutput = False
#nPV = "PV_npvsGood"
nPV = "PV_npvs"

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("Functions_cc.so")

outdir = f"plots/MET_xy"
if not os.path.exists(outdir):
   print(f"Creating {outdir}")
   os.makedirs(outdir)

chain = ROOT.TChain("Events")
chain.Add("/afs/cern.ch/work/y/yofeng/public/outputroot/Data.root")
rdf_data = ROOT.ROOT.RDataFrame(chain)

branches_data = list(rdf_data.GetColumnNames())

chainMC = ROOT.TChain("Events")
chainMC.Add("/afs/cern.ch/work/y/yofeng/public/outputroot/DY.root")
rdf_MC = ROOT.ROOT.RDataFrame(chainMC)

branches_MC = list(rdf_MC.GetColumnNames())

def prepareVars(rdf):
   rdf = rdf.Define("pT_muons", "Z_pt").Define("phi_muons", "Z_phi")
   rdf = rdf.Define("PF_x", "MET_pt*TMath::Cos(MET_phi)") \
            .Define("PF_y", "MET_pt*TMath::Sin(MET_phi)") \
            .Define("PF_phi", "MET_phi") \
            .Define("PUPPI_x", "PuppiMET_pt*TMath::Cos(PuppiMET_phi)") \
            .Define("PUPPI_y", "PuppiMET_pt*TMath::Sin(PuppiMET_phi)") \
            .Define("PUPPI_phi", "PuppiMET_phi") \
            .Define("DeepMET_x", "DeepMETResolutionTune_pt*TMath::Cos(DeepMETResolutionTune_phi)") \
            .Define("DeepMET_y", "DeepMETResolutionTune_pt*TMath::Sin(DeepMETResolutionTune_phi)") \
            .Define("DeepMET_phi", "DeepMETResolutionTune_phi")
   return rdf

rdf_data = prepareVars(rdf_data)
rdf_MC = prepareVars(rdf_MC)


recoils = ["PF", "PUPPI", "DeepMET"]

colors = {
            "PF": 1,
            "PUPPI": 2,
            "GEN": 6,
            "DeepMET": 4,
            "DeepMETCorr": 8,
         }

labels = {
            "TK": "TK",
            "PF": "PF",
            "PUPPI": "PUPPI",
            "GEN": "GEN",
            "TKPHO": "TK+Photon",
            "DeepMET": "DeepMET",
            "DeepMETCorr": "DeepMET",
         }

markers = {
            "PF": 20,
            "PUPPI": 21,
            "GEN": 22,
            "DeepMET": 22,
            "DeepMETCorr": 23,
}

linestyles = {}

for itype in recoils:
   colors[itype + "_MC"] = colors[itype]

   linestyles[itype] = 1
   linestyles[itype + "_MC"] = 2
   

def GetLineStyles(hdict):
    return [1 if "_MC" not in itype else 2 for itype in hdict.keys()]

def GetMarkers(hdict):
   return [markers[itype] if "_MC" not in itype else 1 for itype in hdict.keys()]

def GetDrawOptions(hdict):
   return ["EP" if "_MC" not in itype else "HIST" for itype in hdict.keys()]
   
   
def deriveXYCorr(oname = "MET_xy_dataMC.root", suffix = "PreXYCorr", writeOutput = False):
   """
   derive the XY corrections for the MET
   and save the output in a root file
   
   profile the MET in x and y components
   and fit the profiles with a pol1
   as a function of nVtx
   """

   hprofs_data = OrderedDict()
   hprofs_MC = OrderedDict()
   
   hphis_data = OrderedDict()
   hphis_MC = OrderedDict()

   for ax in ["x", "y"]:
      hprofs_data[ax] = OrderedDict()
      hprofs_MC[ax] = OrderedDict()
   
   for itype in recoils:
      for ax in ["x", "y"]:
         hname = f"hprof_{itype}_{ax}_vs_nVtx_{suffix}"
         hname_data = f"{hname}_data"
         hname_MC = f"{hname}_MC"
         hprofs_data[ax][itype] = rdf_data.Profile1D((hname_data, hname_data, 70, 0, 70), nPV, itype + "_" + ax)
         hprofs_MC[ax][itype] = rdf_MC.Profile1D((hname_MC, hname_MC, 70, 0, 70), nPV, itype + "_" + ax)
         
      hphi_name = f"hphi_{itype}_{suffix}"
      hphi_name_data = f"{hphi_name}_data"
      hphi_name_MC = f"{hphi_name}_MC"
      hphis_data[itype] = rdf_data.Histo1D((hphi_name_data, f"{itype}_phi", 100, -3.2, 3.2), itype + "_phi")
      hphis_MC[itype] = rdf_MC.Histo1D((hphi_name_MC, f"{itype}_phi", 100, -3.2, 3.2), itype + "_phi")

   # fit the profiles with a pol1
   tfs_to_save = []
   for itype in recoils:
      for ax in ["x", "y"]:
         fname = f"tf1_{itype}_{ax}_vs_nVtx"
         fname_data = f"{fname}_data"
         fname_MC = f"{fname}_MC"

         tf1_data = ROOT.TF1(fname_data, "pol1", 0, 70)
         hprofs_data[ax][itype].Fit(tf1_data, "R")
         tf1_MC = ROOT.TF1(fname_MC, "pol1", 0, 70)
         hprofs_MC[ax][itype].Fit(tf1_MC, "R")
         
         hname = f"hprof_{itype}_{ax}_vs_nVtx_{suffix}"

         DrawHistos([hprofs_data[ax][itype], hprofs_MC[ax][itype], tf1_data, tf1_MC], [labels[itype] + " Data", labels[itype] + " MC"], 0, 70, "nPV", -10.0, 10.0, f"<MET_{ax}> [GeV]", hname, drawashist=True, dology=False, legendPos=[0.20, 0.70, 0.50, 0.90], mycolors=[colors[itype]]*4, linestyles = [linestyles[itype], linestyles[itype + "_MC"], 6, 6], noLumi=noLumi, outdir=outdir)

         tfs_to_save.append(tf1_data)
         tfs_to_save.append(tf1_MC)
         
      hphi_name = f"hphi_{itype}_{suffix}"
      
      DrawHistos([hphis_data[itype], hphis_MC[itype]], [labels[itype] + " Data", labels[itype] + " MC"], -3.2, 3.2, "#phi", 0, 0.027, "a.u.", hphi_name, drawashist=True, dology=False, legendPos=[0.20, 0.70, 0.50, 0.90], mycolors=[colors[itype], colors[itype + "_MC"]], linestyles = [linestyles[itype], linestyles[itype + "_MC"]], noLumi=noLumi, outdir=outdir, donormalize=True)
   
   
   
   
   h_toDraws = OrderedDict()
   for itype in recoils:
      h_toDraws[itype] = hphis_data[itype]
   for itype in recoils:
      h_toDraws[itype + "_MC"] = hphis_MC[itype]
   
   hcolors = [colors[itype] for itype in recoils]*2
   hmarkers = GetMarkers(h_toDraws)
   hlinestyles = GetLineStyles(h_toDraws)
   hdrawoptions = GetDrawOptions(h_toDraws)
   
   
   args = {
      "mycolors": hcolors,
      "markerstyles": hmarkers,
      "linestyles": hlinestyles,
      "drawoptions": hdrawoptions,
      "outdir": outdir,
      "noLumi": noLumi,
      "dology": False,
      "drawashist": False,
      "donormalize": True,
      "legendPos": [0.20, 0.70, 0.50, 0.90],
      "lheader": "Before XY correction",
   }
   
   DrawHistos(h_toDraws.values(), [labels[itype] for itype in recoils], -3.2, 3.2, "#phi", 0, 0.0299, "a.u.", f"hphi_{suffix}", **args)
            
   #DrawHistos([hphis_data[itype] for itype in recoils] + [hphis_MC[itype] for itype in recoils], [labels[itype] for itype in recoils], -3.2, 3.2, "#phi", 0, 0.027, "a.u.", f"hphi_{suffix}", drawashist=True, dology=False, legendPos=[0.20, 0.70, 0.50, 0.90], mycolors=[colors[itype] for itype in recoils]*2, linestyles = [linestyles[itype] for itype in recoils] + [linestyles[itype + "_MC"] for itpe in recoils], noLumi=noLumi, outdir=outdir, donormalize=True)

   # save the output fits
   if writeOutput:
      print(f"Saving the output in {oname}")
      ofile = ROOT.TFile(oname, "recreate")
      for itf in tfs_to_save:
         itf.Write()
      ofile.Close()
   
   return
   
   
def applyXYCorr(corrname, suffix = "PostXYCorr", writeOutput = False):
   """
   apply XY corrections to the MET 
   (by subtracting the profiled XY from the MET components)
   and save the corrected output in a new file
   """
   ROOT.gROOT.ProcessLine(f"TFile* ifile = TFile::Open(\"{corrname}\")")
   for itype in recoils:
      for ax in ["x", "y"]:
         ROOT.gROOT.ProcessLine(f"TF1* tf1_{itype}_{ax}_vs_nVtx_data = (TF1*)ifile->Get(\"tf1_{itype}_{ax}_vs_nVtx_data\")")
         ROOT.gROOT.ProcessLine(f"TF1* tf1_{itype}_{ax}_vs_nVtx_MC = (TF1*)ifile->Get(\"tf1_{itype}_{ax}_vs_nVtx_MC\")")
         
   def applyCorr(rdf, dftype):
      for itype in recoils:
         for ax in ["x", "y"]:
            rdf = rdf.Define(f"{itype}_{ax}_XYCorr", f"applyXYCorr({itype}_{ax}, {nPV}, tf1_{itype}_{ax}_vs_nVtx_{dftype})")
         rdf = rdf.Define(f"{itype}_pt_XYCorr", f"TMath::Sqrt({itype}_x_XYCorr*{itype}_x_XYCorr + {itype}_y_XYCorr*{itype}_y_XYCorr)")
         rdf = rdf.Define(f"{itype}_phi_XYCorr", f"TVector2::Phi_mpi_pi(TMath::ATan2({itype}_y_XYCorr, {itype}_x_XYCorr))")
      return rdf
            
   global rdf_data
   global rdf_MC
   rdf_data = applyCorr(rdf_data, "data")
   rdf_MC = applyCorr(rdf_MC, "MC")
   
   hprofs_data = OrderedDict()
   hprofs_MC = OrderedDict()
   
   hphis_data = OrderedDict()
   hphis_MC = OrderedDict()

   for ax in ["x", "y"]:
      hprofs_data[ax] = OrderedDict()
      hprofs_MC[ax] = OrderedDict()
   
   for itype in recoils:
      for ax in ["x", "y"]:
         hname = f"hprof_{itype}_{ax}_vs_nVtx_{suffix}"
         hname_data = f"{hname}_data"
         hname_MC = f"{hname}_MC"
         hprofs_data[ax][itype] = rdf_data.Profile1D((hname_data, hname_data, 70, 0, 70), nPV, itype + "_" + ax + "_XYCorr")
         hprofs_MC[ax][itype] = rdf_MC.Profile1D((hname_MC, hname_MC, 70, 0, 70), nPV, itype + "_" + ax + "_XYCorr")
         
      hphi_name = f"hphi_{itype}_{suffix}"
      hphi_name_data = f"{hphi_name}_data"
      hphi_name_MC = f"{hphi_name}_MC"
      hphis_data[itype] = rdf_data.Histo1D((hphi_name_data, f"{itype}_phi", 100, -3.2, 3.2), itype + "_phi_XYCorr")
      hphis_MC[itype] = rdf_MC.Histo1D((hphi_name_MC, f"{itype}_phi", 100, -3.2, 3.2), itype + "_phi_XYCorr")

   # fit the profiles with a pol1
   for itype in recoils:
      for ax in ["x", "y"]:
         hname = f"hprof_{itype}_{ax}_vs_nVtx_{suffix}"
         DrawHistos([hprofs_data[ax][itype], hprofs_MC[ax][itype]], [labels[itype] + " Data", labels[itype] + " MC"], 0, 70, "nPV", -10.0, 10.0, f"<MET_{ax}> [GeV]", hname, drawashist=True, dology=False, legendPos=[0.20, 0.70, 0.50, 0.90], mycolors=[colors[itype]]*2, linestyles = [linestyles[itype], linestyles[itype + "_MC"]], noLumi=noLumi, outdir=outdir)

      hphi_name = f"hphi_{itype}_{suffix}"
      DrawHistos([hphis_data[itype], hphis_MC[itype]], [labels[itype] + " Data", labels[itype] + " MC"], -3.2, 3.2, "#phi", 0, 0.019, "a.u.", hphi_name, drawashist=True, dology=False, legendPos=[0.20, 0.70, 0.50, 0.90], mycolors=[colors[itype], colors[itype + "_MC"]], linestyles = [linestyles[itype], linestyles[itype + "_MC"]], noLumi=noLumi, outdir=outdir, donormalize=True)
      
   h_toDraws = OrderedDict()
   for itype in recoils:
      h_toDraws[itype] = hphis_data[itype]
   for itype in recoils:
      h_toDraws[itype + "_MC"] = hphis_MC[itype]
   
   hcolors = [colors[itype] for itype in recoils]*2
   hmarkers = GetMarkers(h_toDraws)
   hlinestyles = GetLineStyles(h_toDraws)
   hdrawoptions = GetDrawOptions(h_toDraws)
   
   
   args = {
      "mycolors": hcolors,
      "markerstyles": hmarkers,
      "linestyles": hlinestyles,
      "drawoptions": hdrawoptions,
      "outdir": outdir,
      "noLumi": noLumi,
      "dology": False,
      "drawashist": False,
      "donormalize": True,
      "legendPos": [0.20, 0.70, 0.50, 0.90],
      "lheader": "After XY correction",
   }
   
   DrawHistos(h_toDraws.values(), [labels[itype] for itype in recoils], -3.2, 3.2, "#phi", 0, 0.019, "a.u.", f"hphi_{suffix}", **args)
      
   #DrawHistos([hphis_data[itype] for itype in recoils] + [hphis_MC[itype] for itype in recoils], [labels[itype] for itype in recoils], -3.2, 3.2, "#phi", 0, 0.019, "a.u.", f"hphi_{suffix}", drawashist=True, dology=False, legendPos=[0.20, 0.70, 0.50, 0.90], mycolors=[colors[itype] for itype in recoils]*2, linestyles = [linestyles[itype] for itype in recoils] + [linestyles[itype + "_MC"] for itpe in recoils], noLumi=noLumi, outdir=outdir, donormalize=True)
      
   if writeOutput:
      def regularizeNames(rdf):
         rdf = rdf.Define("MET_XYCorr_pt", "PF_pt_XYCorr") \
                  .Define("MET_XYCorr_phi", "PF_phi_XYCorr") \
                  .Define("PuppiMET_XYCorr_pt", "PUPPI_pt_XYCorr") \
                  .Define("PuppiMET_XYCorr_phi", "PUPPI_phi_XYCorr") \
                  .Define("DeepMETResolutionTune_XYCorr_pt", "DeepMET_pt_XYCorr") \
                  .Define("DeepMETResolutionTune_XYCorr_phi", "DeepMET_phi_XYCorr")
         return rdf
   
      rdf_data = regularizeNames(rdf_data)
      rdf_MC = regularizeNames(rdf_MC)
      
      global branches_data
      global branches_MC
      branches_data += ["MET_XYCorr_pt", "MET_XYCorr_phi", "PuppiMET_XYCorr_pt", "PuppiMET_XYCorr_phi", "DeepMETResolutionTune_XYCorr_pt", "DeepMETResolutionTune_XYCorr_phi"]
      branches_MC += ["MET_XYCorr_pt", "MET_XYCorr_phi", "PuppiMET_XYCorr_pt", "PuppiMET_XYCorr_phi", "DeepMETResolutionTune_XYCorr_pt", "DeepMETResolutionTune_XYCorr_phi"]
    
      odir = "/afs/cern.ch/work/y/yofeng/public/outputroot_withXYCorr"
      if not os.path.exists(odir):
         print(f"Creating {odir}")
         os.makedirs(odir)
         
      rdf_data.Snapshot("Events", f"{odir}/Data.root", branches_data)
      rdf_MC.Snapshot("Events", f"{odir}/DY.root", branches_MC)
         
   
   return
   
   
if __name__ == "__main__":
   print("Running the XY corrections")
   corrname = "results/MET_xy_dataMC.root"
   
   if 1:
      print(f"Deriving XY corrections for MET and saving them in {corrname}")
      deriveXYCorr(corrname, writeOutput=writeOutput)
      print("\n\n\n")
      
   if 1:
      print(f"Applying XY corrections to MET and saving the corrected output in {corrname}")
      applyXYCorr(corrname, writeOutput=writeOutput)
      print("\n\n\n")
      
   print("Done")