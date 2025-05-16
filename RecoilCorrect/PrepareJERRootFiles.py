import ROOT
import numpy as np
import sys

ROOT.ROOT.EnableImplicitMT()

def convert_jer_sf_to_ROOT():
    table = []
    f_name = "corrections/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt"
    with open(f_name) as f:
        for line in f:
            if line.startswith("#") or line.startswith("{") or not line.strip():
                continue
            toks = line.strip().split()
            eta_min, eta_max = float(toks[0]), float(toks[1])
            sf_nom = float(toks[3])
            sf_down = float(toks[4])
            sf_up = float(toks[5])
            table.append((eta_min, eta_max, sf_nom, sf_down, sf_up))
            
    # make a th1d
    binnings = np.array([t[0] for t in table] + [table[-1][1]])
    print("binnings:", binnings)
    h_sfs_norm = ROOT.TH1D("h_sfs_norm", "JER SF", len(binnings)-1, binnings)
    h_sfs_down = ROOT.TH1D("h_sfs_down", "JER SF down", len(binnings)-1, binnings)
    h_sfs_up = ROOT.TH1D("h_sfs_up", "JER SF up", len(binnings)-1, binnings)
    for i, t in enumerate(table):
        h_sfs_norm.SetBinContent(i+1, t[2])
        h_sfs_down.SetBinContent(i+1, t[3])
        h_sfs_up.SetBinContent(i+1, t[4])
        
    return h_sfs_norm, h_sfs_down, h_sfs_up

def convert_jer_resol_to_ROOT():
    table = []
    f_name = "corrections/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt"
    with open(f_name) as f:
        for line in f:
            if line.startswith("#") or line.startswith("{") or not line.strip():
                continue
            toks = line.strip().split()
            eta_min, eta_max = float(toks[0]), float(toks[1])
            rho_min, rho_max = float(toks[2]), float(toks[3])
            par0, par1, par2, par3 = float(toks[7]), float(toks[8]), float(toks[9]), float(toks[10])
            table.append((eta_min, eta_max, rho_min, rho_max, par0, par1, par2, par3))
            
    # make a th2d
    n_rho_bins = 7
    binnings_eta = []
    for i in range(0, len(table), n_rho_bins):
        binnings_eta.append(table[i][0])
    binnings_eta.append(table[-1][1])
    print("binnings_eta:", binnings_eta)
    
    binnings_rho = []
    for i in range(n_rho_bins):
        binnings_rho.append(table[i][2])
    binnings_rho.append(table[-1][3])
    print("binnings_rho:", binnings_rho)
    
    h_resol_par0 = ROOT.TH2D("h_resol_par0", "JER Resolution par0", len(binnings_eta)-1, np.array(binnings_eta), len(binnings_rho)-1, np.array(binnings_rho))
    h_resol_par1 = ROOT.TH2D("h_resol_par1", "JER Resolution par1", len(binnings_eta)-1, np.array(binnings_eta), len(binnings_rho)-1, np.array(binnings_rho))
    h_resol_par2 = ROOT.TH2D("h_resol_par2", "JER Resolution par2", len(binnings_eta)-1, np.array(binnings_eta), len(binnings_rho)-1, np.array(binnings_rho))
    h_resol_par3 = ROOT.TH2D("h_resol_par3", "JER Resolution par3", len(binnings_eta)-1, np.array(binnings_eta), len(binnings_rho)-1, np.array(binnings_rho))
    
    for ieta in range(1, len(binnings_eta)):
        for irho in range(1, len(binnings_rho)):
            t = table[(ieta-1)*n_rho_bins + irho - 1]
            assert t[0] == binnings_eta[ieta-1] and t[1] == binnings_eta[ieta], f"Error in eta binning, ieta {ieta} table: {t}, binnings: {binnings_eta[ieta-1]}, {binnings_eta[ieta]}"
            assert t[2] == binnings_rho[irho-1] and t[3] == binnings_rho[irho], f"Error in rho binning, irho {irho} table: {t}, binnings: {binnings_rho[irho-1]}, {binnings_rho[irho]}"
            h_resol_par0.SetBinContent(ieta, irho, t[4])
            h_resol_par1.SetBinContent(ieta, irho, t[5])
            h_resol_par2.SetBinContent(ieta, irho, t[6])
            h_resol_par3.SetBinContent(ieta, irho, t[7])
    
    return h_resol_par0, h_resol_par1, h_resol_par2, h_resol_par3

def saveJERInfo_to_ROOT():
    # Save the JER SF and resolution information to ROOT files
    th1d_sfs_norm, th1d_sfs_down, th1d_sfs_up = convert_jer_sf_to_ROOT()
    th2d_resol_par0, th2d_resol_par1, th2d_resol_par2, th2d_resol_par3 = convert_jer_resol_to_ROOT()
    
    # TF1s
    tf1s = []
    for i in range(26):
        for j in range(7):
            tf1 = ROOT.TF1(f"h_resol_pol3_{i}_{j}", "sqrt([0]*abs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])", 15.0, 3000)
            tf1.SetParameter(0, th2d_resol_par0.GetBinContent(i+1, j+1))
            tf1.SetParameter(1, th2d_resol_par1.GetBinContent(i+1, j+1))
            tf1.SetParameter(2, th2d_resol_par2.GetBinContent(i+1, j+1))
            tf1.SetParameter(3, th2d_resol_par3.GetBinContent(i+1, j+1))
            
            tf1s.append(tf1)

    output_file = ROOT.TFile("corrections/jer.root", "RECREATE")
    th1d_sfs_norm.Write()
    th1d_sfs_down.Write()
    th1d_sfs_up.Write()
    th2d_resol_par0.Write()
    th2d_resol_par1.Write()
    th2d_resol_par2.Write()
    th2d_resol_par3.Write()
    for tf1 in tf1s:
        tf1.Write()
    output_file.Close()
    
    
def ApplyJER():
    load_code = """
    // Load JER SF and resolution information from ROOT files
    TFile* f = TFile::Open("corrections/jer.root");
    TH1D* h_sfs_norm = (TH1D*)f->Get("h_sfs_norm");
    TH1D* h_sfs_down = (TH1D*)f->Get("h_sfs_down");
    TH1D* h_sfs_up = (TH1D*)f->Get("h_sfs_up");
    TH2D* h_resol_par0 = (TH2D*)f->Get("h_resol_par0");

    TF1* h_resols[26][7];
    for (int i = 0; i < 26; ++i) {{
        for (int j = 0; j < 7; ++j) {{
            h_resols[i][j] = (TF1*)f->Get(Form("h_resol_pol3_%d_%d", i, j));
        }}
    }}
    """
    ROOT.gInterpreter.ProcessLine(load_code)

    # Inject tables as C++ constants
    ROOT.gInterpreter.Declare(f"""
    #include <cmath>
    #include <vector>
    #include <tuple>
    #include "TVector2.h"
    using namespace ROOT::VecOps;


    float getJERSF(float eta, int norm=1) {{
        int ibin = h_sfs_norm->FindBin(eta);
        if (ibin < 1 || ibin > h_sfs_norm->GetNbinsX()) {{
            return 1.0;
        }} 
        if (norm == 1) {{
            return h_sfs_norm->GetBinContent(ibin);
        }} else if (norm == 0) {{
            return h_sfs_down->GetBinContent(ibin);
        }} else if (norm == 2) {{
            return h_sfs_up->GetBinContent(ibin);
        }} else {{
            std::cerr << "Invalid norm value: " << norm << std::endl;
            return 1.0;
        }}
    }}

    float getJERResolution(float pt, float eta, float rho) {{
        int ibin_eta = h_resol_par0->GetXaxis()->FindBin(eta);
        int ibin_rho = h_resol_par0->GetYaxis()->FindBin(rho);
        if (ibin_eta < 1 || ibin_eta > h_resol_par0->GetNbinsX() || ibin_rho < 1 || ibin_rho > h_resol_par0->GetNbinsY()) {{
            return 0.0;
        }}
        return h_resols[ibin_eta-1][ibin_rho-1]->Eval(pt);
    }}

    RVec<float> smearJER(RVec<float> pt, RVec<float> eta, RVec<float> genpt, float rho) {{
        RVec<float> out;
        for (size_t i = 0; i < pt.size(); ++i) {{
            float sf = getJERSF(eta[i]);
            float res = getJERResolution(pt[i], eta[i], rho);
            float gpt = genpt[i];
            float smeared = pt[i];
            if (gpt > 0) {{
                smeared = gpt + sf * (pt[i] - gpt);
            }} else {{
                float sigma = res * sqrt(std::max(sf*sf - 1.f, 0.f));
                smeared = pt[i] * (1.0 + sigma * gRandom->Gaus(0, 1));
            }}
            out.push_back(smeared);
        }}
        return out;
    }}

    TVector2 recomputeMET(float met_pt, float met_phi, RVec<float> pt_old, RVec<float> pt_new, RVec<float> phi) {{
        float met_px = met_pt * cos(met_phi);
        float met_py = met_pt * sin(met_phi);
        for (size_t i = 0; i < pt_old.size(); ++i) {{
            float dpx = (pt_old[i] - pt_new[i]) * cos(phi[i]);
            float dpy = (pt_old[i] - pt_new[i]) * sin(phi[i]);
            met_px += dpx;
            met_py += dpy;
        }}
        TVector2 met_vec(met_px, met_py);
        return met_vec;
    }}
    """)

    # --- Load your NanoAOD file with RDataFrame ---
    fname = "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/DY.root"
    rdf = ROOT.RDataFrame("Events", fname)
    branches = rdf.GetColumnNames()
    # look at the first 10 events
    #rdf = rdf.Filter("rdfentry_ < 1000")
    rdf = rdf.Define("Jet_pt_smeared", "smearJER(Jet_pt_tight, Jet_eta_tight, Jet_truthpt_tight,fixedGridRhoFastjetAll)")

    rdf = rdf.Define("vMET_smeared", "recomputeMET(MET_pt, MET_phi, Jet_pt_tight, Jet_pt_smeared, Jet_phi_tight)") \
        .Define("MET_pt_smeared", "vMET_smeared.Mod()") \
        .Define("MET_phi_smeared", "vMET_smeared.Phi()")

    # save all the branches
    branches += ["Jet_pt_smeared", "MET_pt_smeared", "MET_phi_smeared"]
    rdf.Snapshot("Events", "/home/yongbinfeng/Desktop/DeepMET/data/outputroot_withjets/DY_smeared.root", branches)
    
if __name__ == "__main__":
    # Convert JER SF and resolution information to ROOT format
    saveJERInfo_to_ROOT()
    
    # Apply JER to the NanoAOD file
    ApplyJER()