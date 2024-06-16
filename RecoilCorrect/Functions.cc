#ifndef FUNCTIONS
#define FUNCTIONS

#include "TLorentzVector.h"
#include "TVector2.h"
#include "TF1.h"
#include "TMath.h"
#include "TH1F.h"
#include "TList.h"
#include <iostream>

TVector2 SetMagPhi(float pt, float phi)
{
    // in the default SetMagPhi in TVector2 in root,
    // mag is set to positive without changing the direction
    // https://github.com/root-project/root/blob/master/math/physics/src/TVector2.cxx#L124
    if (pt < 0)
        phi = phi + TMath::Pi();
    TVector2 v2;
    v2.SetMagPhi(pt, phi);
    return v2;
}

TLorentzVector VVec(float lep0_pt, float lep0_eta, float lep0_phi, float lep0_e,
                    float lep1_pt, float lep1_eta, float lep1_phi, float lep1_e)
{
    TLorentzVector vlep0;
    vlep0.SetPtEtaPhiE(lep0_pt, lep0_eta, lep0_phi, lep0_e);
    TLorentzVector vlep1;
    vlep1.SetPtEtaPhiE(lep1_pt, lep1_eta, lep1_phi, lep1_e);

    return (vlep0 + vlep1);
}

TLorentzVector VVecM(float lep0_pt, float lep0_eta, float lep0_phi, float lep0_m,
                     float lep1_pt, float lep1_eta, float lep1_phi, float lep1_m)
{
    TLorentzVector vlep0;
    vlep0.SetPtEtaPhiM(lep0_pt, lep0_eta, lep0_phi, lep0_m);
    TLorentzVector vlep1;
    vlep1.SetPtEtaPhiM(lep1_pt, lep1_eta, lep1_phi, lep1_m);

    return (vlep0 + vlep1);
}

TVector2 UVec(float leps_pt, float leps_phi,
              float met_pt, float met_phi)
{
    TVector2 vleps = SetMagPhi(leps_pt, leps_phi);
    TVector2 vmet = SetMagPhi(met_pt, met_phi);

    return (vleps + vmet) * (-1.0);
}

TVector2 METVec(float leps_pt, float leps_phi,
                float uparal, float uperp)
{
    TVector2 vleps = SetMagPhi(leps_pt, leps_phi);
    TVector2 vparal = SetMagPhi(uparal, leps_phi + TMath::Pi());
    TVector2 vperp = SetMagPhi(uperp, leps_phi + TMath::Pi() / 2.0 * 3.0);

    return (vleps + vparal + vperp) * (-1.0);
}

TVector2 METVec(float lep_pt, float lep_phi,
                float uparal, float uperp,
                float V_phi)
{
    TVector2 vlep = SetMagPhi(lep_pt, lep_phi);
    TVector2 vparal = SetMagPhi(uparal, V_phi + TMath::Pi());
    TVector2 vperp = SetMagPhi(uperp, V_phi + TMath::Pi() / 2.0 * 3.0);
    return (vlep + vparal + vperp) * (-1.0);
}

Float_t UCorrection(float u, float zpt, TF1 *f_mc_mean, TF1 *f_mc_sigma, TF1 *f_data_mean, TF1 *f_data_sigma)
{
    float zptmax = f_mc_mean->GetXmax();
    if (zpt > zptmax)
        return u;
    float u_mc_mean = f_mc_mean->Eval(zpt);
    float u_mc_sigma = f_mc_sigma->Eval(zpt);
    float u_data_mean = f_data_mean->Eval(zpt);
    float u_data_sigma = f_data_sigma->Eval(zpt);
    float u_corr = (u - u_mc_mean) / (u_mc_sigma + 0.0001) * u_data_sigma + u_data_mean;
    return u_corr;
}

Float_t UCorrection(float u, float zpt, TF1 *f_mc_sigma, TF1 *f_data_sigma)
{
    float zptmax = f_mc_sigma->GetXmax();
    if (zpt > zptmax)
        return u;
    float u_mc_sigma = f_mc_sigma->Eval(zpt);
    float u_data_sigma = f_data_sigma->Eval(zpt);
    float u_corr = u / (u_mc_sigma + 0.0001) * u_data_sigma;
    return u_corr;
}

Float_t UCorrection_Quant(float u, float zpt, TH1F *hptbins, TList *tfs_Data, TList *tfs_MC, float qcut = 0.005, bool isU1Diff = false)
{
    float zptmax = hptbins->GetXaxis()->GetXmax();
    float zptmin = hptbins->GetXaxis()->GetXmin();
    // outlier in the z pt bin
    if (zpt > zptmax || zpt < zptmin)
        return u;
    int index = hptbins->FindBin(zpt) - 1;
    TF1 *tf1_MC = (TF1 *)tfs_MC->At(index);
    TF1 *tf1_Data = (TF1 *)tfs_Data->At(index);

    if (u > tf1_MC->GetXmax() || u < tf1_MC->GetXmin())
        return u;

    float u_temp = u;
    if (isU1Diff)
    {
        u_temp -= zpt;
    }
    float qt = tf1_MC->Eval(u_temp); // quantile, cdf
    // no correction on outliers
    if (qcut > 0)
    {
        if (qt < qcut || qt > 1 - qcut)
            return u;
    }
    float u_corr = tf1_Data->GetX(qt);
    if (isU1Diff)
    {
        u_corr += zpt;
    }
    return u_corr;
}

Float_t UCorrection_Quant(float u, int njet, float zpt, TH1F *hjetbins, TH1F *hptbins, TList *tfs_Data, TList *tfs_MC, float qcut = 0.005, bool isU1Diff = false)
{
    float zptmax = hptbins->GetXaxis()->GetXmax();
    float njetmax = hjetbins->GetXaxis()->GetXmax();
    float zptmin = hptbins->GetXaxis()->GetXmin();
    float njetmin = hjetbins->GetXaxis()->GetXmin();
    // outlier
    if (zpt > zptmax || njet > njetmax || zpt < zptmin || njet < njetmin)
        return u;
    int index = hptbins->FindBin(zpt) - 1;
    int njetbin = hjetbins->FindBin(njet) - 1;
    TF1 *tf1_MC = (TF1 *)((TList *)tfs_MC->At(njetbin))->At(index);
    TF1 *tf1_Data = (TF1 *)((TList *)tfs_Data->At(njetbin))->At(index);

    if (u > tf1_MC->GetXmax() || u < tf1_MC->GetXmin())
        return u;

    float u_temp = u;
    if (isU1Diff)
    {
        u_temp -= zpt;
    }
    float qt = tf1_MC->Eval(u_temp); // quantile, cdf
    // no correction on outliers
    if (qcut > 0)
    {
        if (qt < qcut || qt > 1 - qcut)
            return u;
    }
    float u_corr = tf1_Data->GetX(qt);
    if (isU1Diff)
    {
        u_corr += zpt;
    }
    if (isnan(u_corr))
    {
        std::cout << "u_corr is nan" << std::endl;
        std::cout << "u: " << u << " zpt: " << zpt << " njet: " << njet << std::endl;
        std::cout << "qt: " << qt << std::endl;
        std::cout << "tf1_MC: " << tf1_MC->GetName() << std::endl;
        std::cout << "tf1_Data: " << tf1_Data->GetName() << std::endl;
        // return original value
        return u;
    }
    return u_corr;
}

Float_t ZptReWeight(float zpt, TH1D *h_zpt_ratio_data_vs_MC, bool isData = false)
{
    if (isData)
        return 1.0;
    float zptmax = h_zpt_ratio_data_vs_MC->GetXaxis()->GetXmax();
    if (zpt > zptmax)
        return 1.0;
    return h_zpt_ratio_data_vs_MC->GetBinContent(h_zpt_ratio_data_vs_MC->FindBin(zpt));
}

// calculate mT based on lepton px, py and recoil ux, uy
float calMT_fromUT_XY(float lep_px, float lep_py, float ux, float uy)
{
    float lep_pt = TMath::Sqrt(lep_px * lep_px + lep_py * lep_py);
    float u_pt = TMath::Sqrt(ux * ux + uy * uy);
    float sum_pt = TMath::Sqrt((lep_px + ux) * (lep_px + ux) + (lep_py + uy) * (lep_py + uy));
    float mT2 = 2 * (lep_pt * (lep_pt + sum_pt) + lep_px * ux + lep_py * uy);
    return TMath::Sqrt(mT2);
}

float calMT_fromUT_PtPhi(float lep_pt, float lep_phi, float ut, float uphi)
{
    return calMT_fromUT_XY(lep_pt * TMath::Cos(lep_phi), lep_pt * TMath::Sin(lep_phi), ut * TMath::Cos(uphi), ut * TMath::Sin(uphi));
}

float calMT_fromMET_XY(float lep_px, float lep_py, float met_x, float met_y)
{
    float ux = -(lep_px + met_x);
    float uy = -(lep_py + met_y);
    return calMT_fromUT_XY(lep_px, lep_py, ux, uy);
}

float calMT_fromMET_PtPhi(float lep_pt, float lep_phi, float met_pt, float met_phi)
{
    return calMT_fromMET_XY(lep_pt * TMath::Cos(lep_phi), lep_pt * TMath::Sin(lep_phi), met_pt * TMath::Cos(met_phi), met_pt * TMath::Sin(met_phi));
}

float calUT_XY(float lep_px, float lep_py, float met_x, float met_y)
{
    float ux = -(lep_px + met_x);
    float uy = -(lep_py + met_y);
    return TMath::Sqrt(ux * ux + uy * uy);
}

float calUT_PtPhi(float lep_pt, float lep_phi, float met_pt, float met_phi)
{
    return calUT_XY(lep_pt * TMath::Cos(lep_phi), lep_pt * TMath::Sin(lep_phi), met_pt * TMath::Cos(met_phi), met_pt * TMath::Sin(met_phi));
}

Int_t trueWIndex(int trueWn)
{
    if (trueWn == 1)
    {
        return 0;
    }
    else if (trueWn == 2)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

float applyXYCorr(float met_x, float npv, TF1 *f1)
{
    return met_x - f1->Eval(npv);
}

#endif
