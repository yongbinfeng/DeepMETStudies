#ifndef FUNCTIONS
#define FUNCTIONS

/*
 * A few useful functions for the mT fit studies
 *
 */
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "ROOT/RVec.hxx"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

// reweight w pt
float wpt_slope_weight(float wpt, float offset, float slope, float cut = 20.0)
{
  if (wpt > cut)
    return (offset + slope * cut);
  // if(wpt > cut) return 1.0;
  float weight = offset + slope * wpt;
  return weight;
}

float rapidity(float pt, float eta, float mass)
{
  // calculate rapidity based on pt, eta and mass
  // https://en.wikipedia.org/wiki/Pseudorapidity
  float p = pt * TMath::CosH(eta);
  float e = TMath::Sqrt(p * p + mass * mass);
  float pz = pt * TMath::SinH(eta);
  return 0.5 * TMath::Log((e + pz) / (e - pz));
}

float wpt_alphaF_weight(float wpt, bool apply)
{

  if (!apply)
    return 1.0;

  float weights[] = {1.0135344827586206, 1.0123275862068966, 1.008706896551724, 1.005948275862069, 1.0037068965517242, 1.0026724137931033, 1.0033620689655172, 1.002844827586207, 1.0042241379310344, 1.0052586206896552, 1.005948275862069, 1.0075, 1.0083620689655173, 1.0085344827586207, 1.0095689655172413, 1.0085344827586207, 1.0083620689655173, 1.009741379310345, 1.005603448275862, 1.0037068965517242, 1.0038793103448276, 1.0033620689655172, 1.0057758620689654, 1.0050862068965516, 1.0075, 1.0064655172413792, 1.0080172413793103, 1.0078448275862069, 1.0078448275862069, 1.008706896551724, 1.0078448275862069, 1.008706896551724, 1.0088793103448275, 1.0081896551724139, 1.008706896551724, 1.009396551724138, 1.0092241379310345, 1.0088793103448275, 1.0106034482758621, 1.0106034482758621, 1.009396551724138, 1.0083620689655173};

  // one GeV per bin
  int bin = (int)wpt - 1;
  // overflow, use the max instead
  int nbins = sizeof(weights) / sizeof(weights[0]);
  if (bin >= nbins)
    bin = nbins - 1;
  return weights[bin];
}

float wpt_alphaF_weight_Linear(float wpt, float apply)
{
  if (fabs(apply) < 0.5)
    return 1.0;
  if (wpt < 8.0)
    return -0.0108 / 8.0 * wpt + 1.013535;
  if (wpt >= 8.0 && wpt < 25.0)
    return 0.0058 / 17.0 * (wpt - 8.0) + 1.0027;
  else
  {
    return 1.0085;
  }
}

float wpt_alphaF_weight_Pol4(float wpt, float apply)
{
  if (fabs(apply) < 0.5)
    return 1.0;
  double p[] = {1.01527, -0.0028649, 0.000226093, -6.44028e-06, 6.18964e-08};
  if (wpt > 40.0)
    wpt = 40.0;
  return p[0] + p[1] * wpt + p[2] * wpt * wpt + p[3] * wpt * wpt * wpt + p[4] * wpt * wpt * wpt * wpt;
}

float wpt_mc_weight(float wpt, float apply)
{
  if (fabs(apply) < 0.5)
    return 1.0;
  double p[] = {1.00613, -0.00227578, 0.00026437, -1.05851e-05};
  if (wpt > 10.0)
    wpt = 10.0;
  return p[0] + p[1] * wpt + p[2] * wpt * wpt + p[3] * wpt * wpt * wpt;
}

float wpt_LOPDF_weight(float wpt, float apply)
{
  if (fabs(apply) < 0.5)
    return 1.0;
  // double p[] = { 0.978286, 0.00384662, -0.000227073, 6.04121e-06,  -5.94203e-08 };
  double p[] = {0.978167, 0.00384693, -0.000225248, 5.93053e-06, -5.77007e-08};
  if (wpt > 40.0)
    wpt = 40.0;
  return p[0] + p[1] * wpt + p[2] * wpt * wpt + p[3] * wpt * wpt * wpt;
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

float calM_fromMET_PtEtaPhi(float lep_pt, float lep_eta, float lep_phi, float met_pt, float met_eta, float met_phi)
{
  // m2 = 2*pt1*pt2*(cosh(eta1-eta2)-cos(phi1-phi2))
  // https://en.wikipedia.org/wiki/Invariant_mass
  return TMath::Sqrt(2 * lep_pt * met_pt * (TMath::CosH(lep_eta - met_eta) - TMath::Cos(lep_phi - met_phi)));
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

Eigen::TensorFixedSize<int, Eigen::Sizes<2>> prefsrLeptons(const ROOT::VecOps::RVec<int> &status,
                                                           const ROOT::VecOps::RVec<int> &statusFlags, const ROOT::VecOps::RVec<int> &pdgId, const ROOT::VecOps::RVec<int> &motherIdx)
{

  const std::size_t ngenparts = status.size();

  std::array<std::size_t, 2> photos_idxs;
  std::array<std::size_t, 2> all_idxs;
  std::size_t nphotos = 0;
  std::size_t nall = 0;
  for (std::size_t i = 0; i < ngenparts; ++i)
  {
    const int &istatus = status[i];
    const int &istatusFlags = statusFlags[i];
    const int &ipdgId = pdgId[i];
    const int &imotherIdx = motherIdx[i];

    const int absPdgId = std::abs(ipdgId);

    const bool is_lepton = absPdgId >= 11 && absPdgId <= 16;
    const bool is_status746 = istatus == 746;
    const bool is_status23 = istatus == 23;
    const int &motherPdgId = pdgId[imotherIdx];
    const bool is_motherV = motherPdgId == 23 || std::abs(motherPdgId) == 24;

    const bool is_photos = is_lepton && is_status746 && is_motherV;

    const bool is_fromHardProcess = istatusFlags & (1 << 8);

    // TODO: Is there a way to relax the fromHardProcess condition?
    const bool is_other = is_lepton && (is_motherV || is_status23) && is_fromHardProcess;

    // If there are status = 746 leptons, they came from photos and are pre-FSR
    // (but still need to check the mother in case photos was applied to other particles in the
    // event, and in case of radiation from taus and their decay products)
    const bool is_all = is_photos || is_other;

    if (is_photos)
    {
      if (nphotos < 2)
      {
        photos_idxs[nphotos] = i;
      }
      ++nphotos;
    }

    if (is_all)
    {
      if (nall < 2)
      {
        all_idxs[nall] = i;
      }
      ++nall;
    }
  }

  std::array<std::size_t, 2> selected_idxs;
  if (nphotos == 2)
  {
    selected_idxs = photos_idxs;
  }
  else if (nall == 2)
  {
    selected_idxs = all_idxs;
  }
  else
  {
    throw std::range_error("Expected to find 2 pre-FSR leptons, but found " + std::to_string(nall) + ", " + std::to_string(nphotos));
  }

  std::array<int, 2> selected_pdgids = {pdgId[selected_idxs[0]], pdgId[selected_idxs[1]]};
  const bool partIdx = selected_pdgids[0] > 0;
  Eigen::TensorFixedSize<int, Eigen::Sizes<2>> out;
  out(0) = selected_idxs[!partIdx];
  out(1) = selected_idxs[partIdx];
  return out;
}

#endif
