'''
some useful functions to sync the pt and nvtx bins,
and prepare some variables (uparal, uperp, responses, and paral_diff) 
for the RDataFrame
'''


def doPAS():
    return False


def getqTRange():
    return (0, 350.0)


def getVtxRange():
    return (4.5, 40.5)


def getNPVString():
    return "PV_npvs"


def getqTLabel():
    return "q_{T} [GeV]"


def getnVtxLabel():
    return "N_{vtx}"


def getUparalLabel():
    return "#sigma (u_{#parallel}) [GeV]"


def getUperpLabel():
    return "#sigma (u_{#perp } ) [GeV]"


def getResponseLabel():
    return "-<u_{#parallel}>/<q_{T}>"


def getpTBins():
    import numpy as np
    xbins_qT = np.zeros(19)
    for ibin in range(0, 4):
        xbins_qT[ibin] = 5*ibin
    for ibin in range(4, 8):
        xbins_qT[ibin] = 10*ibin-20
    for ibin in range(8, 13):
        xbins_qT[ibin] = 20*ibin-100
    for ibin in range(13, 18):
        xbins_qT[ibin] = 30*ibin-220
    xbins_qT[18] = 350.0
    # xbins_qT[18:] = (350, 400, 450, 500)
    return xbins_qT


def getpTResponseBins():
    import numpy as np
    xbins_qT_resp = np.zeros(2)
    xbins_qT_resp[0] = 150.0
    xbins_qT_resp[1] = 1000.0
    return xbins_qT_resp


def getnVtxBins():
    import numpy as np
    # xbins_nVtx = np.array([0.,5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21.5, 23.5, 25.5, 27.5, 29.5, 31.5, 33.5, 38.5, 43.5, 50.5])
    # xbins_nVtx = np.array([0., 4.5, 7.5, 9.5, 11.5, 13.5, 15.5,
    #                      17.5, 19.5, 21.5, 23.5, 27.5, 31.5, 35.5, 40.5, 50.5, 60.5])
    xbins_nVtx = np.array([4.5, 7.5, 9.5, 11.5, 13.5,
                          15.5, 17.5, 19.5, 21.5, 23.5, 27.5, 31.5, 35.5, 40.6])
    return xbins_nVtx


def prepVars(rdf, u, utruth):
    """
    prepare the u_paral, perp, response and paral_diff variables
    """
    f_uparal = "( {U}_x * {UTRUTH}_x + {U}_y * {UTRUTH}_y ) / {UTRUTH}_pt".format(U=u, UTRUTH=utruth)
    f_uperp = "(-{U}_x * {UTRUTH}_y + {U}_y * {UTRUTH}_x ) / {UTRUTH}_pt".format(U=u, UTRUTH=utruth)
    f_ureponse = "{U}_paral / {UTRUTH}_pt".format(U=u, UTRUTH=utruth)
    f_uparal_diff = "{U}_paral - {UTRUTH}_pt".format(U=u, UTRUTH=utruth)

    rdf = rdf.Define("{U}_paral"     .format(U=u),     f_uparal) \
             .Define("{U}_perp"      .format(U=u),     f_uperp) \
             .Define("{U}_response"  .format(U=u),     f_ureponse) \
             .Define("{U}_paral_diff".format(U=u),     f_uparal_diff)
    return rdf


def prepRecoilVars(rdf, v, metnames):
    """
    prepare the recoil x, y, and pt variables given V and MET pt
    """
    if isinstance(metnames, str):
        metnames = [metnames]

    for metname in metnames:
        f_ux = "-({V}_pt * TMath::Cos({V}_phi) + {MET}_pt * TMath::Cos({MET}_phi))".format(V=v, MET=metname)
        f_uy = "-({V}_pt * TMath::Sin({V}_phi) + {MET}_pt * TMath::Sin({MET}_phi))".format(V=v, MET=metname)
        f_ut = "TMath::Sqrt(u_{MET}_x * u_{MET}_x + u_{MET}_y * u_{MET}_y)".format(MET=metname)

        print(f"MET name {metname}, f_ux {f_ux}, f_uy {f_uy}, f_ut {f_ut}")

        rdf = rdf.Define("u_{MET}_x".format(MET=metname), f_ux) \
                 .Define("u_{MET}_y".format(MET=metname), f_uy) \
                 .Define("u_{MET}_pt".format(MET=metname), f_ut)
    return rdf


get_response_code = '''
float get_response(float qt, TProfile* h_response)
{
    int maxbin = h_response->GetNbinsX();
    int ibin = h_response->FindBin( qt );
    if( ibin>maxbin ){
        ibin = maxbin;
    }
    return h_response->GetBinContent( ibin );
}

float get_response(float qt, TH1* h_response)
{
    int maxbin = h_response->GetNbinsX();
    int ibin = h_response->FindBin( qt );
    if( ibin>maxbin ){
        ibin = maxbin;
    }
    return h_response->GetBinContent( ibin );
}
'''
