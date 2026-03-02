
"""
psychro_aux_vec.py — Versiones vectoriales (NumPy-first) de VIRCOE, ENHAC, MSTAIR
- Aceptan escalar o ndarray y devuelven escalar/array según la entrada.
- Preparado para Numba en un paso posterior (kernels puros, sin objetos).

Requisitos: numpy
Dependen de: iapws_backend.pws_T (para pws) y, opcionalmente, tws_P si se usa.
"""
import numpy as np
from .iapws_backend import pws_T, water_kT_and_v, water_enthalpy_liq

# ---- Constantes compartidas ----
Rm = 8.314472   # kJ/kmol.K
Ma = 28.966
Mw = 18.015268
MR = Mw / Ma

# IF97 region 1 coef (para ISOTCOMP rama T>=273.15)
Ii = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,
               1.0,1.0,2.0,2.0,2.0,2.0,2.0,3.0,3.0,3.0,4.0,
               4.0,4.0,5.0,8.0,8.0,21.0,23.0,29.0,30.0,31.0,
               32.0])
Ji = np.array([-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,-9.0,-7.0,-1.0,
               0.0,1.0,3.0,-3.0,0.0,1.0,3.0,17.0,-4.0,0.0,6.0,
               -5.0,-2.0,10.0,-8.0,-11.0,-6.0,-29.0,-31.0,-38.0,
               -39.0,-40.0,-41.0])
ni = np.array([0.14632971213167,-0.84548187169114,-3.756360367204,
               3.3855169168385,-0.95791963387872,0.15772038513228,
               -0.016616417199501,0.00081214629983568,0.00028319080123804,
               -0.00060706301565874,-0.018990068218419,-0.032529748770505,
               -0.021841717175414,-0.00005283835796993,
               -0.00047184321073267,-0.00030001780793026,
               0.000047661393906987,-0.0000044141845330846,
               -0.00000000000000072694996297594,-0.000031679644845054,
               -0.0000028270797985312,-0.00000000085205128120103,
               -0.0000022425281908,-0.00000065171222895601,
               -0.00000000000014341729937924,-0.00000040516996860117,
               -0.0000000012734301741641,-0.00000000017424871230634,
               -6.8762131295531E-19, 1.4478307828521E-20,
               2.6335781662795E-23,-1.1947622640071E-23,
               1.8228094581404E-24,-9.3537087292458E-26])

# IAPWS-95 (sólido/líquido bajo 273.15 K)
g0 = np.array([-632020.233449497, 0.655022213658955,
               -0.0000000189369929326131, 0.00000000000000339746123271053,
               -5.56464869058991E-22])
s0 = -3327.33756492168
t1 = 0.0368017112855051 + 0.0510878114959572j
r1 = 44.7050716285388 + 65.6876847463481j
t2 = 0.337315741065416 + 0.335449415919319j
r20 = -72.597457432922 - 78.100842711287j
r21 = -0.0000557107698030123 + 0.0000464578634580806j
r22 = 0.0000000000234801409215913 - 0.0000000000285651142904972j

# ---- Helpers ----
def _as_array(x):
    return np.asarray(x, dtype=float)

def _scalarize_like(ref_in, out_arr):
    return float(out_arr[0]) if np.ndim(ref_in) == 0 else out_arr

def _broadcast_same(*args):
    arrs = [_as_array(a) for a in args]
    return np.broadcast_arrays(*arrs)

# ---- VIRCOE (vectorizado sobre T) ----
def VIRCOE_vec(Tdb):
    """Devuelve [Baa,Caaa,Bww,Cwww,Baw,Caaw,Caww,dBaadT,dCaaadT,dBwwdT,dCwwwdT,dBawdT,dCaawdT,dCawwdT]
       vectorizado sobre Tdb (°C)."""
    Tdb = _as_array(Tdb)
    TDBK = Tdb + 273.15

    # Constantes
    D = 1e-10
    Md = 0.0104477
    Tj = 132.6312
    Tc = 647.096
    Mc = 0.322
    Tred = 100.0
    Mcd = Mc / Mw

    # Coef Baa/Caaa
    Nk = np.array([0.118160747229, 0.713116392079, -1.61824192067,
                   0.0714140178971, -0.0865421396646, 0.134211176704,
                   0.0112626704218, -0.0420533228842, 0.0349008431982,
                   0.000164957183186, -0.101365037912, -0.17381369097,
                   -0.0472103183731, -0.0122523554253, -0.146629609713,
                   -0.0316055879821, 0.000233594806142, 0.0148287891978,
                   -0.00938782884667])
    ik = np.array([1.,1.,1.,2.,3.,3.,4.,4.,4.,6.,1.,3.,5.,6.,1.,3.,11.,1.,3.])
    jk = np.array([0.,0.33,1.01,0.,0.,0.15,0.,0.2,0.35,1.35,1.6,0.8,0.95,1.25,3.6,6.,3.25,3.5,15.])
    lk = np.array([0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,2.,2.,2.,3.,3.])

    Tau_j = Tj / TDBK[:, None]  # (N,1)

    idx = np.arange(Nk.size)
    m1 = idx < 10
    m2 = ~m1
    Ik = ik[None, :]
    Jk = jk[None, :]
    Lk = lk[None, :]
    Nk_ = Nk[None, :]

    Dpow_ikm1 = D**(Ik - 1.0)
    Tau_pow = Tau_j**Jk
    Dpow_lk = D**Lk
    exp_term = np.exp(-Dpow_lk)

    F1 = np.sum(Ik[:, m1] * Nk_[:, m1] * Dpow_ikm1[:, m1] * Tau_pow[:, m1], axis=1)
    Y1 = np.sum(Ik[:, m1] * Jk[:, m1] * Nk_[:, m1] * Dpow_ikm1[:, m1] * Tau_j**(Jk[:, m1]-1.0), axis=1)

    F2 = np.sum(Nk_[:, m2] * Dpow_ikm1[:, m2] * Tau_pow[:, m2] * exp_term[:, m2] * (Ik[:, m2] - Lk[:, m2]*Dpow_lk[:, m2]), axis=1)
    Y2 = np.sum(Jk[:, m2] * Nk_[:, m2] * Dpow_ikm1[:, m2] * Tau_j**(Jk[:, m2]-1.0) * exp_term[:, m2] * (Ik[:, m2] - Lk[:, m2]*Dpow_lk[:, m2]), axis=1)

    Baa = (F1 + F2) / Md
    dBaadT = -(Tau_j[:,0]**2) * (Y1 + Y2) / (Md * Tj)

    # Caaa
    Dpow_ikm2 = D**(Ik - 2.0)
    F3 = np.sum(Ik[:, m1]*(Ik[:, m1]-1.0) * Nk_[:, m1] * Dpow_ikm2[:, m1] * Tau_pow[:, m1], axis=1)
    Y3 = np.sum(Ik[:, m1]*(Ik[:, m1]-1.0) * Jk[:, m1] * Nk_[:, m1] * Dpow_ikm2[:, m1] * Tau_j**(Jk[:, m1]-1.0), axis=1)

    Gfac = ((Ik - Lk*Dpow_lk)*(Ik - 1.0 - Lk*Dpow_lk) - (Lk**2)*Dpow_lk)
    F4 = np.sum(Nk_[:, m2] * Dpow_ikm2[:, m2] * Tau_pow[:, m2] * exp_term[:, m2] * Gfac[:, m2], axis=1)
    Y4 = np.sum(Jk[:, m2] * Nk_[:, m2] * Dpow_ikm2[:, m2] * Tau_j**(Jk[:, m2]-1.0) * exp_term[:, m2] * Gfac[:, m2], axis=1)

    Caaa = (F3 + F4) / (Md**2)
    dCaaadT = - (Tau_j[:,0]**2) * (Y3 + Y4) / (Md**2 * Tj)

    # Bww/Cwww
    ci = np.array([0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
                   1.,1.,1.,1.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,2.,
                   2.,2.,2.,2.,2.,2.,3.,3.,3.,3.,4.,6.,6.,6.,6.,0.,0.,0.,
                   3.5,3.5])
    di = np.array([1.,1.,1.,2.,2.,3.,4.,1.,1.,1.,2.,2.,3.,4.,4.,5.,7.,9.,
                   10.,11.,13.,15.,1.,2.,2.,2.,3.,4.,4.,4.,5.,6.,6.,7.,9.,
                   9.,9.,9.,9.,10.,10.,12.,3.,4.,4.,5.,14.,3.,6.,6.,6.,3.,
                   3.,3.,0.85,0.95])
    ti = np.array([-0.5,0.875,1.,0.5,0.75,0.375,1.,4.,6.,12.,1.,5.,4.,2.,13.,9.,3.,4.,11.,4.,13.,1.,
                   7.,1.,9.,10.,10.,3.,7.,10.,10.,6.,10.,10.,1.,2.,3.,4.,8.,6.,9.,8.,16.,22.,23.,23.,
                   10.,50.,44.,46.,50.,0.,1.,4.,0.2,0.2])
    ni56 = np.array([0.012533547935523,7.8957634722828,-8.7803203303561,
                   0.31802509345418,-0.26145533859358,-0.0078199751687981,
                   0.0088089493102134,-0.66856572307965,0.20433810950965,
                   -6.6212605039687E-05,-0.19232721156002,-0.25709043003438,
                   0.16074868486251,-0.040092828925807,3.9343422603254E-07,
                   -7.5941377088144E-06,0.00056250979351888,-1.5608652257135E-05,
                   1.1537996422951E-09,3.6582165144204E-07,-1.3251180074668E-12,
                   -6.2639586912454E-10,-0.10793600908921,0.017611491008752,
                   0.22132295167546,-0.40247669763528,0.58083399985759,
                   0.0049969146990806,-0.031358700712549,-0.74315929710341,
                   0.4780732991548,0.020527940895948,-0.13636435110343,
                   0.014180634400617,0.0083326504880713,-0.029052336009585,
                   0.038615085574206,-0.020393486513704,-0.0016554050063734,
                   0.0019955571979541,0.00015870308324157,-1.638856834253E-05,
                   0.043613615723811,0.034994005463765,-0.076788197844621,
                   0.022446277332006,-6.2689710414685E-05,-5.571111856564E-10,
                   -0.19905718354408,0.3177497330738,-0.11841182425981,
                   -31.306260323435,31.546140237781,-2521.3154341695,
                   -0.14874640856724,0.31806110878444])
    alfai = np.array([20.,20.,20.,28.,32.])
    betai = np.array([150.,150.,250.,700.,800.])
    gammai = np.array([1.21,1.21,1.25,0.32,0.32])
    etai = np.array([1.,1.,1.,0.3,0.3])

    Tau_c = Tc / TDBK

    # A: i<7
    mA = np.arange(56) < 7
    mB = (np.arange(56) >= 7) & (np.arange(56) < 51)
    mC = (np.arange(56) >= 51) & (np.arange(56) < 54)
    mD = (np.arange(56) >= 54)

    # Precomputos básicos
    F5 = np.zeros_like(TDBK); Y5 = np.zeros_like(TDBK)
    if np.any(mA):
        F5 = np.sum(ni56[mA] * di[mA] * D**(di[mA]-1.0) * (Tau_c[:,None]**ti[mA]), axis=1)
        Y5 = np.sum(ni56[mA] * di[mA] * ti[mA] * D**(di[mA]-1.0) * (Tau_c[:,None]**(ti[mA]-1.0)), axis=1)

    F6 = np.zeros_like(TDBK); Y6 = np.zeros_like(TDBK)
    if np.any(mB):
        exp_ci = np.exp(-D**ci[mB])
        term = (di[mB] - ci[mB] * D**ci[mB])
        F6 = np.sum(ni56[mB] * exp_ci * (D**(di[mB]-1.0)) * (Tau_c[:,None]**ti[mB]) * term, axis=1)
        Y6 = np.sum(ni56[mB] * ti[mB] * (D**(di[mB]-1.0)) * (Tau_c[:,None]**(ti[mB]-1.0)) * term * exp_ci, axis=1)

    F7 = np.zeros_like(TDBK); Y7 = np.zeros_like(TDBK)
    if np.any(mC):
        ic = np.nonzero(mC)[0]
        for k in ic:
            a = alfai[k-51]; b = betai[k-51]; g = gammai[k-51]
            Fexp = np.exp(-a*(D - etai[k-51])**2 - b*(Tau_c - g)**2)
            F7 += ni56[k] * (D**di[k]) * (Tau_c**ti[k]) * Fexp * (di[k]/D - 2*a*(D - etai[k-51]))
            Y7 += ni56[k] * (D**di[k]) * (Tau_c**ti[k]) * Fexp * ( (ti[k]/Tau_c) - 2*b*(Tau_c - g) ) * (di[k]/D - 2*a*(D - etai[k-51]))

    F8 = np.zeros_like(TDBK); Y8 = np.zeros_like(TDBK)
    if np.any(mD):
        idd = np.nonzero(mD)[0]
        for k in idd:
            g = gammai[k-51]; e = etai[k-51]; a = alfai[k-51]; b = betai[k-51]
            titai = (1 - Tau_c) + g*((D-1)**2)**(1/(2*e))
            deltai = titai**2 + ti[k]*((D-1)**2)**ci[k]
            psii = np.exp(-a*(D-1)**2 - b*(Tau_c-1)**2)
            dpsi_dD = -2*a*(D-1)*psii
            dpsi_dTau = -2*b*(Tau_c-1)*psii
            d2psi_dDTau = 4*a*b*(D-1)*(Tau_c-1)*psii
            ddelta_dD = (D-1)*( g * titai * 2 / e * (((D-1)**2)**(1/(2*e) - 1)) + 2*ti[k]*ci[k]*(((D-1)**2)**(ci[k]-1)) )
            ddeltai_bidD = di[k]*deltai**(di[k]-1) * ddelta_dD
            ddeltai_bidTau = -2 * titai * di[k] * deltai**(di[k]-1)
            d2deltai_bidDTau = -g * di[k] * (2/e) * (D-1) * (((D-1)**2)**(1/(2*e)-1)) \
                               - 2 * titai * di[k]*(di[k]-1) * deltai**(di[k]-2) * ddelta_dD

            F8 += ni56[k] * ( deltai**di[k] * (psii + D*dpsi_dD) + ddeltai_bidD * D * psii )
            Y8 += ni56[k] * ( deltai**di[k] * (dpsi_dTau + D*d2psi_dDTau) + D*ddeltai_bidD*dpsi_dTau
                              + ddeltai_bidTau * (psii + D*dpsi_dD) + d2deltai_bidDTau * D * psii )

    Bww = (F5 + F6 + F7 + F8) / Mcd
    dBwwdT = -(Tau_c**2) * (Y5 + Y6 + Y7 + Y8) / (Mcd * Tc)

    # Cwww (análogos)
    F9 = np.zeros_like(TDBK); Y9 = np.zeros_like(TDBK)
    if np.any(mA):
        F9 = np.sum(ni56[mA] * di[mA]*(di[mA]-1.0) * D**(di[mA]-2.0) * (Tau_c[:,None]**ti[mA]), axis=1)
        Y9 = np.sum(ni56[mA] * di[mA]*(di[mA]-1.0) * ti[mA] * D**(di[mA]-2.0) * (Tau_c[:,None]**(ti[mA]-1.0)), axis=1)

    F10 = np.zeros_like(TDBK); Y10 = np.zeros_like(TDBK)
    if np.any(mB):
        exp_ci = np.exp(-D**ci[mB])
        term = ((di[mB] - ci[mB]*D**ci[mB]) * (di[mB]-1.0 - ci[mB]*D**ci[mB]) - (ci[mB]**2) * D**ci[mB])
        F10 = np.sum(ni56[mB] * D**(di[mB]-2.0) * (Tau_c[:,None]**ti[mB]) * exp_ci * term, axis=1)
        Y10 = np.sum(ni56[mB] * ti[mB] * D**(di[mB]-2.0) * (Tau_c[:,None]**(ti[mB]-1.0)) * exp_ci * term, axis=1)

    F11 = np.zeros_like(TDBK); Y11 = np.zeros_like(TDBK)
    if np.any(mC):
        ic = np.nonzero(mC)[0]
        for k in ic:
            a = alfai[k-51]; b = betai[k-51]; g = gammai[k-51]
            Fexp = np.exp(-a*(D - etai[k-51])**2 - b*(Tau_c - g)**2)
            core = (-2*a*D**di[k] + 4*(a**2)*D**di[k]*(D - etai[k-51])**2 - 4*di[k]*a*D**(di[k]-1)*(D - etai[k-51]) + di[k]*(di[k]-1)*D**(di[k]-2))
            F11 += ni56[k] * (Tau_c**ti[k]) * Fexp * core
            Y11 += ni56[k] * (Tau_c**ti[k]) * Fexp * ((ti[k]/Tau_c) - 2*b*(Tau_c - g)) * core

    F12 = np.zeros_like(TDBK); Y12 = np.zeros_like(TDBK)
    if np.any(mD):
        idd = np.nonzero(mD)[0]
        for k in idd:
            g = gammai[k-51]; e = etai[k-51]; a = alfai[k-51]; b = betai[k-51]
            titai = (1 - Tau_c) + g*((D-1)**2)**(1/(2*e))
            deltai = titai**2 + ti[k]*((D-1)**2)**ci[k]
            psii = np.exp(-a*(D-1)**2 - b*(Tau_c-1)**2)
            dpsi_dD = -2*a*(D-1)*psii
            dpsi_dTau = -2*b*(Tau_c-1)*psii
            d2psi_dD = (2*a*((D-1)**2) - 1) * 2*a * psii
            d2psi_dDTau = 4*a*b*(D-1)*(Tau_c-1)*psii
            d3psi_d2DdTau = -2*(4*a*((D-1)**2) - 2) * a*b*(Tau_c-1)*psii

            ddelta_dD = (D-1)*( g * titai * 2 / e * (((D-1)**2)**(1/(2*e) - 1)) + 2*ti[k]*ci[k]*(((D-1)**2)**(ci[k]-1)) )
            d2delta_dD = (1/(D-1))*ddelta_dD + (D-1)**2 * (
                (g**2) * (2/(e**2)) * (((D-1)**2)**(1/(2*e)-1))**2 +
                g*titai*(4/e)*(((D-1)**2)**(1/(2*e)-2)) + 4*ti[k]*ci[k]*(ci[k]-1)*(((D-1)**2)**(ci[k]-2))
            )
            ddeltai_bidD = di[k]*deltai**(di[k]-1) * ddelta_dD
            d2deltai_bidD = di[k]*(deltai**(di[k]-1)*d2delta_dD + (di[k]-1)*deltai**(di[k]-2)*ddelta_dD**2)
            ddeltai_bidTau = -2 * titai * di[k] * deltai**(di[k]-1)
            d2deltai_bidDTau = -g * di[k] * (2/e) * (D-1) * (((D-1)**2)**(1/(2*e)-1)) - 2*titai*di[k]*(di[k]-1)*deltai**(di[k]-2)*ddelta_dD
            d3deltai_bid2DdTau = di[k]*( d2deltai_bidDTau*d2delta_dD )

            F12 += ni56[k] * ( deltai**di[k]*(2*dpsi_dD + D*d2psi_dD) + 2*ddeltai_bidD*(psii + D*dpsi_dD) + d2deltai_bidD*D*psii )
            Y12 += ni56[k] * ( deltai**di[k]*(2*d2psi_dDTau + D*d3psi_d2DdTau) + ddeltai_bidTau*(2*dpsi_dD + D*d2psi_dD)
                               + 2*d2deltai_bidDTau*(psii + D*dpsi_dD) + 2*ddeltai_bidD*(dpsi_dTau + D*d2psi_dDTau) + d3deltai_bid2DdTau*D*psii + d2deltai_bidD*D*dpsi_dTau )

    Cwww = (F9 + F10 + F11 + F12) / (Mcd**2)
    dCwwwdT = -(Tau_c**2) * (Y9 + Y10 + Y11 + Y12) / (Mcd**2 * Tc)

    # Baw/Caaw/Caww
    ai = np.array([66.5687,-238.834,-176.755])
    bi = np.array([-0.237,-1.048,-3.183])
    tita = TDBK / Tred
    # Horner-like: potencia elemento a elemento
    F13 = ai[0]*tita**bi[0] + ai[1]*tita**bi[1] + ai[2]*tita**bi[2]
    Y13 = ai[0]*bi[0]*tita**(bi[0]-1) + ai[1]*bi[1]*tita**(bi[1]-1) + ai[2]*bi[2]*tita**(bi[2]-1)
    Baw = F13
    dBawdT = Y13 / Tred

    cii = np.array([482.737, 105678.0, -65639400.0, 29444200000.0, -3193170000000.0])
    dii = np.array([-10.72887, 3478.04, -383383.0, 33406000.0])

    # Caaw
    pow_neg = np.stack([TDBK**(-i) for i in range(cii.size)], axis=1)
    F14 = pow_neg @ cii
    idx = np.arange(cii.size)
    # FIX broadcasting: outer power por columnas
    Y14 = ((TDBK[:, None] ** (-(idx + 1))) * (cii * (-idx))[None, :]).sum(axis=1)
    Caaw = F14
    dCaawdT = Y14

    # Caww
    pow_neg_d = np.stack([TDBK**(-i) for i in range(dii.size)], axis=1)
    F15_sum = pow_neg_d @ dii
    Caww = - (1.0 / (0.001**2)) * np.exp(F15_sum)
    j = np.arange(dii.size)
    dF15 = ((TDBK[:, None] ** (-(j + 1))) * (dii * (-j))[None, :]).sum(axis=1)
    Y15 = Caww * dF15
    dCawwdT = Y15

    out = np.stack([Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww, dBaadT, dCaaadT, dBwwdT, dCwwwdT, dBawdT, dCaawdT, dCawwdT], axis=1)
    return _scalarize_like(Tdb, out)

# ---- ISOTCOMP vectorial ----
def ISOTCOMP_vec(Tdb, PATM_kPa):
    Tdb, PATM_kPa = _broadcast_same(Tdb, PATM_kPa)
    TDBK = Tdb + 273.15

    Taux = 1386.0
    paux = 16530.0
    p0 = 101.325
    pt0 = 0.611657
    Tt = 273.16
    R97 = 0.461526

    kT = np.empty_like(TDBK)
    vws = np.empty_like(TDBK)

    mask_liq = TDBK >= 273.15
    mask_ice = ~mask_liq

    if np.any(mask_liq):
        Tau = Taux / TDBK[mask_liq]
        pi = PATM_kPa[mask_liq] / paux
        tau_pow = (Tau[:, None] - 1.222)**Ji
        seven = (7.1 - pi)[:, None]
        dgammadpi = np.sum(-ni * Ii * (seven**(Ii - 1.0)) * tau_pow, axis=1)
        d2gammadpi = np.sum( ni * Ii * (Ii - 1.0) * (seven**(Ii - 2.0)) * tau_pow, axis=1)
        PPa = PATM_kPa[mask_liq] * 1000.0
        kT[mask_liq] = - (1.0/PPa) * pi * d2gammadpi / dgammadpi
        vws[mask_liq] = Mw * R97 * TDBK[mask_liq] * pi * dgammadpi / PPa

    if np.any(mask_ice):
        tita = TDBK[mask_ice] / Tt
        pi = PATM_kPa[mask_ice] / pt0
        pi0 = p0 / pt0
        PT = pt0 * 1000.0
        # g0'
        g0p = np.zeros_like(pi)
        for i in range(1,5):
            g0p += g0[i] * (i/PT) * (pi - pi0)**(i-1)
        r2p = r21*(1/PT) + r22*(2/PT)*(pi - pi0)
        aux1 = r2p * ((t2 - tita)*np.log(t2 - tita) + (t2 + tita)*np.log(t2 + tita) - 2*t2*np.log(t2) - tita**2 / t2)
        dgdp = g0p + Tt * aux1.real
        # g0''
        g0pp = np.zeros_like(pi)
        for i in range(2,5):
            g0pp += g0[i] * i * ((i-1)/(PT**2)) * (pi - pi0)**(i-2)
        r2pp = r22 * 2 / (PT**2)
        aux2 = r2pp * ((t2 - tita)*np.log(t2 - tita) + (t2 + tita)*np.log(t2 + tita) - 2*t2*np.log(t2) - tita**2 / t2)
        d2gdp = g0pp + Tt * aux2.real
        kT[mask_ice] = - d2gdp / dgdp
        vws[mask_ice] = Mw * dgdp

    return np.stack([kT, vws], axis=1)

# ---- HENRYLAW vectorial ----
def HENRYLAW_vec(Tdb, PATM_kPa):
    Tdb, PATM_kPa = _broadcast_same(Tdb, PATM_kPa)
    Tc = 647.096
    psi = np.array([0.7812, 0.2095, 0.0093])
    Ai = np.array([-9.67578, -9.44833, -8.40954])
    Bi = np.array([4.72162, 4.43822, 4.29587])
    Ci = np.array([11.70585, 11.42005, 10.52779])

    TDBK = Tdb + 273.15
    Tr = TDBK / Tc
    Tau = 1 - Tr
    P = PATM_kPa * 1000.0

    pws = pws_T(Tdb)  # Pa
    pws = np.minimum(pws, P)

    betai_vals = []
    for a,b,c in zip(Ai,Bi,Ci):
        betai_vals.append(pws * np.exp(a/Tr + b * (Tau**0.355) / Tr + c * (Tr**-0.41) * np.exp(Tau)))
    betai = np.stack(betai_vals, axis=1)  # (N,3)
    aux = np.sum(psi / betai, axis=1)
    betaa = 1.0 / aux
    bH = 1.0 / (1.01325 * betaa)
    return bH

# ---- ENHAC vectorial ----

def ENHAC_vec(Tdb, PATM_kPa, B_mat=None, iterations=2):
    Tdb, PATM_kPa = _broadcast_same(Tdb, PATM_kPa)
    T_K = Tdb + 273.15
    R = Rm * 1_000_000.0  # J/kmol.K
    P_Pa = PATM_kPa * 1000.0

    if B_mat is None:
        B_mat = VIRCOE_vec(Tdb)
    Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww, dBaadT, dCaaadT, dBwwdT, dCwwwdT, dBawdT, dCaawdT, dCawwdT = B_mat.T

    pws_Pa = pws_T(Tdb)
    pws_Pa = np.minimum(pws_Pa, P_Pa)

    # kT y vws a P y a Ps
    # (fuerza CoolProp IF97/HEOS, evitando 'Incoming out of bound')
    kT_at_P,  v_molar_at_P  = water_kT_and_v(Tdb, PATM_kPa,         prefer="COOLPROP")
    kT_at_Ps, v_molar_at_Ps = water_kT_and_v(Tdb, pws_Pa/1000.0,    prefer="COOLPROP")

    kT = kT_at_P
    vws = v_molar_at_Ps  # m3/kmol a Ps

    # Si estamos saturados, evitar valores no físicos
    kT = np.where(pws_Pa >= P_Pa, 0.0, kT)

    bH = HENRYLAW_vec(Tdb, PATM_kPa)
    X1 = (P_Pa - pws_Pa) / P_Pa
    X1 = np.clip(X1, 0.0, 1.0)
    f = np.ones_like(Tdb)

    Pk = P_Pa/1000.0
    eps = 1e-12
    damp = 0.5

    for _ in range(iterations):
        # --- Términos F1..F15 (sin cambios) ---
        F1  = (((1 + kT * pws_Pa) * (P_Pa - pws_Pa) - kT*(P_Pa**2 - pws_Pa**2)/2.0) * vws) / (R * T_K)

        arg = 1.0 - bH * X1 * Pk            # <-- blindaje de log
        F2  = np.log(np.clip(arg, eps, None))

        F3  = (X1**2 * Pk * Baa) / (R * T_K)
        F4  = -(2 * X1**2 * Pk * Baw) / (R * T_K)
        F5  = -(Bww * (Pk - (X1**2) * Pk)) / (R * T_K)
        F6  = (Caaa * X1**3 * Pk**2) / ((R*T_K)**2)
        F7  = (Caaw * 3 * X1**2 * (1 - 2*X1) * Pk**2) / (2*(R*T_K)**2)
        F8  = -(Caww * 3 * X1**2 * (1 - X1) * Pk**2) / ((R*T_K)**2)
        F9  = -(Cwww * (((1 + 2*X1) * (1 - X1)**2 * Pk**2) - (pws_Pa/1000.0)**2)) / (2*(R*T_K)**2)
        F10 = -(Baa * Bww * X1**2 * (1 - 3*X1) * (1 - X1) * Pk**2) / ((R*T_K)**2)
        F11 = -(Baa * Baw * 2 * X1**3 * (2 - 3*X1) * Pk**2) / ((R*T_K)**2)
        F12 = (Bww * Baw * 6 * X1**2 * (1 - X1)**2 * Pk**2) / ((R*T_K)**2)
        F13 = -(Baa**2 * 3 * X1**4 * Pk**2) / (2*(R*T_K)**2)
        F14 = -(Baw**2 * X1**2 * 2 * (1 - X1) * (1 - 3*X1) * Pk**2) / ((R*T_K)**2)
        F15 = -(Bww**2 * (((pws_Pa/1000.0)**2) - (1 + 3*X1) * (1 - X1)**3 * Pk**2)) / (2*(R*T_K)**2)

        Fsum = (F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9 + F10 + F11 + F12 + F13 + F14 + F15)

        # --- Hardening numérico: saneado y clipping ANTES del exp ---
        Fsum = np.nan_to_num(Fsum, nan=0.0, posinf=0.0, neginf=0.0)
        Fsum = np.clip(Fsum, -50.0, 50.0)

        # --- f_new correctamente definido y luego saneado ---
        f_new = np.exp(Fsum)
        f_new = np.nan_to_num(f_new, nan=1.0, posinf=1.0, neginf=1.0)

        # --- Actualización con damping ---
        f = f + damp * (f_new - f)

        # --- Recalcular X1 y re-clipear ---
        X1 = (P_Pa - f * pws_Pa) / P_Pa
        X1 = np.clip(X1, 0.0, 1.0)

    XAS = X1
    XWS = 1.0 - XAS
    ps  = f * (pws_Pa/1000.0)
    ps  = np.minimum(ps, PATM_kPa - 1e-12)  # evitar división por cero
    WS  = np.where(pws_Pa/1000.0 == PATM_kPa, 1e10, MR * (ps / (PATM_kPa - ps)))

    Ef = np.stack([f, XAS, XWS, WS, pws_Pa/1000.0], axis=1)
    return _scalarize_like(Tdb, Ef)

def MSTAIR_vec(Tdb, PATM_kPa, pws_kPa, XA, XW, W, B_mat, tol=1e-7, maxiter=50):
    Tdb, PATM_kPa, pws_kPa, XA, XW, W = _broadcast_same(Tdb, PATM_kPa, pws_kPa, XA, XW, W)
    T_K = Tdb + 273.15
    R = Rm * 1_000_000.0

    Baa, Caaa, Bww, Cwww, Baw, Caaw, Caww, dBaadT, dCaaadT, dBwwdT, dCwwwdT, dBawdT, dCaawdT, dCawwdT = B_mat.T

    Bm = XA**2 * Baa + 2*XA*XW * Baw + XW**2 * Bww
    Cm = XA**3 * Caaa + 3*XA**2 * XW * Caaw + 3*XA*XW**2 * Caww + XW**3 * Cwww
    dBmdT = XA**2 * dBaadT + 2*XA*XW * dBawdT + XW**2 * dBwwdT
    dCmdT = XA**3 * dCaaadT + 3*XA**2 * XW * dCaawdT + 3*XA*XW**2 * dCawwdT + XW**3 * dCwwwdT

    V = R * T_K / (PATM_kPa * 1000.0)
    for _ in range(maxiter):
        V_new = (R * T_K / (PATM_kPa * 1000.0)) * (1.0 + Bm / V + Cm / (V**2))
        if np.all(np.abs((V_new - V)/V_new) <= tol):
            V = V_new
            break
        V = V_new

    VS = np.where(XA < 1e-7, V / 1000.0, V / (1000.0 * (Ma * XA)))

    # ---- Entalpía ----
    # Aire (Lemmon)
    Ni0 = np.array([6.057194e-8, -2.10274769e-5, -0.000158860716,
                    -13.841928076, 17.275266575, -0.00019536342, 2.490888032,
                    0.791309509, 0.212236768, -0.197938904, 25.36365, 16.90741,
                    87.31279])
    Tj = 132.6312
    Tau = Tj / T_K
    F1 = np.zeros_like(T_K)
    for i in range(5):
        F1 += (i - 3) * Ni0[i] * Tau**(i - 4)
    dalfadTau = F1 + 1.5*Ni0[5]*Tau**0.5 + Ni0[6]*Tau**(-1) \
                 + Ni0[7]*Ni0[10]/(np.exp(Ni0[10]*Tau)-1) \
                 + Ni0[8]*Ni0[11]/(np.exp(Ni0[11]*Tau)-1) \
                 + Ni0[9]*Ni0[12]/((2/3*np.exp(-Ni0[12]*Tau))+1)
    h0lem = -7914.149298
    Rlem = 8.31451
    H1 = h0lem + Rlem * T_K * (1 + Tau * dalfadTau)

    # Agua: IAPWS-97 si T>=273.15, IAPWS-95 si T<273.15
    Tc = 647.096
    Tred = 540.0
    R97u = 8.31451
    R95u = 8.314371

    H2 = np.empty_like(T_K)
    m_liq = T_K >= 273.15
    if np.any(m_liq):
        Ji0 = np.array([0.0,1.0,-5.0,-4.0,-3.0,-2.0,-1.0,2.0,3.0])
        ni097 = np.array([-9.6927686500217, 10.086655968018,
                          -0.005608791128302, 0.071452738081455,
                          -0.40710498223928, 1.4240819171444, -4.383951131945,
                          -0.28408632460772, 0.021268463753307])
        Tau97 = Tred / T_K[m_liq]
        dgammadTau = np.sum(ni097 * Ji0 * (Tau97[:,None]**(Ji0-1.0)), axis=1)
        h097 = -0.01102142797
        H2[m_liq] = h097 + R97u * T_K[m_liq] * Tau97 * dgammadTau
    m_ice = ~m_liq
    if np.any(m_ice):
        ni095 = np.array([-8.3204464837497, 6.6832105275932, 3.00632, 0.012436,
                           0.97315, 1.2795, 0.96956, 0.24873])
        gammai0 = np.array([0.0, 0.0, 0.0, 1.28728967, 3.53734222, 7.74073708,
                             9.24437796, 27.5075105])
        Tau95 = Tc / T_K[m_ice]
        F2 = np.sum(ni095[3:] * gammai0[3:] * (np.exp(gammai0[3:]*Tau95[:,None]) - 1)**-1, axis=1)
        dalfa95dTau = ni095[1] + ni095[2]*Tau95**-1 + F2
        h095 = -0.01102303806
        H2[m_ice] = h095 + R95u * T_K[m_ice] * (1 + Tau95 * dalfa95dTau)

    HS = 2.924425468e-06 + XA*H1 + XW*H2 + Rm * T_K * (((Bm - T_K*dBmdT)/V) + ((Cm - 0.5*T_K*dCmdT)/(V**2)))
    HS = np.where(XA > 1e-7, HS / (Ma * XA), HS)

    out = np.stack([VS, HS], axis=1)
    return _scalarize_like(Tdb, out)
