# -*- coding: utf-8 -*-
"""
psychro_aux.py — Versión fusionada y modernizada (1 → 23)
---------------------------------------------------------
- Integra las funciones originales con la capa vectorizada robusta
  (ENHAC_vec, MSTAIR_vec, VIRCOE_vec) y el backend endurecido
  (IF97/HEOS con fallbacks) ya trabajados.
- Mantiene firmas y retornos compatibles con tu implementación histórica
  siempre que fue posible, pero delega a los módulos vectoriales/ backend
  para mayor estabilidad numérica.

Autor: José Azzaro + M365 Copilot
"""
from __future__ import annotations
import numpy as np
from typing import Tuple
from scipy.optimize import root_scalar

# ----------------- Constantes y backend -----------------
Rm = 8.314472     # kJ/kmol.K
Ma = 28.966       # kg/kmol
Mw = 18.015268    # kg/kmol
MR = Mw / Ma      # Mw/Ma

# Backend de agua y saturación (endurecido)
from .iapws_backend import (pws_T, tws_P, water_enthalpy_liq, water_kT_and_v)
from .psychro_aux_vec import (VIRCOE_vec, ENHAC_vec, MSTAIR_vec, ISOTCOMP_vec, HENRYLAW_vec)

# Vectoriales (endurecidas)
from .psychro_aux_vec import (
    VIRCOE_vec, ENHAC_vec, MSTAIR_vec, ISOTCOMP_vec, HENRYLAW_vec,
)

# ----------------- Utilitarios internos -----------------

def _clip01(x: float, eps: float = 1e-8) -> float:
    return float(np.clip(x, eps, 1.0 - eps))


def _ef_at(Tdb: float, PATM: float):
    """Fila Ef para (Tdb, PATM): [f, XAS, XWS, WS, pws_kPa]."""
    Ef = ENHAC_vec(np.atleast_1d(Tdb), np.atleast_1d(PATM))
    return Ef[0]


def _B_at(Tdb: float):
    """Fila B (viriales y derivadas) para Tdb."""
    B = VIRCOE_vec(np.atleast_1d(Tdb))
    return B[0]


def _mstairs_at(Tdb: float, PATM: float, Ef_row, XA: float, XW: float, W: float, B_row):
    """Retorna (VS, HS) kJ/kgda usando la vectorial robusta."""
    VS_HS = MSTAIR_vec(
        np.atleast_1d(Tdb),
        np.atleast_1d(PATM),
        np.atleast_1d(Ef_row[4]),  # pws_kPa
        np.atleast_1d(XA),
        np.atleast_1d(XW),
        np.atleast_1d(W),
        np.atleast_2d(B_row)
    )
    return VS_HS[0]  # (VS, HS)

# --------------------------------------------------------
# 1) ICAO — presión barométrica vs altitud
# --------------------------------------------------------

def patm_ICAO(H):
    P0, P1, P2 = 101.325, 2.25577, 5.2559
    if H > 11000 or H < -5000:
        return "La Altura no puede ser <-5000 ni >11000"
    return P0 * (1 - P1 * 1e-5 * H) ** P2

# --------------------------------------------------------
# 2) PWS_IAPWS97 — presión de saturación (Pa) (LV, T>=0°C)
# 2b) pws_Aprox — aproximación Magnus-Tetens
# 2c) dpws_Aprox — derivada aproximada
# --------------------------------------------------------

def PWS_IAPWS97(Tdb):
    return float(pws_T(Tdb))


def pws_Aprox(Tdb):
    import math
    return 610.78 * math.exp(17.2694 * Tdb / (Tdb + 237.3))


def dpws_Aprox(patm, tsref):
    import math
    aux = math.exp(17.2694 * tsref / (tsref + 237.3))
    return MR * patm * 2.502904E6 * aux / (((tsref + 273.3) ** 2) * (patm - 610.78 * aux))

# --------------------------------------------------------
# 3) TWS_IAPWS97 — temperatura de saturación (K) vs ps (Pa)
# --------------------------------------------------------

def TWS_IAPWS97(ps):
    return float(tws_P(ps))

# --------------------------------------------------------
# 4) PWS_IAPWS98 — sublimación (hielo, T<0°C)
#    (compatibilidad: devolvemos pws_T si T<0; si T>=0 devolvemos aviso)
# --------------------------------------------------------

def PWS_IAPWS98(Tdb):
    if Tdb < 0.0:
        return float(pws_T(Tdb))
    return "Temperatura fuera de rango >= 273.15°K"

# --------------------------------------------------------
# 5) VIRCOE — coeficientes viriales (fila)
# --------------------------------------------------------

def VIRCOE(Tdb):
    return _B_at(Tdb)

# --------------------------------------------------------
# 6) ISOTCOMP — compresibilidad isotérmica (1/Pa) y vws (m3/kmol a Ps)
# --------------------------------------------------------

def ISOTCOMP(Tdb, PATM):
    # k_T a P (kPa) y volumen molar a Ps(T)
    kT, _ = water_kT_and_v(Tdb, PATM)               # k_T (1/Pa), v_m (m3/kmol) a P
    Ef = _ef_at(Tdb, PATM)
    Ps = Ef[4]  # pws_kPa
    _, v_m_Ps = water_kT_and_v(Tdb, Ps)             # v_m en saturación
    return np.array([float(np.atleast_1d(kT)[0]), float(np.atleast_1d(v_m_Ps)[0])])

# --------------------------------------------------------
# 7) HENRYLAW — constante de Henry (fila-> escalar)
# --------------------------------------------------------

def HENRYLAW(Tdb, PATM):
    return float(HENRYLAW_vec(np.atleast_1d(Tdb), np.atleast_1d(PATM))[0])

# --------------------------------------------------------
# 8) ENHAC — factor de mejoramiento (fila escalar)
# --------------------------------------------------------

def ENHAC(Tdb, PATM, B=None):
    if B is None:
        Ef = ENHAC_vec(np.atleast_1d(Tdb), np.atleast_1d(PATM))
    else:
        Ef = ENHAC_vec(np.atleast_1d(Tdb), np.atleast_1d(PATM), np.atleast_2d(B))
    return Ef[0]

# --------------------------------------------------------
# 9) MSTAIR — (VS, HS) para mezcla
# --------------------------------------------------------

def MSTAIR(Tdb, PATM, pws, XA, XW, W, B):
    VS_HS = MSTAIR_vec(
        np.atleast_1d(Tdb), np.atleast_1d(PATM), np.atleast_1d(pws),
        np.atleast_1d(XA), np.atleast_1d(XW), np.atleast_1d(W), np.atleast_2d(B)
    )
    return VS_HS[0]

# --------------------------------------------------------
# 10) DEWPT — wrappers a versión vectorizada
# --------------------------------------------------------

def DEWPT_vec(W, PATM_kPa, tol=1e-4, maxiter=50):
    W = np.asarray(W, dtype=float)
    P = np.asarray(PATM_kPa, dtype=float)
    W, P = np.broadcast_arrays(W, P)
    out = np.full_like(W, np.nan, dtype=float)

    zero_mask = (W <= 0.0)
    out[zero_mask] = 0.0
    m = ~zero_mask
    if not np.any(m):
        return float(out) if out.shape == () else out

    # Semilla (Peppers / Magnus-Tetens)
    A = np.log(np.clip(P[m] * W[m] / (MR + W[m]), 1e-12, None))
    T0 = (6.54 + 14.526*A + 0.7389*A*A + 0.09486*A*A*A
          + 0.4569 * np.power(np.clip(P[m] * W[m] / (MR + W[m]), 1e-12, None), 0.1984))
    T0 = np.where(T0 < 0.0, 6.09 + 12.608*A + 0.4959*A*A, T0)
    T = T0.copy(); T1 = T0 + 0.5
    T_min, T_max = -100.0, 100.0

    def Wsat_at(T_C, P_kPa):
        Ef = ENHAC_vec(T_C, P_kPa)
        return Ef[:, 3]  # WS

    F0 = Wsat_at(T, P[m]) - W[m]
    F1 = Wsat_at(T1, P[m]) - W[m]
    converged = (np.abs(F1) < tol)

    for _ in range(maxiter):
        if np.all(converged):
            break
        denom = (F1 - F0)
        denom = np.where(np.abs(denom) < 1e-16, 1e-16, denom)
        T_new = T1 - F1 * (T1 - T) / denom
        T_new = np.clip(T_new, T_min, T_max)
        upd = ~converged
        T[upd], T1[upd] = T1[upd], T_new[upd]
        F0[upd], F1[upd] = F1[upd], (Wsat_at(T1[upd], P[m][upd]) - W[m][upd])
        converged = (np.abs(F1) < tol)

    out[m] = T1
    return float(out) if out.shape == () else out


def DEWPT(W, PATM):
    return float(DEWPT_vec(W, PATM))

# --------------------------------------------------------
# 11) HW — entalpía específica agua (kJ/kg), delegada en backend
# --------------------------------------------------------

def HW(Tdb, PATM):
    return float(water_enthalpy_liq(Tdb, PATM))

# --------------------------------------------------------
# 12) WETTP — wrappers a versión vectorizada
# --------------------------------------------------------

def WETTP_vec(Tdb_C, Tdp_C, H_kJkgda, W, PATM_kPa, tol=1e-4, maxiter=50):
    """Bulbo húmedo (°C) resolviendo balance de entalpía.
    - Acepta escalares o arrays; si la llamada fue escalar, retorna float.
    - Mantiene el mismo título/firmas que tu versión actual.
    """
    import numpy as _np

    Tdb = _np.asarray(Tdb_C, dtype=float)
    Tdp = _np.asarray(Tdp_C, dtype=float)
    H   = _np.asarray(H_kJkgda, dtype=float)
    W   = _np.asarray(W, dtype=float)
    P   = _np.asarray(PATM_kPa, dtype=float)

    # Detectar caso escalar puro ANTES del broadcast
    scalar_case = (Tdb.ndim == 0 and Tdp.ndim == 0 and H.ndim == 0 and W.ndim == 0 and P.ndim == 0)

    # Broadcast y coerción a 1-D para usar máscaras sin problemas
    Tdb, Tdp, H, W, P = _np.broadcast_arrays(Tdb, Tdp, H, W, P)
    Tdb = Tdb.reshape(-1)
    Tdp = Tdp.reshape(-1)
    H   = H.reshape(-1)
    W   = W.reshape(-1)
    P   = P.reshape(-1)

    T0 = (Tdb + Tdp) * 0.5
    T1 = T0 + 0.5
    T_min = _np.minimum(Tdp, Tdb) - 5.0
    T_max = _np.maximum(Tdp, Tdb) + 5.0

    def eval_fun(T_C, P_kPa, W_val, H_val):
        # Mantener tu enfoque por bucle; robusto para 1-D y evita indexado escalar
        F = _np.empty_like(T_C, dtype=float)
        for i, (t, p, w, h) in enumerate(zip(T_C, P_kPa, W_val, H_val)):
            B  = VIRCOE(t)
            Ef = ENHAC(t, p, B)
            WS = Ef[3]
            HS = MSTAIR(t, p, Ef[4], Ef[1], Ef[2], WS, B)[1]
            HWWB = HW(t, p)
            F[i] = h - HS - (w - WS) * HWWB
        return F

    F0 = eval_fun(T0, P, W, H)
    F1 = eval_fun(T1, P, W, H)
    converged = _np.abs(F1) < tol

    for _ in range(maxiter):
        if _np.all(converged):
            break
        denom = (F1 - F0)
        denom = _np.where(_np.abs(denom) < 1e-16, 1e-16, denom)
        T_new = T1 - F1 * (T1 - T0) / denom
        T_new = _np.minimum(_np.maximum(T_new, T_min), T_max)
        upd = ~converged
        T0[upd], T1[upd] = T1[upd], T_new[upd]
        F0[upd], F1[upd] = F1[upd], eval_fun(T1[upd], P[upd], W[upd], H[upd])
        converged = (_np.abs(F1) < tol)

    TWB = T1
    return float(TWB[0]) if scalar_case else TWB
def WETTP(Tdb, Tdp, H, W, PATM):
    return float(WETTP_vec(Tdb, Tdp, H, W, PATM))

# --------------------------------------------------------
# 13) WETTP2 — compatibilidad: delega en WETTP
# --------------------------------------------------------

def WETTP2(Tdb, H, W, PATM):
    # compat: estimar Tdp desde W y delegar en WETTP
    Tdp = DEWPT(W, PATM)
    return float(WETTP(Tdb, Tdp, H, W, PATM))

# --------------------------------------------------------
# 14) HFN2 — [XW, XA, W, V] dado (Tdb, PATM, H, XW0)
# --------------------------------------------------------

def HFN2(Tdb, PATM, H, XW):
    XW = _clip01(float(XW))
    Ef_row = _ef_at(Tdb, PATM)
    B_row  = _B_at(Tdb)

    def fun(xw):
        xw = _clip01(xw)
        xa = 1.0 - xw
        w  = MR * xw / xa
        VS, HS = _mstairs_at(Tdb, PATM, Ef_row, xa, xw, w, B_row)
        return HS - H

    x0 = XW
    x1 = _clip01(XW * 1.02 + 1e-4)
    try:
        r = root_scalar(fun, method="secant", x0=x0, x1=x1)
        xw_sol = _clip01(r.root)
    except Exception:
        grid = np.clip(np.linspace(1e-6, 1-1e-6, 21), 1e-6, 1-1e-6)
        vals = [fun(xx) for xx in grid]
        sgn  = np.sign(vals)
        brk = None
        for i in range(len(grid)-1):
            if sgn[i] == 0.0:
                brk = (grid[i], grid[i]); break
            if sgn[i]*sgn[i+1] < 0:
                brk = (grid[i], grid[i+1]); break
        if brk is not None and brk[0] != brk[1]:
            r = root_scalar(fun, method="brentq", bracket=brk)
            xw_sol = _clip01(r.root)
        else:
            xw_sol = float(grid[int(np.argmin(np.abs(vals)))])

    xa_sol = 1.0 - xw_sol
    w_sol  = MR * xw_sol / xa_sol
    VS, HS = _mstairs_at(Tdb, PATM, Ef_row, xa_sol, xw_sol, w_sol, B_row)
    return np.array([xw_sol, xa_sol, w_sol, VS])

# --------------------------------------------------------
# 15) VFNITR — [XW, XA, W, H] dado (Tdb, PATM, V, XW0)
# --------------------------------------------------------

def VFNITR(Tdb, PATM, V, XW):
    XW = _clip01(float(XW))
    XA = 1.0 - XW
    W  = MR * XW / XA

    Ef_row = _ef_at(Tdb, PATM)
    B_row  = _B_at(Tdb)

    def g(xw):
        xw = _clip01(xw)
        xa = 1.0 - xw
        w  = MR * xw / xa
        VS, HS = _mstairs_at(Tdb, PATM, Ef_row, xa, xw, w, B_row)
        return VS - V

    x0 = XW
    x1 = _clip01(XW * 1.02 + 1e-4)
    try:
        r = root_scalar(g, method="secant", x0=x0, x1=x1)
        xw_sol = _clip01(r.root)
    except Exception:
        grid = np.clip(np.linspace(1e-6, 1-1e-6, 21), 1e-6, 1-1e-6)
        vals = [g(xx) for xx in grid]
        sgn = np.sign(vals)
        brk = None
        for i in range(len(grid)-1):
            if sgn[i] == 0.0:
                brk = (grid[i], grid[i]); break
            if sgn[i]*sgn[i+1] < 0:
                brk = (grid[i], grid[i+1]); break
        if brk is not None and brk[0] != brk[1]:
            r = root_scalar(g, method="brentq", bracket=brk)
            xw_sol = _clip01(r.root)
        else:
            xw_sol = float(grid[int(np.argmin(np.abs(vals)))])

    xa_sol = 1.0 - xw_sol
    w_sol  = MR * xw_sol / xa_sol
    VS, HS = _mstairs_at(Tdb, PATM, Ef_row, xa_sol, xw_sol, w_sol, B_row)
    return np.array([xw_sol, xa_sol, w_sol, HS])

# --------------------------------------------------------
# 16) HFUNIT — [XW, XA, W, H, V] dado (Tdb, PATM, FUN, HW, XW0)
# --------------------------------------------------------

def HFUNIT(Tdb, PATM, FUN, HW_val, XW):
    XW = _clip01(float(XW))
    Ef_row = _ef_at(Tdb, PATM)
    B_row  = _B_at(Tdb)

    def fc(xw):
        xw = _clip01(xw)
        xa = 1.0 - xw
        w  = MR * xw / xa
        VS, HS = _mstairs_at(Tdb, PATM, Ef_row, xa, xw, w, B_row)
        return (HS - w * HW_val) - FUN

    x0 = XW
    x1 = _clip01(XW * 1.05 + 1e-4)
    try:
        r = root_scalar(fc, method="secant", x0=x0, x1=x1)
        xw_sol = _clip01(r.root)
    except Exception:
        grid = np.clip(np.linspace(1e-6, 1-1e-6, 21), 1e-6, 1-1e-6)
        vals = [fc(xx) for xx in grid]
        sgn  = np.sign(vals)
        brk = None
        for i in range(len(grid)-1):
            if sgn[i] == 0.0:
                brk = (grid[i], grid[i]); break
            if sgn[i]*sgn[i+1] < 0:
                brk = (grid[i], grid[i+1]); break
        if brk is not None and brk[0] != brk[1]:
            r = root_scalar(fc, method="brentq", bracket=brk)
            xw_sol = _clip01(r.root)
        else:
            xw_sol = float(grid[int(np.argmin(np.abs(vals)))])

    xa_sol = 1.0 - xw_sol
    w_sol  = MR * xw_sol / xa_sol
    VS, HS = _mstairs_at(Tdb, PATM, Ef_row, xa_sol, xw_sol, w_sol, B_row)
    return np.array([xw_sol, xa_sol, w_sol, HS, VS])

# --------------------------------------------------------
# 17) TDBITb — Tdb dado (W, H, PATM)
# --------------------------------------------------------

def TDBITb(W, H, PATM):
    W = float(W)
    T0 = H - 2501.0*W / (1.006 + 1.86*W)
    T0 = float(np.clip(T0, -60.0, 80.0))

    def F(T):
        Ef = _ef_at(T, PATM)
        XW = np.clip(Ef[2], 1e-12, 1-1e-12)
        XA = 1.0 - XW
        VS, HS = _mstairs_at(T, PATM, Ef, XA, XW, W, _B_at(T))
        return H - HS

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

# --------------------------------------------------------
# 18) TDB2ITb — Tdb dado (RH, H, PATM)
# --------------------------------------------------------

def TDB2ITb(RH, H, PATM):
    RH = float(RH)
    T0 = 20.0

    def F(T):
        Ef = _ef_at(T, PATM)
        XW_sat = np.clip(Ef[2], 1e-12, 1-1e-12)
        XW = np.clip(XW_sat * RH/100.0, 1e-12, 1-1e-12)
        XA = 1.0 - XW
        W  = MR * XW / XA
        VS, HS = _mstairs_at(T, PATM, Ef, XA, XW, W, _B_at(T))
        return H - HS

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

# --------------------------------------------------------
# 19) TDB3ITb — Tdb dado (W, V, PATM)
# --------------------------------------------------------

def TDB3ITb(W, V, PATM):
    T0 = PATM*1000.0*V / (Rm*(1.0 + 1.608*W)) - 273.15
    T0 = float(np.clip(T0, -60.0, 80.0))

    def F(T):
        Ef = _ef_at(T, PATM)
        XW = np.clip(W / (MR + W), 1e-12, 1-1e-12)
        XA = 1.0 - XW
        VS, HS = _mstairs_at(T, PATM, Ef, XA, XW, W, _B_at(T))
        return V - VS

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

# --------------------------------------------------------
# 20) TDB4ITb — Tdb dado (RH, V, PATM)
# --------------------------------------------------------

def TDB4ITb(RH, V, PATM):
    T0 = 20.0

    def F(T):
        Ef = _ef_at(T, PATM)
        XW_sat = np.clip(Ef[2], 1e-12, 1-1e-12)
        XW = np.clip(XW_sat * RH/100.0, 1e-12, 1-1e-12)
        XA = 1.0 - XW
        W  = MR * XW / XA
        VS, HS = _mstairs_at(T, PATM, Ef, XA, XW, W, _B_at(T))
        return V - VS

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

# --------------------------------------------------------
# 21) TDB5ITb — Tdb dado (W, Twb, PATM)
# --------------------------------------------------------

def TDB5ITb(W, Twb, PATM):
    T0 = float(Twb)

    def F(T):
        Ef = _ef_at(T, PATM)
        XW = np.clip(W / (MR + W), 1e-12, 1-1e-12)
        XA = 1.0 - XW
        VS, HS = _mstairs_at(T, PATM, Ef, XA, XW, W, _B_at(T))
        Tdp = DEWPT(W, PATM)
        Twb_calc = WETTP(T, Tdp, HS, W, PATM)
        return Twb - Twb_calc

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

# --------------------------------------------------------
# 22) TDB6ITb — Tdb dado (RH, Twb, PATM)
# --------------------------------------------------------

def TDB6ITb(RH, Twb, PATM):
    T0 = float(Twb)

    def F(T):
        Ef = _ef_at(T, PATM)
        XW_sat = np.clip(Ef[2], 1e-12, 1-1e-12)
        XW = np.clip(XW_sat * RH/100.0, 1e-12, 1-1e-12)
        XA = 1.0 - XW
        W  = MR * XW / XA
        VS, HS = _mstairs_at(T, PATM, Ef, XA, XW, W, _B_at(T))
        Tdp = DEWPT(W, PATM)
        Twb_calc = WETTP(T, Tdp, HS, W, PATM)
        return Twb - Twb_calc

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

# --------------------------------------------------------
# 23) TDB7ITb — Tdb dado (W, RH, PATM)
# --------------------------------------------------------

def TDB7ITb(W, RH, PATM):
    try:
        _ = W * PATM / ((RH / 100.0) * (MR + W))  # solo para validar
        T0 = 20.0
    except Exception:
        T0 = 20.0

    def F(T):
        Ef = _ef_at(T, PATM)
        f  = Ef[0]
        pws_kPa = Ef[4]
        lhs = f * pws_kPa
        rhs = W * PATM / ((RH/100.0) * (MR + W))
        return lhs - rhs

    x0 = T0; x1 = x0 + 0.5
    r = root_scalar(F, method="secant", x0=x0, x1=x1)
    return float(r.root)

if __name__ == "__main__":
    # Pequeña smoke-test opcional
    Tdb, P = 25.0, 101.325
    Ef = ENHAC(Tdb, P)
    print("ENHAC:", Ef)
    print("DEWPT(0.008,101.325):", DEWPT(0.008, 101.325))
    HS_test = MSTAIR(Tdb, P, Ef[4], Ef[1], Ef[2], Ef[3], VIRCOE(Tdb))[1]
    print("WETTP(25, Tdp, H, W, P):", WETTP(25.0, DEWPT(0.008, P), HS_test, 0.008, P))
