# -*- coding: utf-8 -*-
"""
Funciones principales para calcular propiedades psicrométricas
----------------------------------------------------------------
Vectorizadas y optimizadas para operar con grandes volúmenes de datos.

Salidas estándar (columnas):
 [Tdb, Twb, RH, W*1000, v, h, Tdp, pw]
  - Tdb (°C)
  - Twb (°C)
  - RH (%)
  - W*1000 (g/kgda)
  - v (m³/kg)
  - h (kJ/kgda)
  - Tdp (°C)
  - pw (Pa)

Autor original: José E. Azzaro
Adaptación vectorizada: M365 Copilot
"""
from __future__ import annotations
import numpy as np

# Núcleo moderno
from .psychro_aux import (
    MR, DEWPT, DEWPT_vec, WETTP, WETTP_vec, HW,
    PWS_IAPWS97, PWS_IAPWS98, HFN2, VFNITR, HFUNIT,
)
from .psychro_aux_vec import ENHAC_vec, MSTAIR_vec, VIRCOE_vec

# ----------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------

def _arr1d(*xs):
    """Convierte entradas a np.ndarray y aplica broadcast a 1-D."""
    arrs = [np.asarray(x, dtype=float) for x in xs]
    b = np.broadcast_arrays(*arrs)
    return [np.ravel(a) for a in b]


def _scalar_case(*xs) -> bool:
    return all(np.ndim(x) == 0 for x in xs)


def _stack_result(Tdb, Twb, RH, W, VS, HS, Tdp, pw):
    out = np.column_stack([
        Tdb, Twb, RH, W * 1000.0, VS, HS, Tdp, pw
    ])
    return out

# ----------------------------------------------------------------------------------
# 1. Psy_TdbRH — propiedades dadas Tdb (°C), RH (%) y PATM (kPa)
# ----------------------------------------------------------------------------------

def Psy_TdbRH(Tdb, RH, PATM):
    scalar = _scalar_case(Tdb, RH, PATM)
    Tdb, RH, P = _arr1d(Tdb, RH, PATM)

    # Vectorizado: viriales, enhancement y mezcla
    B  = VIRCOE_vec(Tdb)
    Ef = ENHAC_vec(Tdb, P, B)  # [f, XAS, XWS, WS, pws_kPa]

    XW = np.clip(Ef[:, 2] * (RH / 100.0), 0.0, 1.0 - 1e-12)
    XA = 1.0 - XW
    W  = MR * XW / np.maximum(XA, 1e-12)
    pw = XW * P * 1000.0  # Pa

    VS_HS = MSTAIR_vec(Tdb, P, Ef[:, 4], XA, XW, W, B)
    VS, HS = VS_HS[:, 0], VS_HS[:, 1]

    Tdp = DEWPT_vec(W, P)
    Twb = WETTP_vec(Tdb, Tdp, HS, W, P)

    out = _stack_result(Tdb, Twb, RH, W, VS, HS, Tdp, pw)
    return out[0] if scalar else out

# ----------------------------------------------------------------------------------
# 2. Psy_TdbW — propiedades dadas Tdb (°C), W (g/kgda) y PATM (kPa)
# ----------------------------------------------------------------------------------

def Psy_TdbW(Tdb, W_gpkg, PATM):
    scalar = _scalar_case(Tdb, W_gpkg, PATM)
    Tdb, W_input, P = _arr1d(Tdb, W_gpkg, PATM)

    W = W_input / 1000.0  # a kg/kgda
    B  = VIRCOE_vec(Tdb)
    Ef = ENHAC_vec(Tdb, P, B)
    WS = Ef[:, 3]

    invalid = W > WS + 1e-12
    if scalar and np.any(invalid):
        return "Datos Incroguentes W > Ws"

    XW = np.clip(W / (W + MR), 0.0, 1.0 - 1e-12)
    XA = 1.0 - XW

    VS_HS = MSTAIR_vec(Tdb, P, Ef[:, 4], XA, XW, W, B)
    VS, HS = VS_HS[:, 0], VS_HS[:, 1]

    Tdp = DEWPT_vec(W, P)
    Twb = WETTP_vec(Tdb, Tdp, HS, W, P)

    RH = np.clip((XW / np.maximum(Ef[:, 2], 1e-12)) * 100.0, 0.0, 100.0)
    pw = XW * P * 1000.0  # Pa

    # Para filas inválidas, llenar con NaN (sólo modo vector)
    if np.any(invalid) and not scalar:
        fill = np.nan
        Tdb[invalid] = Tdb[invalid]  # sin cambio
        Twb[invalid] = fill
        RH[invalid]  = fill
        W[invalid]   = fill
        VS[invalid]  = fill
        HS[invalid]  = fill
        Tdp[invalid] = fill
        pw[invalid]  = fill

    out = _stack_result(Tdb, Twb, RH, W, VS, HS, Tdp, pw)
    return out[0] if scalar and not isinstance(out, str) else out

# ----------------------------------------------------------------------------------
# 3. Psy_Tdbh — propiedades dadas Tdb (°C), h (kJ/kgda) y PATM (kPa)
#      Usa HFN2 (inversor escalar); vectoriza por comprensión para performance razonable.
# ----------------------------------------------------------------------------------

def Psy_Tdbh(Tdb, h, PATM):
    scalar = _scalar_case(Tdb, h, PATM)
    Tdb, H, P = _arr1d(Tdb, h, PATM)

    # Estimación XW inicial desde ENHAC
    Ef = ENHAC_vec(Tdb, P)
    XW0 = np.clip(Ef[:, 2], 1e-6, 1.0 - 1e-6)

    rows = []
    for t, hh, p, xw_init in zip(Tdb, H, P, XW0):
        out = HFN2(t, p, hh, xw_init)  # [XW, XA, W, V]
        XW, XA, W, VS = out
        if (W < 0.0) or (W > 1.0):
            rows.append([np.nan]*8)
            continue
        Tdp = DEWPT(W, p)
        HS  = hh
        Twb = WETTP(t, Tdp, HS, W, p)
        RH  = (XW / max(1e-12, ENHAC_vec(np.atleast_1d(t), np.atleast_1d(p))[0,2])) * 100.0
        pw  = XW * p * 1000.0
        rows.append([t, Twb, RH, W*1000.0, VS, HS, Tdp, pw])

    out = np.array(rows, dtype=float)
    return out[0] if scalar else out

# ----------------------------------------------------------------------------------
# 4. Psy_Tdbv — propiedades dadas Tdb (°C), v (m³/kg) y PATM (kPa)
#      Usa VFNITR (inversor escalar); vectoriza por comprensión.
# ----------------------------------------------------------------------------------

def Psy_Tdbv(Tdb, v, PATM):
    scalar = _scalar_case(Tdb, v, PATM)
    Tdb, V, P = _arr1d(Tdb, v, PATM)

    # XW inicial desde ENHAC
    Ef = ENHAC_vec(Tdb, P)
    XW0 = np.clip(Ef[:, 2], 1e-6, 1.0 - 1e-6)

    rows = []
    for t, vv, p, xw_init in zip(Tdb, V, P, XW0):
        out = VFNITR(t, p, vv, xw_init)  # [XW, XA, W, H]
        if isinstance(out, str):
            rows.append([np.nan]*8)
            continue
        XW, XA, W, HS = out
        if (W < 0.0) or (W > 1.0):
            rows.append([np.nan]*8)
            continue
        Tdp = DEWPT(W, p)
        Twb = WETTP(t, Tdp, HS, W, p)
        RH  = (XW / max(1e-12, ENHAC_vec(np.atleast_1d(t), np.atleast_1d(p))[0,2])) * 100.0
        pw  = XW * p * 1000.0
        rows.append([t, Twb, RH, W*1000.0, vv, HS, Tdp, pw])

    out = np.array(rows, dtype=float)
    return out[0] if scalar else out

# ----------------------------------------------------------------------------------
# 5. Psy_TdbTwb — propiedades dadas Tdb (°C), Twb (°C) y PATM (kPa)
#      Usa HFUNIT (inversor escalar); vectoriza por comprensión.
# ----------------------------------------------------------------------------------

def Psy_TdbTwb(Tdb, Twb, PATM):
    scalar = _scalar_case(Tdb, Twb, PATM)
    Tdb, Twb, P = _arr1d(Tdb, Twb, PATM)

    rows = []
    for t, twb, p in zip(Tdb, Twb, P):
        # Estado a Twb para construir FUN (H - W*HW)
        Bwb = VIRCOE_vec(np.atleast_1d(twb))[0]
        Efwb = ENHAC_vec(np.atleast_1d(twb), np.atleast_1d(p), np.atleast_2d(Bwb))[0]
        PSAT = (PWS_IAPWS97(twb) if twb >= 0.0 else PWS_IAPWS98(twb)) / 1000.0  # kPa
        WW = Efwb[0] * PSAT
        WS = MR * WW / max(1e-12, (p - WW))
        VS_HS = MSTAIR_vec(np.atleast_1d(twb), np.atleast_1d(p), np.atleast_1d(Efwb[4]),
                           np.atleast_1d(Efwb[1]), np.atleast_1d(Efwb[2]), np.atleast_1d(WS),
                           np.atleast_2d(Bwb))[0]
        HS_wb = VS_HS[1]
        HWWB  = HW(twb, p)
        FUN   = HS_wb - WS * HWWB

        # XW inicial: promedio XWS(db, wb)
        Bdb = VIRCOE_vec(np.atleast_1d(t))[0]
        Efdb = ENHAC_vec(np.atleast_1d(t), np.atleast_1d(p), np.atleast_2d(Bdb))[0]
        XW0  = max(1e-6, min(1.0 - 1e-6, 0.5 * (Efdb[2] + Efwb[2])))

        HF = HFUNIT(t, p, FUN, HWWB, XW0)  # [XW, XA, W, H, V]
        if isinstance(HF, str):
            rows.append([np.nan]*8)
            continue
        XW, XA, W, HS, VS = HF
        Tdp = DEWPT(W, p)
        RH  = (XW / max(1e-12, Efdb[2])) * 100.0
        pw  = XW * p * 1000.0
        rows.append([t, twb, RH, W*1000.0, VS, HS, Tdp, pw])

    out = np.array(rows, dtype=float)
    return out[0] if scalar else out

if __name__ == '__main__':
    # Smoke tests rápidos
    print(Psy_TdbRH(23.2547847501, 60.0, 101.325))
    print(Psy_TdbW(27.0, 10.69377158499, 101.325))
