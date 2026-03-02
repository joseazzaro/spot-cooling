# -*- coding: utf-8 -*-
# Coeficientes de las regresiones T_a(0.5) = a * T_mr + b
# extraidos de la Tabla 1 (Azer, ASHRAE Part 1).
# Las pendientes 'a' son negativas (T_a(0.5) desciende al crecer T_mr).
# Estructura: TABLE[Icl][M][V] = (a, b)

from typing import Dict, Tuple

ICL_LEVELS = [0.3, 0.6, 0.9]
M_LEVELS = [87.0, 174.0, 232.0]
V_LEVELS = [0.5, 1.0, 1.5, 2.0]

TABLE: Dict[float, Dict[float, Dict[float, Tuple[float, float]]]] = {
    0.3: {
        87.0: {0.5: (-0.293, 48.1), 1.0: (-0.185, 46.2), 1.5: (-0.143, 45.7), 2.0: (-0.118, 45.2)},
        174.0: {0.5: (-0.263, 35.9), 1.0: (-0.211, 39.4), 1.5: (-0.136, 38.4), 2.0: (-0.118, 39.2)},
        232.0: {0.5: (-0.385, 32.6), 1.0: (-0.193, 32.8), 1.5: (-0.118, 32.9), 2.0: (-0.118, 35.0)},
    },
    0.6: {
        87.0: {0.5: (-0.277, 45.9), 1.0: (-0.158, 43.3), 1.5: (-0.129, 43.3), 2.0: (-0.118, 43.3)},
        174.0: {0.5: (-0.291, 32.0), 1.0: (-0.139, 31.5), 1.5: (-0.118, 32.8), 2.0: (-0.118, 34.9)},
        232.0: {0.5: (-0.505, 30.0), 1.0: (-0.248, 28.7), 1.5: (-0.185, 29.8), 2.0: (-0.118, 28.5)},
    },
    0.9: {
        87.0: {0.5: (-0.248, 42.5), 1.0: (-0.191, 43.3), 1.5: (-0.118, 40.9), 2.0: (-0.100, 41.1)},
        174.0: {0.5: (-0.467, 34.0), 1.0: (-0.227, 31.2), 1.5: (-0.139, 29.9), 2.0: (-0.118, 30.2)},
        232.0: {1.0: (-0.361, 25.9), 1.5: (-0.229, 24.4), 2.0: (-0.216, 26.4)},
    }
}

def _bracket(value, grid):
    """Devuelve (lo, hi, t) donde t en [0,1] para interpolacion lineal.
    Si value esta fuera del rango, se satura en extremos (t=0).
    """
    if value <= grid[0]:
        return grid[0], grid[0], 0.0
    if value >= grid[-1]:
        return grid[-1], grid[-1], 0.0
    for i in range(len(grid)-1):
        lo, hi = grid[i], grid[i+1]
        if lo <= value <= hi:
            t = 0.0 if hi == lo else (value - lo) / (hi - lo)
            return lo, hi, t
    return grid[-2], grid[-1], 1.0

def regress_Ta50_coeffs(Icl: float, M: float, V: float):
    ic_lo, ic_hi, t_ic = _bracket(Icl, ICL_LEVELS)
    m_lo, m_hi, t_m = _bracket(M, M_LEVELS)
    v_lo, v_hi, t_v = _bracket(V, V_LEVELS)

    def ab_at(ic, mm, vv):
        table_ic = TABLE.get(ic, {})
        table_m = table_ic.get(mm, {})
        if vv in table_m:
            return table_m[vv]
        # vecino mas cercano en V
        if not table_m:
            raise ValueError('Faltan datos en la tabla para esta combinacion')
        vv2, ab = min(table_m.items(), key=lambda kv: abs(kv[0]-vv))
        return ab

    def lerp(a, b, t):
        return a + (b-a)*t

    # Ic = ic_lo
    a_ll, b_ll = ab_at(ic_lo, m_lo, v_lo)
    a_lh, b_lh = ab_at(ic_lo, m_lo, v_hi)
    a_hl, b_hl = ab_at(ic_lo, m_hi, v_lo)
    a_hh, b_hh = ab_at(ic_lo, m_hi, v_hi)
    a_lm = lerp(a_ll, a_lh, t_v); b_lm = lerp(b_ll, b_lh, t_v)
    a_hm = lerp(a_hl, a_hh, t_v); b_hm = lerp(b_hl, b_hh, t_v)
    a_im = lerp(a_lm, a_hm, t_m); b_im = lerp(b_lm, b_hm, t_m)

    # Ic = ic_hi
    a_ll, b_ll = ab_at(ic_hi, m_lo, v_lo)
    a_lh, b_lh = ab_at(ic_hi, m_lo, v_hi)
    a_hl, b_hl = ab_at(ic_hi, m_hi, v_lo)
    a_hh, b_hh = ab_at(ic_hi, m_hi, v_hi)
    a_lm2 = lerp(a_ll, a_lh, t_v); b_lm2 = lerp(b_ll, b_lh, t_v)
    a_hm2 = lerp(a_hl, a_hh, t_v); b_hm2 = lerp(b_hl, b_hh, t_v)
    a_im2 = lerp(a_lm2, a_hm2, t_m); b_im2 = lerp(b_lm2, b_hm2, t_m)

    a = lerp(a_im, a_im2, t_ic)
    b = lerp(b_im, b_im2, t_ic)
    return a, b
