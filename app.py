# -*- coding: utf-8 -*-
# Spot Cooling Designer - single file
import sys, math
from PyQt5 import QtWidgets, QtCore, QtGui, QtPrintSupport
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

# Import psychrometric functions from psychro module
from psychro.psychro_fun import Psy_TdbRH

# psychro utils
import math

def psat_mmhg(T):
    return math.exp(18.6686 - 4030.183/(T + 235.0))

def pv_from_rh_T(rh, T):
    return rh * psat_mmhg(T)

def rh_from_pv_T(Pv, T):
    return Pv / psat_mmhg(T)

MMHG_TO_KPA = 0.133322368

def humidity_ratio_from_pv(Pv_mmhg, p_atm_kpa=101.325):
    Pv_kpa = Pv_mmhg * MMHG_TO_KPA
    return 0.62198 * Pv_kpa / max(1e-6, (p_atm_kpa - Pv_kpa))

def pv_from_humidity_ratio(w, p_atm_kpa=101.325):
    Pv_kpa = (w * p_atm_kpa) / max(1e-9, (0.62198 + w))
    return Pv_kpa / MMHG_TO_KPA

def moist_air_enthalpy_kJkg(T, w):
    return 1.006*T + w*(2501.0 + 1.86*T)

def air_density_approx(T, w, p_atm_kpa=101.325):
    T_k = T + 273.15
    P_pa = p_atm_kpa*1000.0
    R_d = 287.055
    return P_pa/(R_d*T_k*(1.0 + 1.607*w))

# regression table
ICL_LEVELS = [0.3, 0.6, 0.9]
M_LEVELS   = [87.0, 174.0, 232.0]
V_LEVELS   = [0.5, 1.0, 1.5, 2.0]

TABLE = {
    0.3:{
        87.0:  {0.5:(-0.293,48.1), 1.0:(-0.185,46.2), 1.5:(-0.143,45.7), 2.0:(-0.118,45.2)},
        174.0: {0.5:(-0.263,35.9), 1.0:(-0.211,39.4), 1.5:(-0.136,38.4), 2.0:(-0.118,39.2)},
        232.0: {0.5:(-0.385,32.6), 1.0:(-0.193,32.8), 1.5:(-0.118,32.9), 2.0:(-0.118,35.0)}
    },
    0.6:{
        87.0:  {0.5:(-0.277,45.9), 1.0:(-0.158,43.3), 1.5:(-0.129,43.3), 2.0:(-0.118,43.3)},
        174.0: {0.5:(-0.291,32.0), 1.0:(-0.139,31.5), 1.5:(-0.118,32.8), 2.0:(-0.118,34.9)},
        232.0: {0.5:(-0.505,30.0), 1.0:(-0.248,28.7), 1.5:(-0.185,29.8), 2.0:(-0.118,28.5)}
    },
    0.9:{
        87.0:  {0.5:(-0.248,42.5), 1.0:(-0.191,43.3), 1.5:(-0.118,40.9), 2.0:(-0.100,41.1)},
        174.0: {0.5:(-0.467,34.0), 1.0:(-0.227,31.2), 1.5:(-0.139,29.9), 2.0:(-0.118,30.2)},
        232.0: {1.0:(-0.361,25.9), 1.5:(-0.229,24.4), 2.0:(-0.216,26.4)}
    }
}

def _bracket(value, grid):
    if value <= grid[0]: return grid[0], grid[0], 0.0
    if value >= grid[-1]: return grid[-1], grid[-1], 0.0
    for i in range(len(grid)-1):
        lo, hi = grid[i], grid[i+1]
        if lo <= value <= hi:
            t = 0.0 if hi==lo else (value-lo)/(hi-lo)
            return lo, hi, t
    return grid[-2], grid[-1], 1.0

def regress_Ta50_coeffs(Icl, M, V):
    ic_lo, ic_hi, t_ic = _bracket(Icl, ICL_LEVELS)
    m_lo, m_hi, t_m  = _bracket(M,   M_LEVELS)
    v_lo, v_hi, t_v  = _bracket(V,   V_LEVELS)
    def ab_at(ic, mm, vv):
        tab_ic = TABLE.get(ic, {})
        tab_m  = tab_ic.get(mm, {})
        if vv in tab_m: return tab_m[vv]
        if not tab_m: raise ValueError('Missing table data')
        vv2, ab = min(tab_m.items(), key=lambda kv: abs(kv[0]-vv))
        return ab
    def lerp(a,b,t): return a+(b-a)*t
    a_ll,b_ll = ab_at(ic_lo,m_lo,v_lo); a_lh,b_lh = ab_at(ic_lo,m_lo,v_hi)
    a_hl,b_hl = ab_at(ic_lo,m_hi,v_lo); a_hh,b_hh = ab_at(ic_lo,m_hi,v_hi)
    a_lm,b_lm = lerp(a_ll,a_lh,t_v), lerp(b_ll,b_lh,t_v)
    a_hm,b_hm = lerp(a_hl,a_hh,t_v), lerp(b_hl,b_hh,t_v)
    a_im,b_im = lerp(a_lm,a_hm,t_m), lerp(b_lm,b_hm,t_m)
    a_ll,b_ll = ab_at(ic_hi,m_lo,v_lo); a_lh,b_lh = ab_at(ic_hi,m_lo,v_hi)
    a_hl,b_hl = ab_at(ic_hi,m_hi,v_lo); a_hh,b_hh = ab_at(ic_hi,m_hi,v_hi)
    a_lm2,b_lm2 = lerp(a_ll,a_lh,t_v), lerp(b_ll,b_lh,t_v)
    a_hm2,b_hm2 = lerp(a_hl,a_hh,t_v), lerp(b_hl,b_hh,t_v)
    a_im2,b_im2 = lerp(a_lm2,a_hm2,t_m), lerp(b_lm2,b_hm2,t_m)
    a = lerp(a_im,a_im2,t_ic); b = lerp(b_im,b_im2,t_ic)
    return a,b

# core

def compute_factors(Vj, Icl, Tmr):
    hc = 8.3*(max(1e-6,Vj)**0.6)
    hr = 3.87 + 0.031*Tmr
    h  = hc + hr
    fcl = 1.0 + 0.2*Icl
    Fcl  = 1.0/(1.0 + 0.155*fcl*h*Icl)
    Fpcl = 1.0/(1.0 + 0.143*hc*Icl)
    return hc,hr,h,fcl,Fcl,Fpcl

def ta50_from_regression(Tmr,Icl,M,Vj):
    a,b = regress_Ta50_coeffs(Icl,M,Vj)
    return a*Tmr + b

def acceptable_line(Icl, M, Vj, Tmr):
    # ASHRAE TABLE 1: Computed m and C Values from Equations 5-7
    # Direct interpolation from ASHRAE Table 1 for higher precision
    vj_table = [0.5, 1.0, 1.5, 2.0]
    m_table = [0.709, 0.722, 0.731, 0.737]
    c_table = [43.03, 48.66, 51.60, 52.85]
    
    # Linear interpolation
    if Vj <= vj_table[0]:
        m = m_table[0]
        C = c_table[0]
    elif Vj >= vj_table[-1]:
        m = m_table[-1]
        C = c_table[-1]
    else:
        for i in range(len(vj_table)-1):
            if vj_table[i] <= Vj <= vj_table[i+1]:
                w = (Vj - vj_table[i]) / (vj_table[i+1] - vj_table[i])
                m = m_table[i] + w * (m_table[i+1] - m_table[i])
                C = c_table[i] + w * (c_table[i+1] - c_table[i])
                break
    
    Ta50 = ta50_from_regression(Tmr,Icl,M,Vj)
    return m,C,Ta50

def jet_ratios(X0,D0, include_buoyancy, TA, TO, V0_guess=10.0):
    r = (X0/D0) + 2.572
    if not include_buoyancy:
        return 1.464/r, 4.539/r
    beta=1.0/273.15; g=9.81; dT=max(0.0,TA-TO)
    Ar=g*beta*dT*D0/max(1e-6,V0_guess**2)
    Vratio=(1.464/r)*((1.0+0.21*Ar*r*r)**(1.0/3.0))
    Tratio=(4.539/r)*((1.0+0.21*Ar*r*r)**(1.0/3.0))
    return Vratio,Tratio

def solve_case(TA,RH_A,Tmr,Vj,M,Icl,D0,X0, rh0=0.95, p_atm_kpa=101.325, include_buoyancy=False):
    PA = pv_from_rh_T(RH_A,TA)
    m,C,Ta50 = acceptable_line(Icl,M,Vj,Tmr)
    # ASHRAE Equations 25 and 27: Find T0 where psychrometric P0 intersects physiological P0
    # Eq 25: P0 = rh0 * Exp[18.6686 - 4030.183/(T0+235)]
    # Eq 27: P0 = -m*T0 + 45.32  (physiological relationship at nozzle)
    # Note: Eq 26 uses C=52.85 for target area, Eq 27 uses 45.32 for nozzle
    C_nozzle = 45.32  # From Eq 27
    
    def residual_T0(T0):
        P0_psychro = rh0*psat_mmhg(T0)           # Eq 25: psychrometric P0
        P0_physio = -m*T0 + C_nozzle             # Eq 27: physiological P0
        return P0_psychro - P0_physio
    
    lo,hi=-5.0,40.0
    def f(x): return residual_T0(x)
    for (a,b) in [(-5.0,40.0),(-20.0,40.0),(-10.0,50.0),(0.0,60.0)]:
        if f(a)*f(b)<=0: lo,hi=a,b; break
    for _ in range(100):
        mid=0.5*(lo+hi); fm=f(mid)
        if abs(fm)<1e-4 or (hi-lo)<1e-4: T0=mid; break
        if f(lo)*fm<=0: hi=mid
        else: lo=mid
    else:
        T0=mid
    P0 = rh0*psat_mmhg(T0)

    # Jet velocity ratio from spread model
    Vratio, Tratio_jet = jet_ratios(X0,D0,include_buoyancy,TA,T0)

    # Target state from intersection of operating line (A->0) with chosen acceptable line (Vj)
    # This keeps report values consistent with charted ASHRAE lines.
    if abs(T0 - TA) < 1e-9:
        Tj = T0
        Pj = P0
    else:
        a_mix = (P0 - PA) / (T0 - TA)
        b_mix = PA - a_mix * TA
        denom = a_mix + m
        if abs(denom) < 1e-9:
            # Fallback to jet-ratio estimate if lines are nearly parallel
            Tj = TA - Tratio_jet*(TA - T0)
            Pj = PA - Tratio_jet*(PA - P0)
        else:
            Tj = (C - b_mix) / denom
            Pj = a_mix * Tj + b_mix

    # Effective thermal ratio consistent with final Tj/T0 used
    if abs(TA - T0) < 1e-9:
        Tratio = Tratio_jet
    else:
        Tratio = (TA - Tj) / (TA - T0)

    RH_j = rh_from_pv_T(Pj,Tj)  # Relative humidity at target area
    V0 = Vj/max(1e-9,Vratio)
    area0 = math.pi*(D0**2)/4.0
    Q0 = V0*area0
    Qj = Q0/max(1e-9,Tratio)
    Qe = max(0.0, Qj-Q0)
    wA = humidity_ratio_from_pv(PA, p_atm_kpa)
    w0 = humidity_ratio_from_pv(P0, p_atm_kpa)
    hA = moist_air_enthalpy_kJkg(TA,wA)
    h0 = moist_air_enthalpy_kJkg(T0,w0)
    rho0 = air_density_approx(T0,w0,p_atm_kpa)
    m_dot0 = rho0*Q0
    Q_total = m_dot0*(hA - h0)
    Q_sens  = m_dot0*1.006*(TA - T0)
    return dict(m=m, C=C, Ta50=Ta50, PA=PA, Pj=Pj, P0=P0, Tj=Tj, T0=T0,
                RH_j=RH_j, RH_0=rh0, Vratio=Vratio, Tratio=Tratio,
                V0=V0, Q0=Q0, Qj=Qj, Qe=Qe, m_dot0=m_dot0, Q_total=Q_total, Q_sens=Q_sens)

def solve_case_from_selected_target(TA, RH_A, Tmr, Vj, M, Icl, D0, X0,
                                    Tj_target, Wj_target,
                                    p_atm_kpa=101.325, include_buoyancy=False):
    PA = pv_from_rh_T(RH_A, TA)
    Pj = pv_from_humidity_ratio(Wj_target, p_atm_kpa)
    m, C, Ta50 = acceptable_line(Icl, M, Vj, Tmr)

    if include_buoyancy:
        T0 = Tj_target
        for _ in range(50):
            Vratio, Tratio = jet_ratios(X0, D0, True, TA, T0)
            T0_new = TA - (TA - Tj_target) / max(1e-9, Tratio)
            if abs(T0_new - T0) < 1e-5:
                T0 = T0_new
                break
            T0 = T0_new
        Vratio, Tratio = jet_ratios(X0, D0, True, TA, T0)
    else:
        Vratio, Tratio = jet_ratios(X0, D0, False, TA, Tj_target)
        T0 = TA - (TA - Tj_target) / max(1e-9, Tratio)

    P0 = PA - (PA - Pj) / max(1e-9, Tratio)
    RH_j = rh_from_pv_T(Pj, Tj_target)
    RH_0 = rh_from_pv_T(P0, T0)

    V0 = Vj / max(1e-9, Vratio)
    area0 = math.pi * (D0**2) / 4.0
    Q0 = V0 * area0
    Qj = Q0 / max(1e-9, Tratio)
    Qe = max(0.0, Qj - Q0)

    wA = humidity_ratio_from_pv(PA, p_atm_kpa)
    w0 = humidity_ratio_from_pv(P0, p_atm_kpa)
    hA = moist_air_enthalpy_kJkg(TA, wA)
    h0 = moist_air_enthalpy_kJkg(T0, w0)
    rho0 = air_density_approx(T0, w0, p_atm_kpa)
    m_dot0 = rho0 * Q0
    Q_total = m_dot0 * (hA - h0)
    Q_sens = m_dot0 * 1.006 * (TA - T0)

    return dict(m=m, C=C, Ta50=Ta50, PA=PA, Pj=Pj, P0=P0, Tj=Tj_target, T0=T0,
                RH_j=RH_j, RH_0=RH_0, Vratio=Vratio, Tratio=Tratio,
                V0=V0, Q0=Q0, Qj=Qj, Qe=Qe, m_dot0=m_dot0, Q_total=Q_total, Q_sens=Q_sens)

# Psychrometric chart helper functions
def get_acceptable_line_points(Vj, T_range=(0, 50), p_atm_kpa=101.325):
    """Generate (T, W) points for acceptable line at given Vj"""
    vj_table = [0.5, 1.0, 1.5, 2.0]
    m_table = [0.709, 0.722, 0.731, 0.737]
    c_table = [43.03, 48.66, 51.60, 52.85]
    
    if Vj <= vj_table[0]:
        m, C = m_table[0], c_table[0]
    elif Vj >= vj_table[-1]:
        m, C = m_table[-1], c_table[-1]
    else:
        for i in range(len(vj_table)-1):
            if vj_table[i] <= Vj <= vj_table[i+1]:
                w = (Vj - vj_table[i]) / (vj_table[i+1] - vj_table[i])
                m = m_table[i] + w * (m_table[i+1] - m_table[i])
                C = c_table[i] + w * (c_table[i+1] - c_table[i])
                break
    
    T_vals = np.linspace(T_range[0], T_range[1], 100)
    P_vals = -m * T_vals + C
    W_vals = [humidity_ratio_from_pv(max(0.1, p), p_atm_kpa) for p in P_vals]
    return T_vals, np.array(W_vals), m, C

def get_physiological_line_points(T_range=(0, 50), rh0=0.95, p_atm_kpa=101.325):
    """Generate (T, W) points for physiological line (Eq 27: P0 = -m*T0 + 45.32)"""
    m = 0.737  # Use nominal m value
    C_physio = 45.32
    T_vals = np.linspace(T_range[0], T_range[1], 100)
    P_vals = -m * T_vals + C_physio
    W_vals = [humidity_ratio_from_pv(max(0.1, p), p_atm_kpa) for p in P_vals]
    return T_vals, np.array(W_vals)

def get_saturation_line_points(T_range=(0, 50), p_atm_kpa=101.325):
    """Generate (T, W) points for saturation curve (RH = 100%)"""
    T_vals = np.linspace(T_range[0], T_range[1], 100)
    try:
        # Psy_TdbRH returns: [Tdb, Twb, RH, W*1000, v, h, Tdp, pw]
        result = Psy_TdbRH(T_vals, 100.0, p_atm_kpa)
        W_vals = result[:, 3] / 1000.0  # Convert from g/kg to kg/kg
        return T_vals, W_vals
    except Exception as e:
        print(f"Warning: Error in saturation line: {e}")
        return T_vals, np.full_like(T_vals, np.nan)

def get_constant_rh_points(RH, T_range=(0, 50), p_atm_kpa=101.325):
    """Generate (T, W) points for constant RH line"""
    T_vals = np.linspace(T_range[0], T_range[1], 100)
    try:
        # Psy_TdbRH returns: [Tdb, Twb, RH, W*1000, v, h, Tdp, pw]
        result = Psy_TdbRH(T_vals, RH, p_atm_kpa)
        W_vals = result[:, 3] / 1000.0  # Convert from g/kg to kg/kg
        return T_vals, W_vals
    except Exception as e:
        print(f"Warning: Error in RH={RH}% line: {e}")
        return T_vals, np.full_like(T_vals, np.nan)

# Matplotlib chart widget
class PsychroChart(FigureCanvas):
    def __init__(self, parent=None, Icl=0.6, M=87, Tmr=45, p_atm_kpa=101.325):
        self.fig = Figure(figsize=(8, 6), dpi=90, tight_layout=True)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        self.parent_window = parent
        self.Icl = Icl
        self.M = M
        self.Tmr = Tmr
        self.p_atm_kpa = p_atm_kpa
        self.selected_Tj = None
        self.selected_Tj_point = None  # Store (T_j, W_j) point
        self.selection_artists = []
        self._sat_T = None
        self._sat_W_g = None
        self._curve_cache = {}
        self.operating_solution = None
        self._layout_initialized = False
        self.mpl_connect('button_press_event', self.on_click)
        self.update_chart()

    def clear_operating_solution(self):
        self.operating_solution = None

    def set_operating_solution(self, res, vj_selected, source='run'):
        self.operating_solution = {
            'source': source,
            'vj_selected': vj_selected,
            'TA': float(self.parent_window.e_TA.value()) if self.parent_window is not None else None,
            'PA': res['PA'],
            'T0': res['T0'],
            'P0': res['P0'],
            'Tj': res['Tj'],
            'Pj': res['Pj'],
        }

    def _get_cached_psychro_curves(self):
        key = round(float(self.p_atm_kpa), 3)
        if key in self._curve_cache:
            return self._curve_cache[key]

        sat = None
        try:
            T_sat, W_sat = get_saturation_line_points(p_atm_kpa=self.p_atm_kpa)
            sat = (T_sat, W_sat)
        except Exception:
            sat = None

        rh_curves = {}
        for rh in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
            try:
                T_rh, W_rh = get_constant_rh_points(rh, p_atm_kpa=self.p_atm_kpa)
                rh_curves[rh] = (T_rh, W_rh)
            except Exception:
                rh_curves[rh] = None

        payload = {'sat': sat, 'rh': rh_curves}
        self._curve_cache[key] = payload
        if len(self._curve_cache) > 8:
            first_key = next(iter(self._curve_cache))
            self._curve_cache.pop(first_key, None)
        return payload

    def _plot_operating_line_and_intersections(self):
        op = self.operating_solution
        if op is None:
            return

        TA = op['TA']
        PA = op['PA']
        T0 = op['T0']
        P0 = op['P0']
        source = op.get('source', 'run')
        if source == 'selected':
            line_color = '#17a2b8'
            line_style = '--'
            line_label = 'Operating line A→0 (Selected point)'
        else:
            line_color = '#6a3d9a'
            line_style = '-.'
            line_label = 'Operating line A→0 (Run)'

        WA_g = humidity_ratio_from_pv(PA, self.p_atm_kpa) * 1000.0
        W0_g = humidity_ratio_from_pv(P0, self.p_atm_kpa) * 1000.0

        self.ax.plot([TA, T0], [WA_g, W0_g], color=line_color, linewidth=2.2,
                     linestyle=line_style, zorder=4, alpha=0.95, label=line_label)
        self.ax.plot(TA, WA_g, marker='s', color=line_color, markersize=7, zorder=5)
        self.ax.text(TA + 0.4, WA_g + 0.35, 'A', color=line_color, fontsize=10, weight='bold', zorder=6)
        self.ax.plot(T0, W0_g, marker='D', color=line_color, markersize=7, zorder=5)
        self.ax.text(T0 + 0.4, W0_g + 0.35, '0', color=line_color, fontsize=10, weight='bold', zorder=6)

        a_mix = (P0 - PA) / (T0 - TA)
        b_mix = PA - a_mix * TA

        if op.get('source') == 'run':
            vj = op.get('vj_selected', 1.0)
            m, C, _ = acceptable_line(self.Icl, self.M, vj, self.Tmr)
            denom = (a_mix + m)
            if abs(denom) >= 1e-9:
                T_int = (C - b_mix) / denom
                P_int = a_mix * T_int + b_mix
                W_int_g = humidity_ratio_from_pv(max(0.1, P_int), self.p_atm_kpa) * 1000.0
                if 0.0 <= T_int <= 50.0 and 0.0 <= W_int_g <= 36.0:
                    self.ax.plot(T_int, W_int_g, marker='o', color='black', markersize=6, zorder=6)
                    self.ax.text(T_int + 0.35, W_int_g - 0.55, f'Int. {vj:.1f} m/s', color='black', fontsize=8, zorder=6)
        else:
            Tj = op.get('Tj')
            Pj = op.get('Pj')
            if Tj is not None and Pj is not None:
                Wj_g = humidity_ratio_from_pv(max(0.1, Pj), self.p_atm_kpa) * 1000.0
                if 0.0 <= Tj <= 50.0 and 0.0 <= Wj_g <= 36.0:
                    self.ax.plot(Tj, Wj_g, marker='o', color='black', markersize=6, zorder=6)
                    self.ax.text(Tj + 0.35, Wj_g - 0.55, 'J seleccionado', color='black', fontsize=8, zorder=6)

    def _clear_selection_artists(self):
        for artist in self.selection_artists:
            try:
                artist.remove()
            except Exception:
                pass
        self.selection_artists = []

    def _draw_selection_overlay(self):
        self._clear_selection_artists()
        if self.selected_Tj is None:
            return
        W_selected_g = self.selected_Tj[1] * 1000  # Convert to g/kg
        marker, = self.ax.plot(self.selected_Tj[0], W_selected_g, 'o', color='red',
                               markersize=12, markeredgecolor='darkred', markeredgewidth=2.5,
                               zorder=5)
        vline = self.ax.axvline(x=self.selected_Tj[0], color='red', linestyle='--', alpha=0.4, linewidth=1, zorder=1)
        hline = self.ax.axhline(y=W_selected_g, color='red', linestyle='--', alpha=0.4, linewidth=1, zorder=1)
        self.selection_artists = [marker, vline, hline]

    def clear_selected_overlay(self):
        self.selected_Tj = None
        self.selected_Tj_point = None
        self._clear_selection_artists()
        self.fig.canvas.draw_idle()
        
    def update_chart(self):
        self.ax.clear()
        curves = self._get_cached_psychro_curves()
        
        # Plot saturation line (RH = 100%) - no label to keep chart clean
        try:
            sat_curve = curves.get('sat')
            if sat_curve is None:
                raise ValueError('Missing saturation curve')
            T_sat, W_sat = sat_curve
            self._sat_T = T_sat
            self._sat_W_g = W_sat * 1000.0
            self.ax.plot(T_sat, W_sat*1000, 'k-', linewidth=1.6, alpha=0.85, zorder=2)
        except Exception as e:
            print(f"Warning: Could not plot saturation line: {e}")
            self._sat_T = None
            self._sat_W_g = None
        
        # Plot constant RH lines (10%, 20%, 30%,...90%) - no label to keep chart clean
        rh_values = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        colors_rh = ['#f0f0f0', '#e8e8e8', '#dcdcdc', '#d0d0d0', '#c0c0c0', '#b0b0b0', '#a0a0a0', '#909090', '#808080']
        
        for rh, color in zip(rh_values, colors_rh):
            try:
                rh_curve = curves['rh'].get(rh)
                if rh_curve is None:
                    raise ValueError(f'Missing RH={rh}% curve')
                T_rh, W_rh = rh_curve
                self.ax.plot(T_rh, W_rh*1000, color=color, linewidth=1.3, linestyle=':',
                            zorder=1, alpha=0.75)
            except Exception as e:
                print(f"Warning: Could not plot RH={rh}% line: {e}")
        
        # Plot acceptable lines for each Vj
        colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']  # Blue, green, orange, red
        vj_list = [0.5, 1.0, 1.5, 2.0]
        line_styles = ['-', '-', '-', '-']
        line_widths = [2.0, 2.2, 2.4, 2.6]
        
        for vj, color, ls, lw in zip(vj_list, colors, line_styles, line_widths):
            T_vals, W_vals, m, C = get_acceptable_line_points(vj, p_atm_kpa=self.p_atm_kpa)
            self.ax.plot(T_vals, W_vals*1000, color=color, linewidth=lw, linestyle=ls, 
                        label=f'V_j = {vj} m/s', zorder=3, alpha=0.85)
            
            # Add rotated text label on left side of saturation curve
            # Find point at T=12°C for positioning
            T_label = 12.0
            # Interpolate W at T_label
            W_at_label = -m * T_label + C
            W_at_label_g = W_at_label * 1000  # Convert to g/kg for display
            
            # Calculate rotation angle based on line slope
            # For the psychrometric line: dW/dT = -m (in kg/kg units)
            # In g/kg units: dW_g/dT ≈ -m * 1000
            # Calculate angle in data coordinates
            T_ref1, T_ref2 = 12.0, 20.0
            W_ref1 = (-m * T_ref1 + C) * 1000
            W_ref2 = (-m * T_ref2 + C) * 1000
            dT = T_ref2 - T_ref1
            dW = W_ref2 - W_ref1
            # Convert to display angle (accounting for axis aspect ratio)
            display_dT = dT / 50.0  # Normalize by xlim
            display_dW = dW / 36.0  # Normalize by ylim (in g/kg)
            angle = math.degrees(math.atan2(display_dW, display_dT))
            
            self.ax.text(T_label, W_at_label_g, f'{vj} m/s', 
                        fontsize=10, color=color, weight='bold', 
                        ha='center', va='center', rotation=angle, zorder=4,
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none'))

        # Plot operating line A->0 and predicted intersections with velocity lines
        self._plot_operating_line_and_intersections()
        
        # Add reference temperature lines (every 5°C)
        for T_ref in range(0, 65, 5):
            self.ax.axvline(x=T_ref, color='gray', linestyle=':', alpha=0.2, linewidth=0.8, zorder=0)
        
        # Formatting
        self.ax.set_xlabel('Dry Bulb Temperature (°C)', fontsize=12, weight='bold', labelpad=10)
        self.ax.set_ylabel('Humidity Ratio (g_water/kg_dry_air)', fontsize=12, weight='bold', labelpad=10)
        self.ax.set_title('Psychrometric Chart - Acceptable Comfort Zones at Different Air Velocities\n(ASHRAE Standard Calculation)', 
                         fontsize=13, weight='bold', pad=15)
        
        # Enhanced grid
        self.ax.grid(True, which='major', alpha=0.5, linestyle='-', linewidth=0.9)
        self.ax.grid(True, which='minor', alpha=0.2, linestyle=':', linewidth=0.5)
        self.ax.minorticks_on()
        
        # Legend (only show Vj lines, physiological, and selected point) - positioned to not block chart
        self.ax.legend(loc='upper right', fontsize=10, framealpha=0.95, edgecolor='black', fancybox=True)
        
        # Set axis limits and ticks
        self.ax.set_xlim(0, 50)  # 0 to 50°C
        self.ax.set_ylim(0, 36)  # 0 to 36 g/kg
        
        self.ax.set_xticks(range(0, 51, 5))
        self.ax.set_xticks(range(0, 51, 1), minor=True)
        
        # Enforce strict axis limits
        self.ax.margins(x=0, y=0)
        
        # Add background color - light background
        self.ax.set_facecolor('#ffffff')
        self.fig.patch.set_facecolor('white')
        
        # Tight layout is expensive; run once to avoid first-interaction stalls
        if not self._layout_initialized:
            self.fig.tight_layout()
            self._layout_initialized = True

        # Draw selected point overlay (if any) without affecting static chart content
        self._draw_selection_overlay()
        
        self.fig.canvas.draw_idle()

    def _validate_selected_point(self, T_sel, W_sel_g):
        # Chart domain limits
        if T_sel < 0.0 or T_sel > 50.0 or W_sel_g < 0.0 or W_sel_g > 36.0:
            return False, 'Point outside chart limits.'

        # Must be to the right/below saturation curve (i.e., not supersaturated)
        if self._sat_T is not None and self._sat_W_g is not None:
            W_sat_g = float(np.interp(T_sel, self._sat_T, self._sat_W_g))
            if W_sel_g > W_sat_g + 0.15:
                return False, 'Invalid point: left of saturation curve.'

        # Must remain within the envelope of expected acceptable lines
        vj_list = [0.5, 1.0, 1.5, 2.0]
        w_lines_g = []
        for vj in vj_list:
            m, C, _ = acceptable_line(self.Icl, self.M, vj, self.Tmr)
            P_line = -m * T_sel + C
            W_line = humidity_ratio_from_pv(max(0.1, P_line), self.p_atm_kpa)
            w_lines_g.append(W_line * 1000.0)

        w_min = min(w_lines_g)
        w_max = max(w_lines_g)
        if W_sel_g < (w_min - 0.35) or W_sel_g > (w_max + 0.35):
            return False, 'Invalid point: outside allowable line band.'

        return True, None
    
    def on_click(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        is_valid, msg = self._validate_selected_point(event.xdata, event.ydata)
        if not is_valid:
            if hasattr(self.parent_window, 'on_invalid_chart_selection'):
                self.parent_window.on_invalid_chart_selection(msg)
            return
        # event.ydata is in g/kg (chart units). Convert to kg/kg for internal storage.
        W_j_g_per_kg = event.ydata
        W_j_kg_per_kg = W_j_g_per_kg / 1000.0 if W_j_g_per_kg is not None else None
        self.selected_Tj = (event.xdata, W_j_kg_per_kg)
        self._draw_selection_overlay()
        self.fig.canvas.draw_idle()
        # Signal parent to update values
        if hasattr(self.parent_window, 'on_chart_selected'):
            # Pass W_j in kg/kg to the parent for consistent units
            self.parent_window.on_chart_selected(event.xdata, W_j_kg_per_kg)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Spot Cooling - ASHRAE (Azer)')
        self.resize(1600, 900)
        self.last_report_html = None  # Store report for PDF export
        cw=QtWidgets.QWidget(); self.setCentralWidget(cw)
        layout=QtWidgets.QHBoxLayout(cw)
        
        # Left panel: ONLY form controls
        form=QtWidgets.QFormLayout()
        form_widget=QtWidgets.QWidget(); form_widget.setLayout(form)
        scroll=QtWidgets.QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(form_widget)
        left_panel = QtWidgets.QWidget()
        left_panel.setLayout(QtWidgets.QVBoxLayout())
        left_panel.layout().addWidget(scroll)
        
        # Right panel: FULL psychrometric chart
        self.chart = PsychroChart(parent=self)
        
        # Use QSplitter to allow resizing between left and right panels
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(left_panel)
        splitter.addWidget(self.chart)
        splitter.setStretchFactor(0, 1)  # Left panel has some stretch
        splitter.setStretchFactor(1, 2)  # Right panel has more stretch (chart)
        splitter.setSizes([350, 1000])   # Initial sizes: 350px for controls, ~1000px for chart
        layout.addWidget(splitter)
        
        # inputs
        self.e_TA=QtWidgets.QDoubleSpinBox(); self.e_TA.setRange(-20,80); self.e_TA.setValue(40.0); self.e_TA.setSuffix(' C')
        self.e_RHA=QtWidgets.QDoubleSpinBox(); self.e_RHA.setRange(0,100); self.e_RHA.setValue(50.0); self.e_RHA.setSuffix(' %')
        self.e_TMR=QtWidgets.QDoubleSpinBox(); self.e_TMR.setRange(-20,120); self.e_TMR.setValue(45.0); self.e_TMR.setSuffix(' C')
        self.e_VJ=QtWidgets.QDoubleSpinBox(); self.e_VJ.setRange(0.2,3.0); self.e_VJ.setSingleStep(0.1); self.e_VJ.setValue(2.0); self.e_VJ.setSuffix(' m/s')
        self.e_M=QtWidgets.QDoubleSpinBox(); self.e_M.setRange(40,400); self.e_M.setValue(87.0); self.e_M.setSuffix(' W/m2')
        self.e_ICL=QtWidgets.QDoubleSpinBox(); self.e_ICL.setRange(0.0,2.0); self.e_ICL.setSingleStep(0.05); self.e_ICL.setValue(0.6); self.e_ICL.setSuffix(' clo')
        self.e_D0=QtWidgets.QDoubleSpinBox(); self.e_D0.setRange(0.03,0.5); self.e_D0.setSingleStep(0.005); self.e_D0.setValue(0.3048); self.e_D0.setSuffix(' m')
        self.e_X0=QtWidgets.QDoubleSpinBox(); self.e_X0.setRange(0.3,4.0); self.e_X0.setSingleStep(0.01); self.e_X0.setValue(3.048); self.e_X0.setSuffix(' m')  # 10 ft (ASHRAE Example 1)
        self.e_RT=QtWidgets.QDoubleSpinBox(); self.e_RT.setRange(0.1,1.0); self.e_RT.setSingleStep(0.01); self.e_RT.setValue(0.3048); self.e_RT.setSuffix(' m')
        self.e_ang=QtWidgets.QDoubleSpinBox(); self.e_ang.setRange(5,45); self.e_ang.setValue(22.0); self.e_ang.setSuffix(' deg')
        self.btn_calcX0=QtWidgets.QPushButton('Calc X0 from Rt & angle')
        self.e_RH0=QtWidgets.QDoubleSpinBox(); self.e_RH0.setRange(80,100); self.e_RH0.setValue(95.0); self.e_RH0.setSuffix(' %')
        self.e_T0_inv=QtWidgets.QDoubleSpinBox(); self.e_T0_inv.setRange(-50,80); self.e_T0_inv.setDecimals(2); self.e_T0_inv.setReadOnly(True); self.e_T0_inv.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons); self.e_T0_inv.setSuffix(' C')
        self.e_P0_inv=QtWidgets.QDoubleSpinBox(); self.e_P0_inv.setRange(0,100); self.e_P0_inv.setDecimals(3); self.e_P0_inv.setReadOnly(True); self.e_P0_inv.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons); self.e_P0_inv.setSuffix(' mmHg')
        self.e_V0_inv=QtWidgets.QDoubleSpinBox(); self.e_V0_inv.setRange(0,100); self.e_V0_inv.setDecimals(3); self.e_V0_inv.setReadOnly(True); self.e_V0_inv.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons); self.e_V0_inv.setSuffix(' m/s')
        self.e_PATM=QtWidgets.QDoubleSpinBox(); self.e_PATM.setRange(70,110); self.e_PATM.setValue(101.325); self.e_PATM.setDecimals(3); self.e_PATM.setSuffix(' kPa')
        self.chk_buoy=QtWidgets.QCheckBox('Include buoyancy (Ar)')
        self.btn_run=QtWidgets.QPushButton('Run & report')
        self.btn_use_selected=QtWidgets.QPushButton('Use selected point')
        self.btn_use_selected.setEnabled(False)
        self.btn_pdf=QtWidgets.QPushButton('Export PDF...')
        self.lbl_selected=QtWidgets.QLabel('Selected from chart: None')

        # Reduce redraw storms while typing values in spinboxes
        self._chart_update_timer = QtCore.QTimer(self)
        self._chart_update_timer.setSingleShot(True)
        self._chart_update_timer.setInterval(220)
        self._chart_update_timer.timeout.connect(self._apply_chart_update)

        spinboxes = [
            self.e_TA, self.e_RHA, self.e_TMR, self.e_VJ, self.e_M, self.e_ICL,
            self.e_D0, self.e_X0, self.e_RT, self.e_ang, self.e_RH0,
            self.e_T0_inv, self.e_P0_inv, self.e_V0_inv, self.e_PATM
        ]
        for spin in spinboxes:
            spin.setKeyboardTracking(False)
        
        # layout
        form.addRow(QtWidgets.QLabel('<b>Ambient</b>'))
        form.addRow('T_A:', self.e_TA); form.addRow('RH_A:', self.e_RHA); form.addRow('T_mr:', self.e_TMR); form.addRow('V_j:', self.e_VJ)
        form.addRow(QtWidgets.QLabel('<b>Work & clothing</b>'))
        form.addRow('M:', self.e_M); form.addRow('I_cl:', self.e_ICL)
        form.addRow(QtWidgets.QLabel('<b>Jet geometry</b>'))
        form.addRow('D_0:', self.e_D0); form.addRow('X_0:', self.e_X0); form.addRow('R_t:', self.e_RT); form.addRow('Angle:', self.e_ang); form.addRow(self.btn_calcX0)
        form.addRow(QtWidgets.QLabel('<b>Coil/nozzle</b>'))
        form.addRow('rh_0:', self.e_RH0)
        form.addRow(QtWidgets.QLabel('<b>From selected point (inverse)</b>'))
        form.addRow('T_0*:', self.e_T0_inv)
        form.addRow('P_0*:', self.e_P0_inv)
        form.addRow('V_0*:', self.e_V0_inv)
        form.addRow(QtWidgets.QLabel('<b>Other</b>'))
        form.addRow('P_atm:', self.e_PATM); form.addRow(self.chk_buoy); form.addRow(self.btn_run); form.addRow(self.btn_use_selected); form.addRow(self.btn_pdf)
        form.addRow(QtWidgets.QLabel('<b>Psychrometric</b>'))
        form.addRow(self.lbl_selected)
        
        # signals
        self.btn_calcX0.clicked.connect(self.calc_x0_from_rt)
        self.btn_run.clicked.connect(self.run_calc)
        self.btn_use_selected.clicked.connect(self.use_selected_point)
        self.btn_pdf.clicked.connect(self.export_pdf)
        self.e_TA.valueChanged.connect(self.on_chart_params_changed)
        self.e_RHA.valueChanged.connect(self.on_chart_params_changed)
        self.e_VJ.valueChanged.connect(self.on_chart_params_changed)
        self.e_RH0.valueChanged.connect(self.on_chart_params_changed)
        self.chk_buoy.toggled.connect(self.on_chart_params_changed)
        self.e_TMR.valueChanged.connect(self.on_chart_params_changed)
        self.e_M.valueChanged.connect(self.on_chart_params_changed)
        self.e_ICL.valueChanged.connect(self.on_chart_params_changed)
        self.e_PATM.valueChanged.connect(self.on_chart_params_changed)

    def _is_physical_solution(self, r, TA=None):
        vals = [r.get('T0'), r.get('Tj'), r.get('P0'), r.get('Pj'), r.get('RH_0'), r.get('RH_j')]
        if any((v is None or not np.isfinite(v)) for v in vals):
            return False
        if r['P0'] <= 0.0 or r['Pj'] <= 0.0:
            return False
        if not (0.0 <= r['RH_0'] <= 1.0 and 0.0 <= r['RH_j'] <= 1.0):
            return False
        # Thermal consistency: for cooling jet, Tj cannot be below T0.
        # Generalized as Tj must lie between TA and T0.
        if TA is not None and np.isfinite(TA):
            t_min = min(TA, r['T0']) - 1e-6
            t_max = max(TA, r['T0']) + 1e-6
            if not (t_min <= r['Tj'] <= t_max):
                return False
        elif r['Tj'] < r['T0'] - 1e-6:
            return False
        return True

    def calc_x0_from_rt(self):
        Rt=self.e_RT.value(); D0=self.e_D0.value(); ang=self.e_ang.value()*math.pi/180.0
        a=0.5/D0; b=math.tan(ang); c=-Rt
        disc=b*b-4*a*c
        if disc<=0:
            QtWidgets.QMessageBox.warning(self,'Geometry','No real solution for X0.')
            return
        x1=(-b+disc**0.5)/(2*a); x2=(-b-disc**0.5)/(2*a)
        X0=max(x1,x2); self.e_X0.setValue(X0)

    def on_chart_params_changed(self):
        """Update chart when input parameters change"""
        if hasattr(self, '_chart_update_timer'):
            self._chart_update_timer.start()

    def _apply_chart_update(self):
        if hasattr(self, 'chart'):
            self.chart.Icl = self.e_ICL.value()
            self.chart.M = self.e_M.value()
            self.chart.Tmr = self.e_TMR.value()
            self.chart.p_atm_kpa = self.e_PATM.value()
            self.chart.clear_operating_solution()
            self.chart.update_chart()
    
    def on_chart_selected(self, T_j, W_j):
        """Handle selection from psychrometric chart"""
        # Here W_j is in kg/kg; store it and show in g/kg for readability
        self.chart.selected_Tj_point = (T_j, W_j)  # Store for use_selected_point
        self.lbl_selected.setText(f'Selected: T_j={T_j:.1f}°C, W_j={W_j*1000:.2f} g/kg')
        self.btn_use_selected.setEnabled(True)

    def on_invalid_chart_selection(self, msg):
        """Handle invalid point selection from chart"""
        self.chart.clear_selected_overlay()
        self.lbl_selected.setText(f'Selected from chart: None ({msg})')
        self.btn_use_selected.setEnabled(False)
        QtWidgets.QMessageBox.warning(self, 'Invalid selection', msg)

    def use_selected_point(self):
        """Calculate with the selected point from the chart"""
        if not hasattr(self.chart, 'selected_Tj_point') or self.chart.selected_Tj_point is None:
            QtWidgets.QMessageBox.warning(self, 'Error', 'No point selected on chart')
            return
        
        T_j_selected, W_j_selected = self.chart.selected_Tj_point
        
        try:
            # Run inverse calculation constrained by selected target point (T_j, W_j)
            res = solve_case_from_selected_target(
                self.e_TA.value(), self.e_RHA.value()/100.0, self.e_TMR.value(),
                self.e_VJ.value(), self.e_M.value(), self.e_ICL.value(),
                self.e_D0.value(), self.e_X0.value(),
                T_j_selected, W_j_selected,
                p_atm_kpa=self.e_PATM.value(),
                include_buoyancy=self.chk_buoy.isChecked()
            )

            if not self._is_physical_solution(res, TA=self.e_TA.value()):
                self.chart.clear_operating_solution()
                self.chart.update_chart()
                QtWidgets.QMessageBox.warning(
                    self,
                    'Invalid selection',
                    'Selected point does not produce a physically feasible solution (e.g., T_j outside TA↔T_0 range). Report was not generated.'
                )
                return

            # Sync UI controls with inverse solution so the state is visible and reusable
            self.e_RH0.setValue(max(self.e_RH0.minimum(), min(self.e_RH0.maximum(), res['RH_0'] * 100.0)))
            self.e_T0_inv.setValue(res['T0'])
            self.e_P0_inv.setValue(res['P0'])
            self.e_V0_inv.setValue(res['V0'])

            # Draw operating line A->0 for selected-point inverse solution
            self.chart.set_operating_solution(res, self.e_VJ.value(), source='selected')
            self.chart.update_chart()

            
            # Show results in dialog
            html = self.build_html(res)
            self.last_report_html = html
            self.show_report_dialog(html)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, 'Error', str(e))

    def run_calc(self):
        try:
            # Run mode should not keep/display previously selected user point
            self.chart.clear_selected_overlay()
            self.lbl_selected.setText('Selected from chart: None')
            self.btn_use_selected.setEnabled(False)

            res = solve_case(self.e_TA.value(), self.e_RHA.value()/100.0, self.e_TMR.value(), self.e_VJ.value(), self.e_M.value(), self.e_ICL.value(), self.e_D0.value(), self.e_X0.value(), rh0=self.e_RH0.value()/100.0, p_atm_kpa=self.e_PATM.value(), include_buoyancy=self.chk_buoy.isChecked())
            if not self._is_physical_solution(res, TA=self.e_TA.value()):
                self.chart.clear_operating_solution()
                self.chart.update_chart()
                QtWidgets.QMessageBox.warning(self, 'Invalid solution', 'Run result is not physically feasible (e.g., T_j outside TA↔T_0 range). Report was not generated.')
                return

            # Draw operating line and intersection with selected velocity line for this run
            self.chart.set_operating_solution(res, self.e_VJ.value(), source='run')
            self.chart.update_chart()

            html = self.build_html(res)
            self.last_report_html = html
            self.show_report_dialog(html)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self,'Error',str(e))

    def build_html(self, r):
        def fmt(x,d=3):
            try: return f"{x:.{d}f}"
            except: return str(x)
        html = """
        <html><head><meta charset='utf-8'>
        <style>body{font-family:Arial;margin:20px}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccc;padding:6px;text-align:right}th{background:#f2f2f2}.left{text-align:left}</style>
        </head><body>
        <h1>Spot Cooling Report</h1>
        <h2>Results</h2>
        <table>
        <tr><th class='left'>Var</th><th>Value</th><th>Units</th></tr>
        """
        rows = [
            ('V_j (used)', fmt(self.e_VJ.value(),3), 'm/s'),
            ('T_0', fmt(r['T0'],3), 'C'),
            ('P_0', fmt(r['P0'],3), 'mmHg'),
            ('RH_0', fmt(r['RH_0']*100,2), '%'),
            ('T_j', fmt(r['Tj'],3), 'C'),
            ('P_j', fmt(r['Pj'],3), 'mmHg'),
            ('RH_j', fmt(r['RH_j']*100,2), '%'),
            ('V_0', fmt(r['V0'],3), 'm/s'),
            ('V_j/V_0', fmt(r['Vratio'],4), '-'),
            ('(T_A-T_j)/(T_A-T_0)', fmt(r['Tratio'],4), '-'),
            ('Q_0', fmt(r['Q0'],4), 'm3/s'),
            ('Q_j', fmt(r['Qj'],4), 'm3/s'),
            ('Entrainment', fmt(r['Qe'],4), 'm3/s'),
            ('m_dot0', fmt(r['m_dot0'],4), 'kg/s'),
            ('Q_total', fmt(r['Q_total'],3), 'kW'),
            ('Q_sensible', fmt(r['Q_sens'],3), 'kW'),
        ]
        for k,v,u in rows:
            html += f"<tr><td class='left'>{k}</td><td>{v}</td><td>{u}</td></tr>"
        html += "</table>"
        html += f"<p>m={fmt(r['m'],4)} mmHg/C; C={fmt(r['C'],3)} mmHg; T_a(0.5)={fmt(r['Ta50'],2)} C.</p>"
        html += "</body></html>"
        return html

    def show_report_dialog(self, html):
        """Show report in a popup dialog with export button"""
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle('Spot Cooling Report')
        dlg.resize(700, 600)
        layout = QtWidgets.QVBoxLayout(dlg)
        view = QtWidgets.QTextBrowser()
        view.setHtml(html)
        layout.addWidget(view)
        btn_export = QtWidgets.QPushButton('Export as PDF')
        btn_export.clicked.connect(lambda: self.export_pdf_from_dialog(html))
        layout.addWidget(btn_export)
        dlg.exec_()
    
    def export_pdf_from_dialog(self, html):
        """Export report dialog to PDF"""
        path,_=QtWidgets.QFileDialog.getSaveFileName(self,'Export to PDF','spot_cooling_report.pdf','PDF (*.pdf)')
        if not path: return
        printer=QtPrintSupport.QPrinter(QtPrintSupport.QPrinter.HighResolution)
        printer.setOutputFormat(QtPrintSupport.QPrinter.PdfFormat)
        printer.setOutputFileName(path)
        printer.setPageMargins(12,12,12,12, QtPrintSupport.QPrinter.Millimeter)
        doc=QtGui.QTextDocument(); doc.setHtml(html); doc.print_(printer)

    def export_pdf(self):
        """Legacy method for export_pdf button"""
        if self.last_report_html is None:
            QtWidgets.QMessageBox.warning(self, 'Warning', 'No report generated yet. Run a calculation first.')
            return
        self.export_pdf_from_dialog(self.last_report_html)

if __name__=='__main__':
    app=QtWidgets.QApplication(sys.argv)
    w=MainWindow(); w.show()
    sys.exit(app.exec_())
