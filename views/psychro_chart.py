# -*- coding: utf-8 -*-
"""
Psychrometric chart widget for PyQt5.
Displays acceptable comfort zones and operating lines.
"""

import math
import numpy as np
from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from models.psychrometric import (
    saturation_vapor_pressure_mmhg, humidity_ratio_from_vapor_pressure
)
from models.cooling_calc import acceptable_line_params


class PsychroChart(FigureCanvas):
    """Matplotlib-based psychrometric chart for spot cooling design."""
    
    def __init__(self, parent=None, icl_clo=0.6, metabolic_rate=87,
                 mean_radiant_temp=45, p_atm_kpa=101.325):
        self.fig = Figure(figsize=(8, 6), dpi=90)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        self.parent_window = parent
        
        self.Icl = icl_clo
        self.M = metabolic_rate
        self.Tmr = mean_radiant_temp
        self.p_atm_kpa = p_atm_kpa
        
        self.selected_Tj = None
        self.selected_Tj_point = None
        self.selected_vj = None
        self.selection_artists = []
        self._sat_T = None
        self._sat_W_g = None
        self._curve_cache = {}
        self.operating_solution = None
        self._startup_full_render_done = False
        
        self.mpl_connect('button_press_event', self.on_click)
        
        # Async render: fast first, then full
        QtCore.QTimer.singleShot(0, self._render_startup_chart)
    
    def _render_startup_chart(self):
        """Fast initial render without heavy psychrometric background"""
        self.update_chart(skip_psychro_background=True)
        QtCore.QTimer.singleShot(10, self._render_full_startup_chart)
    
    def _render_full_startup_chart(self):
        """Full render with psychrometric background"""
        self._startup_full_render_done = True
        self.update_chart(skip_psychro_background=False)
    
    def clear_operating_solution(self):
        """Clear operating line display"""
        self.operating_solution = None
    
    def set_operating_solution(self, result, vj_selected, source='run'):
        """Store solution for display on chart"""
        self.operating_solution = {
            'source': source,
            'vj_selected': vj_selected,
            'TA': float(self.parent_window.e_TA.value()) if self.parent_window is not None else None,
            'PA': result['PA'],
            'T0': result['T0'],
            'P0': result['P0'],
            'Tj': result['Tj'],
            'Pj': result['Pj'],
        }
    
    def _get_cached_psychro_curves(self):
        """Cache psychrometric curves by pressure"""
        key = round(float(self.p_atm_kpa), 3)
        if key in self._curve_cache:
            return self._curve_cache[key]
        
        # Saturation curve
        sat = None
        try:
            T_sat = np.linspace(0, 50, 100)
            W_sat = np.array([humidity_ratio_from_vapor_pressure(
                max(0.1, saturation_vapor_pressure_mmhg(float(T)))
            ) for T in T_sat], dtype=float)
            sat = (T_sat, W_sat)
        except Exception:
            sat = None
        
        # Constant RH curves
        rh_curves = {}
        for rh in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
            try:
                T_rh = np.linspace(0, 50, 100)
                rh_frac = float(rh) / 100.0
                pv_vals = np.array([rh_frac * saturation_vapor_pressure_mmhg(float(T))
                                   for T in T_rh], dtype=float)
                W_rh = np.array([humidity_ratio_from_vapor_pressure(max(0.1, pv))
                                for pv in pv_vals], dtype=float)
                rh_curves[rh] = (T_rh, W_rh)
            except Exception:
                rh_curves[rh] = None
        
        payload = {'sat': sat, 'rh': rh_curves}
        self._curve_cache[key] = payload
        
        # Limit cache size
        if len(self._curve_cache) > 8:
            first_key = next(iter(self._curve_cache))
            self._curve_cache.pop(first_key, None)
        
        return payload
    
    def _plot_operating_line_and_intersections(self):
        """Plot operating line from ambient to nozzle condition"""
        if self.operating_solution is None:
            return
        
        op = self.operating_solution
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
        
        WA_g = humidity_ratio_from_vapor_pressure(PA, self.p_atm_kpa) * 1000.0
        W0_g = humidity_ratio_from_vapor_pressure(P0, self.p_atm_kpa) * 1000.0
        
        self.ax.plot([TA, T0], [WA_g, W0_g], color=line_color, linewidth=2.2,
                    linestyle=line_style, zorder=4, alpha=0.95, label=line_label)
        self.ax.plot(TA, WA_g, marker='s', color=line_color, markersize=7, zorder=5)
        self.ax.text(TA + 0.4, WA_g + 0.35, 'A', color=line_color, fontsize=10,
                    weight='bold', zorder=6)
        self.ax.plot(T0, W0_g, marker='D', color=line_color, markersize=7, zorder=5)
        self.ax.text(T0 + 0.4, W0_g + 0.35, '0', color=line_color, fontsize=10,
                    weight='bold', zorder=6)
        
        # Intersection with velocity line
        a_mix = (P0 - PA) / (T0 - TA) if abs(T0 - TA) > 1e-9 else 0
        b_mix = PA - a_mix * TA
        
        if op.get('source') == 'run':
            vj = op.get('vj_selected', 1.0)
            m, C, _ = acceptable_line_params(self.Icl, self.M, vj, self.Tmr)
            denom = (a_mix + m)
            if abs(denom) >= 1e-9:
                T_int = (C - b_mix) / denom
                P_int = a_mix * T_int + b_mix
                W_int_g = humidity_ratio_from_vapor_pressure(max(0.1, P_int),
                                                            self.p_atm_kpa) * 1000.0
                if 0.0 <= T_int <= 50.0 and 0.0 <= W_int_g <= 36.0:
                    self.ax.plot(T_int, W_int_g, marker='o', color='black',
                               markersize=6, zorder=6)
                    self.ax.text(T_int + 0.35, W_int_g - 0.55, f'Int. {vj:.1f} m/s',
                               color='black', fontsize=8, zorder=6)
        else:
            Tj = op.get('Tj')
            Pj = op.get('Pj')
            if Tj is not None and Pj is not None:
                Wj_g = humidity_ratio_from_vapor_pressure(max(0.1, Pj),
                                                         self.p_atm_kpa) * 1000.0
                if 0.0 <= Tj <= 50.0 and 0.0 <= Wj_g <= 36.0:
                    self.ax.plot(Tj, Wj_g, marker='o', color='black',
                               markersize=6, zorder=6)
                    self.ax.text(Tj + 0.35, Wj_g - 0.55, 'J seleccionado',
                               color='black', fontsize=8, zorder=6)
    
    def _clear_selection_artists(self):
        """Remove selection overlay artists"""
        for artist in self.selection_artists:
            try:
                artist.remove()
            except Exception:
                pass
        self.selection_artists = []
    
    def _draw_selection_overlay(self):
        """Draw point selection overlay"""
        self._clear_selection_artists()
        if self.selected_Tj is None:
            return
        
        W_selected_g = self.selected_Tj[1] * 1000  # kg/kg -> g/kg
        marker, = self.ax.plot(self.selected_Tj[0], W_selected_g, 'o', color='red',
                              markersize=12, markeredgecolor='darkred', markeredgewidth=2.5,
                              zorder=5)
        vline = self.ax.axvline(x=self.selected_Tj[0], color='red', linestyle='--',
                               alpha=0.4, linewidth=1, zorder=1)
        hline = self.ax.axhline(y=W_selected_g, color='red', linestyle='--',
                               alpha=0.4, linewidth=1, zorder=1)
        self.selection_artists = [marker, vline, hline]
    
    def clear_selected_overlay(self):
        """Clear point selection overlay"""
        self.selected_Tj = None
        self.selected_Tj_point = None
        self.selected_vj = None
        self._clear_selection_artists()
        self.fig.canvas.draw_idle()

    def _infer_velocity_from_selected_point(self, T_sel, W_sel_g):
        """Infer continuous Vj from selected point using acceptable-line inversion."""
        v_min = 0.5
        v_max = 2.0

        def residual(vj):
            m, C, _ = acceptable_line_params(self.Icl, self.M, float(vj), self.Tmr)
            P_line = -m * T_sel + C
            W_line_g = humidity_ratio_from_vapor_pressure(max(0.1, P_line),
                                                          self.p_atm_kpa) * 1000.0
            return W_line_g - W_sel_g

        # Coarse scan to locate nearest value and any sign-change bracket.
        v_samples = np.linspace(v_min, v_max, 76)
        f_samples = []
        for v in v_samples:
            try:
                f_samples.append(float(residual(v)))
            except Exception:
                f_samples.append(np.nan)

        valid = [(float(v), float(f)) for v, f in zip(v_samples, f_samples) if np.isfinite(f)]
        if not valid:
            return None

        # Default fallback: nearest sampled value.
        best_v, _ = min(valid, key=lambda vf: abs(vf[1]))

        # Refine with bisection when a crossing exists.
        bracket = None
        for i in range(len(valid) - 1):
            v1, f1 = valid[i]
            v2, f2 = valid[i + 1]
            if f1 == 0.0:
                return v1
            if f1 * f2 <= 0.0:
                bracket = (v1, v2)
                break

        if bracket is None:
            return round(best_v, 3)

        lo, hi = bracket
        flo = residual(lo)
        fhi = residual(hi)
        for _ in range(32):
            mid = 0.5 * (lo + hi)
            fmid = residual(mid)
            if abs(fmid) < 1e-6 or (hi - lo) < 1e-4:
                return round(float(mid), 3)
            if flo * fmid <= 0.0:
                hi = mid
                fhi = fmid
            else:
                lo = mid
                flo = fmid

        return round(float(0.5 * (lo + hi)), 3)
    
    def update_chart(self, skip_psychro_background=False):
        """Redraw psychrometric chart"""
        self.ax.clear()
        curves = None if skip_psychro_background else self._get_cached_psychro_curves()
        
        # Saturation curve
        if curves is not None:
            try:
                sat_curve = curves.get('sat')
                if sat_curve is None:
                    raise ValueError('Missing saturation curve')
                T_sat, W_sat = sat_curve
                self._sat_T = T_sat
                self._sat_W_g = W_sat * 1000.0
                self.ax.plot(T_sat, W_sat * 1000, color='#808080', linestyle='-',
                           linewidth=1.5, alpha=0.75, zorder=2)
            except Exception as e:
                print(f"Warning: Could not plot saturation line: {e}")
                self._sat_T = None
                self._sat_W_g = None
        else:
            self._sat_T = None
            self._sat_W_g = None
        
        # Constant RH lines
        rh_values = [10, 20, 30, 40, 50, 60, 70, 80, 90]
        if curves is not None:
            for rh in rh_values:
                try:
                    rh_curve = curves['rh'].get(rh)
                    if rh_curve is None:
                        raise ValueError(f'Missing RH={rh}% curve')
                    T_rh, W_rh = rh_curve
                    self.ax.plot(T_rh, W_rh * 1000, color='#808080', linewidth=1.3,
                               linestyle=':', zorder=1, alpha=0.75)
                except Exception as e:
                    print(f"Warning: Could not plot RH={rh}% line: {e}")
        
        # Acceptable lines for each Vj
        colors = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728']
        vj_list = [0.5, 1.0, 1.5, 2.0]
        line_styles = ['-', '-', '-', '-']
        line_widths = [2.0, 2.2, 2.4, 2.6]
        
        for vj, color, ls, lw in zip(vj_list, colors, line_styles, line_widths):
            try:
                T_vals = np.linspace(0, 50, 100)
                m, C, _ = acceptable_line_params(self.Icl, self.M, vj, self.Tmr)
                P_vals = -m * T_vals + C
                W_vals = np.array([humidity_ratio_from_vapor_pressure(max(0.1, p),
                                                                     self.p_atm_kpa)
                                  for p in P_vals], dtype=float)
                self.ax.plot(T_vals, W_vals * 1000, color=color, linewidth=lw,
                           linestyle=ls, label=f'V_j = {vj} m/s', zorder=3, alpha=0.85)
                
                # Label
                T_label = 12.0
                W_at_label = -m * T_label + C
                W_at_label_g = W_at_label * 1000
                
                T_ref1, T_ref2 = 12.0, 20.0
                W_ref1 = (-m * T_ref1 + C) * 1000
                W_ref2 = (-m * T_ref2 + C) * 1000
                dT = T_ref2 - T_ref1
                dW = W_ref2 - W_ref1
                display_dT = dT / 50.0
                display_dW = dW / 36.0
                angle = math.degrees(math.atan2(display_dW, display_dT))
                
                self.ax.text(T_label, W_at_label_g, f'{vj} m/s', fontsize=10,
                           color=color, weight='bold', ha='center', va='center',
                           rotation=angle, zorder=4,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                    alpha=0.7, edgecolor='none'))
            except Exception as e:
                print(f"Warning: Could not plot Vj={vj} line: {e}")
        
        # Operating line
        self._plot_operating_line_and_intersections()
        
        # Reference lines
        for T_ref in range(0, 65, 5):
            self.ax.axvline(x=T_ref, color='gray', linestyle=':', alpha=0.2,
                           linewidth=0.8, zorder=0)
        
        # Formatting
        self.ax.set_xlabel('Dry Bulb Temperature (°C)', fontsize=12, weight='bold',
                          labelpad=10)
        self.ax.set_ylabel(r'Humidity Ratio ($\mathbf{g_{w}/kg_{da}}$)', fontsize=12,
                          weight='bold', labelpad=10)
        self.ax.set_title('Psychrometric Chart - Acceptable Comfort Zones\n'
                         '(ASHRAE Standard RP-884)',
                         fontsize=13, weight='bold', pad=15)
        
        self.ax.grid(True, which='major', alpha=0.5, linestyle='-', linewidth=0.9)
        self.ax.grid(True, which='minor', alpha=0.2, linestyle=':', linewidth=0.5)
        self.ax.minorticks_on()
        
        self.ax.legend(loc='upper right', fontsize=10, framealpha=0.95,
                      edgecolor='black', fancybox=True)
        
        self.ax.set_xlim(0, 50)
        self.ax.set_ylim(0, 36)
        self.ax.set_xticks(range(0, 51, 5))
        self.ax.set_xticks(range(0, 51, 1), minor=True)
        self.ax.margins(x=0, y=0)
        
        self.ax.set_facecolor('#ffffff')
        self.fig.patch.set_facecolor('white')
        self.fig.subplots_adjust(left=0.10, right=0.98, bottom=0.11, top=0.90)
        
        self._draw_selection_overlay()
        self.fig.canvas.draw_idle()
    
    def _validate_selected_point(self, T_sel, W_sel_g):
        """Validate selected point is inside chart and not beyond saturation."""
        if T_sel < 0.0 or T_sel > 50.0 or W_sel_g < 0.0 or W_sel_g > 36.0:
            return False, 'Point outside chart limits.'
        
        if self._sat_T is not None and self._sat_W_g is not None:
            W_sat_g = float(np.interp(T_sel, self._sat_T, self._sat_W_g))
            if W_sel_g > W_sat_g + 0.15:
                return False, 'Invalid point: left of saturation curve.'
        
        return True, None
    
    def on_click(self, event):
        """Handle chart click for point selection"""
        if event.inaxes != self.ax or event.xdata is None:
            return
        
        is_valid, msg = self._validate_selected_point(event.xdata, event.ydata)
        if not is_valid:
            if hasattr(self.parent_window, 'on_invalid_chart_selection'):
                self.parent_window.on_invalid_chart_selection(msg)
            return
        
        T_j_c = float(event.xdata)
        W_j_g_per_kg = float(event.ydata)
        W_j_kg_per_kg = W_j_g_per_kg / 1000.0 if W_j_g_per_kg is not None else None
        self.selected_vj = self._infer_velocity_from_selected_point(T_j_c, W_j_g_per_kg)
        self.selected_Tj = (T_j_c, W_j_kg_per_kg)
        self._draw_selection_overlay()
        self.fig.canvas.draw_idle()
        
        if hasattr(self.parent_window, 'on_chart_selected'):
            self.parent_window.on_chart_selected(T_j_c, W_j_kg_per_kg, self.selected_vj)
