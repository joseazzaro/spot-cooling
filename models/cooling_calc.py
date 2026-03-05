# -*- coding: utf-8 -*-
"""
Spot cooling design calculations based on ASHRAE Standard RP-884.
Implements the thermal comfort and jet dynamics models.
"""

import math
import numpy as np
from utils.constants import TABLE, ICL_LEVELS, M_LEVELS, V_LEVELS
from utils.helpers import bracket_value, linear_interpolate
from models.psychrometric import (
    saturation_vapor_pressure_mmhg, vapor_pressure_from_rh_and_temp,
    relative_humidity_from_vapor_pressure, humidity_ratio_from_vapor_pressure,
    vapor_pressure_from_humidity_ratio, enthalpy_moist_air, air_density
)


def regress_ta50_coefficients(icl_clo, metabolic_rate, jet_velocity):
    """
    Interpolate regression coefficients (a, b) for Ta(0.5) calculation.
    Ta(0.5) = a*Tmr + b, where a,b depend on Icl, M, Vj
    
    Returns: tuple (a, b)
    """
    # Bracket values in each dimension
    ic_lo, ic_hi, t_ic = bracket_value(icl_clo, ICL_LEVELS)
    m_lo, m_hi, t_m = bracket_value(metabolic_rate, M_LEVELS)
    v_lo, v_hi, t_v = bracket_value(jet_velocity, V_LEVELS)
    
    def get_ab_at(ic, mm, vv):
        """Get (a,b) from table with fallback to nearest V value"""
        tab_ic = TABLE.get(ic, {})
        tab_m = tab_ic.get(mm, {})
        
        if vv in tab_m:
            return tab_m[vv]
        
        if not tab_m:
            raise ValueError(f'Missing table data for Icl={ic}, M={mm}')
        
        # Fallback: find nearest velocity
        vv2, ab = min(tab_m.items(), key=lambda kv: abs(kv[0] - vv))
        return ab
    
    # 3D interpolation in (Icl, M, V) space
    # Interpolate in V first, then M, then Icl
    a_ll, b_ll = get_ab_at(ic_lo, m_lo, v_lo)
    a_lh, b_lh = get_ab_at(ic_lo, m_lo, v_hi)
    a_hl, b_hl = get_ab_at(ic_lo, m_hi, v_lo)
    a_hh, b_hh = get_ab_at(ic_lo, m_hi, v_hi)
    
    a_lm = linear_interpolate(a_ll, a_lh, t_v)
    b_lm = linear_interpolate(b_ll, b_lh, t_v)
    a_hm = linear_interpolate(a_hl, a_hh, t_v)
    b_hm = linear_interpolate(b_hl, b_hh, t_v)
    
    a_im = linear_interpolate(a_lm, a_hm, t_m)
    b_im = linear_interpolate(b_lm, b_hm, t_m)
    
    # Interpolate for higher Icl
    a_ll2, b_ll2 = get_ab_at(ic_hi, m_lo, v_lo)
    a_lh2, b_lh2 = get_ab_at(ic_hi, m_lo, v_hi)
    a_hl2, b_hl2 = get_ab_at(ic_hi, m_hi, v_lo)
    a_hh2, b_hh2 = get_ab_at(ic_hi, m_hi, v_hi)
    
    a_lm2 = linear_interpolate(a_ll2, a_lh2, t_v)
    b_lm2 = linear_interpolate(b_ll2, b_lh2, t_v)
    a_hm2 = linear_interpolate(a_hl2, a_hh2, t_v)
    b_hm2 = linear_interpolate(b_hl2, b_hh2, t_v)
    
    a_im2 = linear_interpolate(a_lm2, a_hm2, t_m)
    b_im2 = linear_interpolate(b_lm2, b_hm2, t_m)
    
    # Final interpolation in Icl
    a = linear_interpolate(a_im, a_im2, t_ic)
    b = linear_interpolate(b_im, b_im2, t_ic)
    
    return a, b


def acceptable_line_params(icl_clo, metabolic_rate, jet_velocity, mean_radiant_temp):
    """
    Calculate acceptable line parameters (m, C) from ASHRAE equations.
    Acceptable line: P = -m*T + C (pressure vs temperature)
    
    Returns: tuple (m, C, Ta50)
    """
    from models.psychrometric import (
        convective_heat_transfer_coeff, radiative_heat_transfer_coeff,
        clothing_factor, fraction_of_body_covered_by_clothing,
        local_effect_factor, saturation_vapor_pressure_mmhg
    )
    
    Ta50 = regress_ta50_from_equation(mean_radiant_temp, icl_clo, metabolic_rate, jet_velocity)
    
    hc = convective_heat_transfer_coeff(jet_velocity)
    hr = radiative_heat_transfer_coeff(mean_radiant_temp)
    h = hc + hr
    fcl = clothing_factor(icl_clo)
    Fcl = fraction_of_body_covered_by_clothing(icl_clo, jet_velocity, mean_radiant_temp)
    Fpcl = local_effect_factor(jet_velocity, icl_clo, mean_radiant_temp)
    
    m = (fcl * Fcl) / max(1e-12, (1.1 * Fpcl))
    Ps_ta50 = saturation_vapor_pressure_mmhg(Ta50)
    C = m * Ta50 + 0.5 * Ps_ta50
    
    return m, C, Ta50


def regress_ta50_from_equation(mean_radiant_temp, icl_clo, metabolic_rate, jet_velocity):
    """Calculate Ta(0.5) from regression: Ta50 = a*Tmr + b"""
    a, b = regress_ta50_coefficients(icl_clo, metabolic_rate, jet_velocity)
    return a * mean_radiant_temp + b


def jet_velocity_and_temperature_ratios(axial_distance_x0, jet_diameter_d0,
                                        include_buoyancy, ambient_temp, nozzle_temp,
                                        velocity_guess=10.0):
    """
    Calculate velocity ratio Vj/V0 and temperature ratio (TA-Tj)/(TA-T0).
    Based on ASHRAE jet decay model with optional buoyancy correction.
    
    Returns: tuple (velocity_ratio, temperature_ratio)
    """
    r = (axial_distance_x0 / jet_diameter_d0) + 2.572
    
    if not include_buoyancy:
        velocity_ratio = 1.464 / r
        temperature_ratio = 4.539 / r
        return velocity_ratio, temperature_ratio
    
    # With buoyancy (Archimedes number)
    beta = 1.0 / 273.15  # 1/T_kelvin
    dT = max(0.0, ambient_temp - nozzle_temp)
    Ar = (9.81 * beta * dT * jet_diameter_d0) / max(1e-6, velocity_guess ** 2)
    
    buoyancy_factor = (1.0 + 0.21 * Ar * r * r) ** (1.0 / 3.0)
    velocity_ratio = (1.464 / r) * buoyancy_factor
    temperature_ratio = (4.539 / r) * buoyancy_factor
    
    return velocity_ratio, temperature_ratio


def solve_spot_cooling(ambient_temp, ambient_rh_fraction, mean_radiant_temp,
                       target_velocity, metabolic_rate, clothing_icl,
                       jet_diameter, axial_distance,
                       nozzle_rh_fraction=0.95, p_atm_kpa=101.325,
                       include_buoyancy=False):
    """
    Solve spot cooling design problem using ASHRAE Standard RP-884.
    
    Args:
        ambient_temp: Ambient air temperature [°C]
        ambient_rh_fraction: Ambient relative humidity [0-1]
        mean_radiant_temp: Mean radiant temperature [°C]
        target_velocity: Target jet velocity at measurement area [m/s]
        metabolic_rate: Metabolic rate [W/m²]
        clothing_icl: Clothing insulation [clo]
        jet_diameter: Jet nozzle diameter [m]
        axial_distance: Axial distance from nozzle [m]
        nozzle_rh_fraction: Humidity at nozzle exit [0-1], default 0.95
        p_atm_kpa: Atmospheric pressure [kPa]
        include_buoyancy: Include buoyancy effects in jet model
    
    Returns: dict with solution (T0, P0, Tj, Pj, RH_j, V0, Q0, etc.)
    """
    P_ambient = vapor_pressure_from_rh_and_temp(ambient_rh_fraction, ambient_temp)
    m, C, Ta50 = acceptable_line_params(clothing_icl, metabolic_rate, target_velocity, mean_radiant_temp)
    
    # Binary search for T0 satisfying all constraints
    T0_lo = -10.0
    T0_hi = ambient_temp - 0.1
    
    def residual_function(T0_test):
        P0 = nozzle_rh_fraction * saturation_vapor_pressure_mmhg(T0_test)
        velocity_ratio, temperature_ratio = jet_velocity_and_temperature_ratios(
            axial_distance, jet_diameter, include_buoyancy, ambient_temp, T0_test
        )
        Tj = ambient_temp - temperature_ratio * (ambient_temp - T0_test)
        
        # Operating line A->0
        if abs(T0_test - ambient_temp) < 1e-9:
            return 1e9
        
        slope_A0 = (P0 - P_ambient) / (T0_test - ambient_temp)
        Pj_operating = P_ambient + slope_A0 * (Tj - ambient_temp)
        
        # Acceptable line at target Vj
        Pj_acceptable = -m * Tj + C
        
        return Pj_operating - Pj_acceptable
    
    # Find bracket
    f_lo = residual_function(T0_lo)
    f_hi = residual_function(T0_hi)
    
    if f_lo * f_hi > 0:
        # Try alternative brackets
        for (a, b) in [(-20.0, ambient_temp - 0.1), (-5.0, ambient_temp - 0.1),
                       (0.0, ambient_temp - 0.5)]:
            if residual_function(a) * residual_function(b) <= 0:
                T0_lo, T0_hi = a, b
                break
    
    # Bisection
    for _ in range(100):
        T0 = 0.5 * (T0_lo + T0_hi)
        f_T0 = residual_function(T0)
        
        if abs(f_T0) < 1e-4 or (T0_hi - T0_lo) < 1e-4:
            break
        
        if residual_function(T0_lo) * f_T0 <= 0:
            T0_hi = T0
        else:
            T0_lo = T0
    
    # Extract solution
    P0 = nozzle_rh_fraction * saturation_vapor_pressure_mmhg(T0)
    velocity_ratio, temperature_ratio = jet_velocity_and_temperature_ratios(
        axial_distance, jet_diameter, include_buoyancy, ambient_temp, T0
    )
    Tj = ambient_temp - temperature_ratio * (ambient_temp - T0)
    Pj = P_ambient - temperature_ratio * (P_ambient - P0)
    
    RH_j = relative_humidity_from_vapor_pressure(Pj, Tj)
    V0 = target_velocity / max(1e-9, velocity_ratio)
    
    # Flow rates
    area0 = math.pi * (jet_diameter ** 2) / 4.0
    Q0 = V0 * area0
    Qj = Q0 / max(1e-9, temperature_ratio)
    Qe = max(0.0, Qj - Q0)
    
    # Energy
    w_ambient = humidity_ratio_from_vapor_pressure(P_ambient, p_atm_kpa)
    w0 = humidity_ratio_from_vapor_pressure(P0, p_atm_kpa)
    h_ambient = enthalpy_moist_air(ambient_temp, w_ambient)
    h0 = enthalpy_moist_air(T0, w0)
    rho0 = air_density(T0, w0, p_atm_kpa)
    m_dot0 = rho0 * Q0
    Q_total = m_dot0 * (h_ambient - h0)
    Q_sensible = m_dot0 * 1.006 * (ambient_temp - T0)
    
    return {
        'm': m, 'C': C, 'Ta50': Ta50,
        'PA': P_ambient, 'Pj': Pj, 'P0': P0,
        'Tj': Tj, 'T0': T0, 'RH_j': RH_j, 'RH_0': nozzle_rh_fraction,
        'Vratio': velocity_ratio, 'Tratio': temperature_ratio,
        'V0': V0, 'Q0': Q0, 'Qj': Qj, 'Qe': Qe,
        'm_dot0': m_dot0, 'Q_total': Q_total, 'Q_sens': Q_sensible
    }


def solve_spot_cooling_from_target_point(ambient_temp, ambient_rh_fraction,
                                         mean_radiant_temp, target_velocity,
                                         metabolic_rate, clothing_icl,
                                         jet_diameter, axial_distance,
                                         target_temp, target_humidity_ratio,
                                         p_atm_kpa=101.325, include_buoyancy=False):
    """
    Inverse solve: find nozzle conditions to achieve target (T_j, W_j) point.
    
    Returns: dict with solution
    """
    P_ambient = vapor_pressure_from_rh_and_temp(ambient_rh_fraction, ambient_temp)
    Pj = vapor_pressure_from_humidity_ratio(target_humidity_ratio, p_atm_kpa)
    m, C, Ta50 = acceptable_line_params(clothing_icl, metabolic_rate,
                                       target_velocity, mean_radiant_temp)
    
    # Iterative solution with buoyancy
    if include_buoyancy:
        T0 = target_temp
        for _ in range(50):
            velocity_ratio, temperature_ratio = jet_velocity_and_temperature_ratios(
                axial_distance, jet_diameter, True, ambient_temp, T0
            )
            T0_new = ambient_temp - (ambient_temp - target_temp) / max(1e-9, temperature_ratio)
            if abs(T0_new - T0) < 1e-5:
                T0 = T0_new
                break
            T0 = T0_new
        velocity_ratio, temperature_ratio = jet_velocity_and_temperature_ratios(
            axial_distance, jet_diameter, True, ambient_temp, T0
        )
    else:
        velocity_ratio, temperature_ratio = jet_velocity_and_temperature_ratios(
            axial_distance, jet_diameter, False, ambient_temp, target_temp
        )
        T0 = ambient_temp - (ambient_temp - target_temp) / max(1e-9, temperature_ratio)
    
    P0 = P_ambient - (P_ambient - Pj) / max(1e-9, temperature_ratio)
    RH_j = relative_humidity_from_vapor_pressure(Pj, target_temp)
    RH_0 = relative_humidity_from_vapor_pressure(P0, T0)
    
    # Flow rates
    V0 = target_velocity / max(1e-9, velocity_ratio)
    area0 = math.pi * (jet_diameter ** 2) / 4.0
    Q0 = V0 * area0
    Qj = Q0 / max(1e-9, temperature_ratio)
    Qe = max(0.0, Qj - Q0)
    
    # Energy
    w_ambient = humidity_ratio_from_vapor_pressure(P_ambient, p_atm_kpa)
    w0 = humidity_ratio_from_vapor_pressure(P0, p_atm_kpa)
    h_ambient = enthalpy_moist_air(ambient_temp, w_ambient)
    h0 = enthalpy_moist_air(T0, w0)
    rho0 = air_density(T0, w0, p_atm_kpa)
    m_dot0 = rho0 * Q0
    Q_total = m_dot0 * (h_ambient - h0)
    Q_sensible = m_dot0 * 1.006 * (ambient_temp - T0)
    
    return {
        'm': m, 'C': C, 'Ta50': Ta50,
        'PA': P_ambient, 'Pj': Pj, 'P0': P0,
        'Tj': target_temp, 'T0': T0, 'RH_j': RH_j, 'RH_0': RH_0,
        'Vratio': velocity_ratio, 'Tratio': temperature_ratio,
        'V0': V0, 'Q0': Q0, 'Qj': Qj, 'Qe': Qe,
        'm_dot0': m_dot0, 'Q_total': Q_total, 'Q_sens': Q_sensible
    }
