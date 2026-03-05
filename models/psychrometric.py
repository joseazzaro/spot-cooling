# -*- coding: utf-8 -*-
"""
Psychrometric calculations for humid air.
Based on ASHRAE Standard RP-884 and standard psychrometric formulas.
"""

import math
from utils.constants import (
    MMHG_TO_KPA, PSY_COEFF_A, PSY_COEFF_B, PSY_COEFF_C,
    SPECIFIC_HEAT_MOIST_AIR_kJ_kg_K, LATENT_HEAT_WATER,
    WATER_VAPOR_COEFF, HUMIDITY_RATIO_COEFF, GAS_CONSTANT_DRY_AIR,
    GRAVITY, ABS_ZERO_CELSIUS
)


def saturation_vapor_pressure_mmhg(T_celsius):
    """
    Saturation vapor pressure of water at temperature T [°C].
    Returns pressure in mmHg.
    Uses ASHRAE formula: Psat = exp(A - B/(C+T))
    """
    return math.exp(PSY_COEFF_A - PSY_COEFF_B / (T_celsius + PSY_COEFF_C))


def vapor_pressure_from_rh_and_temp(relative_humidity_fraction, T_celsius):
    """
    Vapor pressure from relative humidity (0-1) and temperature.
    Returns pressure in mmHg.
    """
    return relative_humidity_fraction * saturation_vapor_pressure_mmhg(T_celsius)


def relative_humidity_from_vapor_pressure(vapor_pressure_mmhg, T_celsius):
    """
    Relative humidity (0-1) from vapor pressure and temperature.
    """
    return vapor_pressure_mmhg / saturation_vapor_pressure_mmhg(T_celsius)


def humidity_ratio_from_vapor_pressure(vapor_pressure_mmhg, p_atm_kpa=101.325):
    """
    Humidity ratio w [kg_water/kg_dry_air] from vapor pressure.
    Uses ASHRAE definition: w = 0.62198 * Pv / (P_atm - Pv)
    """
    vapor_pressure_kpa = vapor_pressure_mmhg * MMHG_TO_KPA
    denominator = max(1e-9, p_atm_kpa - vapor_pressure_kpa)
    return HUMIDITY_RATIO_COEFF * vapor_pressure_kpa / denominator


def vapor_pressure_from_humidity_ratio(humidity_ratio, p_atm_kpa=101.325):
    """
    Vapor pressure [mmHg] from humidity ratio w [kg_water/kg_dry_air].
    Inverse of humidity_ratio_from_vapor_pressure.
    """
    vapor_pressure_kpa = (humidity_ratio * p_atm_kpa) / max(1e-9, (HUMIDITY_RATIO_COEFF + humidity_ratio))
    return vapor_pressure_kpa / MMHG_TO_KPA


def enthalpy_moist_air(T_celsius, humidity_ratio):
    """
    Enthalpy of moist air [kJ/kg_dry_air] at temperature T and humidity ratio w.
    Formula: h = cp*T + w*(L + cp_v*T)
    """
    return (SPECIFIC_HEAT_MOIST_AIR_kJ_kg_K * T_celsius + 
            humidity_ratio * (LATENT_HEAT_WATER + WATER_VAPOR_COEFF * T_celsius))


def air_density(T_celsius, humidity_ratio, p_atm_kpa=101.325):
    """
    Density of moist air [kg/m³] at temperature T, humidity ratio w, and pressure P.
    Uses ideal gas law with correction for moisture.
    """
    T_kelvin = T_celsius + ABS_ZERO_CELSIUS
    p_pa = p_atm_kpa * 1000.0
    return p_pa / (GAS_CONSTANT_DRY_AIR * T_kelvin * (1.0 + 1.607 * humidity_ratio))


def convective_heat_transfer_coeff(jet_velocity=1.0):
    """
    Convective heat transfer coefficient [W/m²K].
    Uses simplified boundary layer formula: hc = 8.3 * V^0.5
    """
    return 8.3 * math.sqrt(max(1e-6, jet_velocity))


def radiative_heat_transfer_coeff(mean_radiant_temp=45.0):
    """
    Radiative heat transfer coefficient [W/m²K].
    Uses linearized radiation: hr = 3.87 + 0.031*Tmr
    """
    return 3.87 + 0.031 * mean_radiant_temp


def total_heat_transfer_coeff(jet_velocity, mean_radiant_temp):
    """Total heat transfer coefficient = convective + radiative"""
    hc = convective_heat_transfer_coeff(jet_velocity)
    hr = radiative_heat_transfer_coeff(mean_radiant_temp)
    return hc + hr


def clothing_factor(icl_clo):
    """Clothing factor 1 + k*Icl"""
    return 1.0 + 0.2 * icl_clo


def fraction_of_body_covered_by_clothing(icl_clo, jet_velocity, mean_radiant_temp):
    """
    Fraction of body covered by clothing (effective).
    Formula: Fcl = 1 / (1 + 0.155*fcl*h*Icl)
    """
    fcl = clothing_factor(icl_clo)
    h = total_heat_transfer_coeff(jet_velocity, mean_radiant_temp)
    return 1.0 / (1.0 + 0.155 * fcl * h * icl_clo)


def local_effect_factor(jet_velocity, icl_clo, mean_radiant_temp):
    """
    Local effect factor for unclothed parts.
    Formula: Fpcl = 1 / (1 + 0.143*hc*Icl)
    """
    hc = convective_heat_transfer_coeff(jet_velocity)
    return 1.0 / (1.0 + 0.143 * hc * icl_clo)
