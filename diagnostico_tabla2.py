# -*- coding: utf-8 -*-
"""
Script de diagnóstico: comparar resultados con Tabla 2 ASHRAE
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import (
    psat_mmhg, pv_from_rh_T, rh_from_pv_T, humidity_ratio_from_pv,
    moist_air_enthalpy_kJkg, air_density_approx, compute_factors,
    ta50_from_regression, acceptable_line, jet_ratios, solve_case
)
import math

print("="*80)
print("DIAGNÓSTICO: COMPARACIÓN CON TABLA 2 ASHRAE")
print("="*80)

# Parámetros que usamos (asumiendo que son los mismos del ejemplo)
TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.127
X0 = 1.259

print("\nParámetros de entrada (ASUMIDOS):")
print(f"  TA = {TA}°C")
print(f"  RH_A = {RH_A*100}%")
print(f"  Tmr = {Tmr}°C")
print(f"  Vj = {Vj} m/s")
print(f"  M = {M} W/m²")
print(f"  Icl = {Icl} clo")

print("\n" + "="*80)
print("TABLA 2 ASHRAE (VALORES ESPERADOS):")
print("="*80)
print("\nrh₀    T₀, °C   Tⱼ, °C   Pⱼ, mm Hg")
print("-" * 40)
print("0.90   27.58    35.48    26.33")
print("0.95   26.96    35.25    26.87")
print("1.00   26.38    35.04    27.03")

print("\n" + "="*80)
print("NUESTROS RESULTADOS:")
print("="*80)

for rh0 in [0.90, 0.95, 1.00]:
    print(f"\nrh₀ = {rh0}:")
    
    try:
        result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0, p_atm_kpa=101.325)
        
        T0 = result['T0']
        Ti = result['Tj']
        Pi = result['Pj']
        
        print(f"  T₀ = {T0:.2f}°C  (esperado: ~26.96°C, error: {T0-26.96:.2f}°C)")
        print(f"  Tⱼ = {Ti:.2f}°C  (esperado: ~35.25°C, error: {Ti-35.25:.2f}°C)")
        print(f"  Pⱼ = {Pi:.2f} mmHg (esperado: ~26.87 mmHg, error: {Pi-26.87:.2f} mmHg)")
        
    except Exception as e:
        print(f"  ERROR: {e}")

print("\n" + "="*80)
print("ANÁLISIS DEL PROBLEMA:")
print("="*80)

print("""
POSIBLES CAUSAS:

1. PARAMETROS INCORRECTOS:
   → ¿TA, RH_A, Tmr, M, Icl son realmente 40, 50%, 45, 87, 0.6?
   → El ejemplo ASHRAE podría tener valores diferentes

2. ALGORITMO DE BISECCIÓN:
   → El residual podría tener un signo incorrecto
   → La función podría converger al punto equivocado
   → Falta de convergencia a la solución correcta

3. ECUACIÓN DE LA LINEA ACEPTABLE:
   → Los valores m y C podrían estar incorrectos
   → O la forma de evaluar si está en la línea es incorrecta

4. INTERPRETACION DE VARIABLES:
   → Pⱼ podría significar algo diferente
   → Tⱼ podría ser otra variable termodinámica

ACCIÓN REQUERIDA:
→ Confirmar exactamente qué parámetros usa el ejemplo de Tabla 2
→ Verificar si hay otras condiciones o restricciones
→ Revisar la ecuación de residuos en solve_case()
""")

print("\n" + "="*80)
print("DEBUG: Analizar la función residual")
print("="*80)

# Intentar entender qué está pasando
PA = pv_from_rh_T(RH_A, TA)
print(f"\nPresión de vapor ambiente: PA = {PA:.2f} mmHg")

m, C, Ta50 = acceptable_line(Icl, M, Vj, Tmr)
print(f"Línea aceptable: Pv = {m:.6f} * Ta + {C:.3f}")

# Probar varios T0
print(f"\nEvaluando función residual en varios puntos:")
print("T0 (°C)  | Presión esperada | Presión en línea | Residual")
print("-" * 70)

for T0_test in [20, 25, 26.96, 30, 35, 40]:
    # Calcular como lo hace el código
    P_sat_T0 = psat_mmhg(T0_test)
    P0 = 0.95 * P_sat_T0
    
    r = (X0/D0) + 2.572
    Vratio = 1.464 / r
    Tratio = 4.539 / r
    
    Ti = TA - Tratio * (TA - T0_test)
    Pi = PA - Tratio * (PA - P0)
    
    # ¿Está en la línea?
    Pi_aceptable = m * Ti + C
    residual = Pi - Pi_aceptable
    
    print(f"{T0_test:6.2f}  | {Pi:16.2f} | {Pi_aceptable:16.2f} | {residual:10.2f}")

print("\nVER: ¿Dónde cambia de signo el residual?")
print("     Si cambia cerca de 27°C, eso explica el resultado esperado.")
print("     Si cambia cerca de 40°C, hay un problema en el algoritmo.")
