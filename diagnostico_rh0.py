# -*- coding: utf-8 -*-
"""
Diagnóstico: ¿Por qué RH0=0.90 y 1.00 no funcionan?
Analizar el residual en diferentes rangos para cada RH0
"""

import sys
import math
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import (psat_mmhg, jet_ratios, acceptable_line, pv_from_rh_T)

print("="*95)
print("DIAGNÓSTICO: Análisis del residual para cada RH0")
print("="*95)

# Parámetros constantes
TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.3048  # 1 ft
X0 = 3.048   # 10 ft

# Parámetros que se calculan una sola vez
PA = pv_from_rh_T(RH_A, TA)
m, C, Ta50 = acceptable_line(Icl, M, Vj, Tmr)

print(f"\nParámetros comunes:")
print(f"  PA = {PA:.2f} mm Hg")
print(f"  m = {m:.6f} mm Hg/°C")
print(f"  C = {C:.2f} mm Hg")

print(f"\n" + "="*95)
print("ANÁLISIS POR RH0:")
print("="*95)

# Literatura TABLA 2
tabla2 = {
    0.90: {"T0": 19.40, "Tj": 32.85, "Pj": 21.50},
    0.95: {"T0": 26.96, "Tj": 35.25, "Pj": 26.87},
    1.00: {"T0": 34.98, "Tj": 37.28, "Pj": 32.02},
}

for rh0 in [0.90, 0.95, 1.00]:
    print(f"\n" + "-"*95)
    print(f"RH0 = {rh0*100:.0f}% (Esperado: T0 = {tabla2[rh0]['T0']}°C)")
    print("-"*95)
    
    def residual_T0(T0):
        P0 = rh0*psat_mmhg(T0)
        Vratio, Tratio = jet_ratios(X0, D0, False, TA, T0)
        Ti = TA - Tratio*(TA - T0)
        Pi = PA - Tratio*(PA - P0)
        return Pi - (-m*Ti + C)
    
    # Probar en diferentes puntos
    T0_test_points = [10, 15, 19.4, 20, 26.96, 30, 34.98, 35, 40]
    
    print(f"\n{'T0 (°C)':<12} {'P0':<12} {'Vratio':<12} {'Tratio':<12} {'Ti':<12} {'Pi':<12} {'Residual':<12}")
    print("-"*95)
    
    for T0_test in T0_test_points:
        P0 = rh0*psat_mmhg(T0_test)
        Vratio, Tratio = jet_ratios(X0, D0, False, TA, T0_test)
        Ti = TA - Tratio*(TA - T0_test)
        Pi = PA - Tratio*(PA - P0)
        res = residual_T0(T0_test)
        
        print(f"{T0_test:<12.2f} {P0:<12.2f} {Vratio:<12.6f} {Tratio:<12.6f} {Ti:<12.2f} {Pi:<12.2f} {res:<12.4f}")
    
    # Buscar dónde cambia de signo
    print(f"\nBúsqueda de cambios de signo:")
    
    ranges_to_test = [
        (-5.0, 40.0),
        (-20.0, 40.0),
        (-10.0, 50.0),
        (0.0, 60.0),
        (-30.0, 50.0),
        (0.0, 40.0),
    ]
    
    found = False
    for a, b in ranges_to_test:
        f_a = residual_T0(a)
        f_b = residual_T0(b)
        sign_change = f_a * f_b <= 0
        print(f"  Rango [{a:6.1f}, {b:6.1f}]: f({a:6.1f})={f_a:8.2f}, f({b:6.1f})={f_b:8.2f} → Cambio de signo: {sign_change}")
        if sign_change and not found:
            found = True
            print(f"    ✓ ENCONTRADO cambio de signo en este rango")
