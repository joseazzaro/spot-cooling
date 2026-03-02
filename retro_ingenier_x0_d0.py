# -*- coding: utf-8 -*-
"""
Retro-ingeniería: Encontrar qué X0 y D0 producen T0=26.96°C
"""

import sys
import math
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import (solve_case, jet_ratios, acceptable_line, pv_from_rh_T)

print("="*80)
print("RETRO-INGENIERÍA: ¿Cuáles son X0 y D0?")
print("="*80)

TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
rh0 = 0.95

# Resultado esperado
T0_expected = 26.96

# Intentar diferentes combinaciones X0/D0
print(f"\nBuscando X0 y D0 para obtener T0 = {T0_expected}°C:\n")
print(f"{'X0 (m)':<10} {'D0 (m)':<10} {'X0/D0':<10} {'r':<10} {'Vratio':<12} {'Tratio':<12} {'T0_calc':<12} {'Error':<12}")
print("-" * 110)

# Datos posibles
# El usuario mencionó: "radius Rc = 0.3048 m (1 ft)"
# Si Rc es el radio, el diámetro sería 2*Rc = 0.6096 m
# O si Rc es ya el radio de la zona de enfriamiento...

# Probar diferentes valores
test_cases = [
    # Asumiendo variantes de 0.3048 m
    (0.3048, 0.3048, "D0=X0=1ft"),
    (0.6096, 0.3048, "D0=2ft, X0=1ft"),
    (0.3048, 0.6096, "D0=1ft, X0=2ft"),
    (1.5240, 0.3048, "D0=5ft, X0=1ft"),
    (0.3048, 1.5240, "D0=1ft, X0=5ft"),
    (2.4384, 0.3048, "D0=8ft, X0=1ft"),
    (0.3048, 2.4384, "D0=1ft, X0=8ft"),
    # Otras combinaciones
    (0.1524, 0.3048, "D0=0.5ft, X0=1ft"),
    (0.3048, 0.1524, "D0=1ft, X0=0.5ft"),
]

best_match = None
best_error = float('inf')

for X0, D0, label in test_cases:
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
    T0_calc = result['T0']
    error = abs(T0_calc - T0_expected)
    
    r = (X0/D0) + 2.572
    Vratio = 1.464/r
    Tratio = 4.539/r
    
    print(f"{X0:<10.4f} {D0:<10.4f} {X0/D0:<10.4f} {r:<10.4f} {Vratio:<12.6f} {Tratio:<12.6f} {T0_calc:<12.2f} {error:<12.2f}")
    
    if error < best_error:
        best_error = error
        best_match = (X0, D0, T0_calc, label)

print("\n" + "="*80)
if best_match and best_match[2]:
    X0_best, D0_best, T0_best, label_best = best_match
    print(f"MEJOR COINCIDENCIA: {label_best}")
    print(f"  X0 = {X0_best:.4f} m")
    print(f"  D0 = {D0_best:.4f} m")
    print(f"  T0 calculado = {T0_best:.2f}°C")
    print(f"  T0 esperado = {T0_expected:.2f}°C")
    print(f"  Error = {abs(T0_best - T0_expected):.2f}°C")
    
    # Ver si esto da los otros valores correctos también
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0_best, X0_best, rh0=rh0)
    print(f"\n  Tj = {result['Tj']:.2f}°C (esperado 35.25°C)")
    print(f"  Pj = {result['Pj']:.2f} mm Hg (esperado 26.87 mm Hg)")
