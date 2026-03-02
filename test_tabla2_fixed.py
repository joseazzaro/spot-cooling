# -*- coding: utf-8 -*-
"""
Test: Verificar si ahora obtenemos los resultados correctos de TABLA 2 ASHRAE
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import solve_case

print("="*80)
print("TEST: TABLA 2 ASHRAE - Ejemplo 1 (RH0 = 0.95)")
print("="*80)

# Parámetros del ejemplo
TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.3048  # 1 ft
X0 = 0.3048  # 1 ft (asumiendo que es igual a D0)

print(f"\nEntradas:")
print(f"  TA = {TA}°C")
print(f"  RH_A = {RH_A*100}%")
print(f"  Tmr = {Tmr}°C")
print(f"  Vj = {Vj} m/s")
print(f"  M = {M} W/m²")
print(f"  Icl = {Icl} clo")
print(f"  D0 = {D0} m")
print(f"  X0 = {X0} m")

# Casos: TABLA 2 tiene RH0 = 0.90, 0.95, 1.00
print(f"\n" + "="*80)
print("RESULTADOS ESPERADOS (TABLA 2 ASHRAE):")
print("="*80)

expected = {
    0.90: {"T0": 19.40, "Tj": 32.85, "Pj": 21.50},
    0.95: {"T0": 26.96, "Tj": 35.25, "Pj": 26.87},
    1.00: {"T0": 34.98, "Tj": 37.28, "Pj": 32.02},
}

for rh0, values in expected.items():
    print(f"\nRH0 = {rh0*100}%:")
    print(f"  T0 = {values['T0']}°C (esperado)")
    print(f"  Tj = {values['Tj']}°C (esperado)")
    print(f"  Pj = {values['Pj']} mm Hg (esperado)")

# Ahora calcular con nuestro código corregido
print(f"\n" + "="*80)
print("RESULTADOS CON CÓDIGO CORREGIDO:")
print("="*80)

for rh0, expected_vals in expected.items():
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
    
    T0 = result['T0']
    Tj = result['Tj']
    Pj = result['Pj']
    
    print(f"\nRH0 = {rh0*100}%:")
    print(f"  T0 = {T0:.2f}°C (calculado) vs {expected_vals['T0']}°C (esperado) → Error: {abs(T0-expected_vals['T0']):.2f}°C")
    print(f"  Tj = {Tj:.2f}°C (calculado) vs {expected_vals['Tj']}°C (esperado) → Error: {abs(Tj-expected_vals['Tj']):.2f}°C")
    print(f"  Pj = {Pj:.2f} mm Hg (calculado) vs {expected_vals['Pj']} mm Hg (esperado) → Error: {abs(Pj-expected_vals['Pj']):.2f} mm Hg")

print(f"\n" + "="*80)
print("DIAGNÓSTICO:")
print("="*80)

result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=0.95)
print(f"\nDetalles para RH0 = 0.95:")
print(f"  m = {result['m']:.6f} mm Hg/°C (correcto: 0.737)")
print(f"  C = {result['C']:.2f} mm Hg (correcto: 52.85)")
print(f"  Ta(0.5) = {result['Ta50']:.2f}°C")
print(f"  PA = {result['PA']:.2f} mm Hg")
print(f"  Pj (presión en jet) = {result['Pj']:.2f} mm Hg")
print(f"  P0 (presión ambiente) = {result['P0']:.2f} mm Hg")
