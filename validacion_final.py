# -*- coding: utf-8 -*-
"""
VALIDACIÓN FINAL: Comparar código corregido contra TABLA 2 ASHRAE completa
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import solve_case

print("\n" + "="*95)
print("VALIDACIÓN FINAL: TABLA 2 ASHRAE - Ejemplo 1 Completo")
print("="*95)

# Parámetros del ejemplo (del documento ASHRAE)
TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.3048  # 1 ft
X0 = 3.048   # 10 ft (CORRECTO)

print(f"\nParámetros del ejemplo ASHRAE:")
print(f"  Ambient: TA = {TA}°C, RH_A = {RH_A*100}%, Tmr = {Tmr}°C")
print(f"  Jet: Vj = {Vj} m/s, D0 = {D0} m (1 ft)")
print(f"  Work: M = {M} W/m², Icl = {Icl} clo")
print(f"  Geometry: X0 = {X0} m (10 ft) - altura vertical del jet")

# TABLA 2 ASHRAE: Resultados esperados para tres valores de RH0
tabla2_expected = {
    0.90: {"T0": 19.40, "Tj": 32.85, "Pj": 21.50},
    0.95: {"T0": 26.96, "Tj": 35.25, "Pj": 26.87},
    1.00: {"T0": 34.98, "Tj": 37.28, "Pj": 32.02},
}

print(f"\n" + "="*95)
print("TABLA 2 ASHRAE: RESULTADOS ESPERADOS vs CALCULADOS")
print("="*95)

print(f"\n{'RH0':<8} {'Parámetro':<12} {'Esperado':<15} {'Calculado':<15} {'Error':<15} {'Status':<8}")
print("-" * 95)

all_ok = True

for rh0 in [0.90, 0.95, 1.00]:
    expected = tabla2_expected[rh0]
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
    
    T0_calc = result['T0']
    Tj_calc = result['Tj']
    Pj_calc = result['Pi']
    
    err_T0 = abs(T0_calc - expected['T0'])
    err_Tj = abs(Tj_calc - expected['Tj'])
    err_Pj = abs(Pj_calc - expected['Pj'])
    
    # Tolerancia aceptable ±0.1 en temperatura, ±0.1 en presión
    status_T0 = "✓" if err_T0 < 0.1 else "✗"
    status_Tj = "✓" if err_Tj < 0.1 else "✗"
    status_Pj = "✓" if err_Pj < 0.1 else "✗"
    
    if err_T0 >= 0.1 or err_Tj >= 0.1 or err_Pj >= 0.1:
        all_ok = False
    
    print(f"{rh0*100:<8.0f} {'T0':<12} {expected['T0']:<15.2f} {T0_calc:<15.2f} {err_T0:<15.4f} {status_T0:<8}")
    print(f"{'': <8} {'Tj':<12} {expected['Tj']:<15.2f} {Tj_calc:<15.2f} {err_Tj:<15.4f} {status_Tj:<8}")
    print(f"{'': <8} {'Pj':<12} {expected['Pj']:<15.2f} {Pj_calc:<15.2f} {err_Pj:<15.4f} {status_Pj:<8}")
    print()

print("="*95)
if all_ok:
    print("✓ VALIDACIÓN EXITOSA: Todos los valores están dentro de tolerancia ±0.1")
else:
    print("⚠ VALIDACIÓN PARCIAL: Algunos valores están fuera de tolerancia")
print("="*95)

# Detalles técnicos para RH0=0.95
print(f"\nDETALLES TÉCNICOS (RH0 = 0.95):")
print("="*95)
result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=0.95)

print(f"\nParámetros calculados:")
print(f"  m (pendiente línea aceptable) = {result['m']:.6f} mm Hg/°C")
print(f"  C (intercepto línea aceptable) = {result['C']:.2f} mm Hg")
print(f"  Ta(0.5) = {result['Ta50']:.2f}°C")
print(f"\nCondiciones ambiente:")
print(f"  PA (presión vapor ambiente) = {result['PA']:.2f} mm Hg")
print(f"  P0 (presión vapor en orificio) = {result['P0']:.2f} mm Hg")
print(f"  RH0 = {result['RH_0']*100:.0f}%")
print(f"\nCondiciones en zona de impacto (target area):")
print(f"  T0 (temperatura en orificio) = {result['T0']:.2f}°C")
print(f"  Tj = Tj (temperatura en zona impacto) = {result['Tj']:.2f}°C")
print(f"  Pj = Pj (presión vapor en zona impacto) = {result['Pj']:.2f} mm Hg")
print(f"  RHj = {result['RH_j']*100:.1f}%")
print(f"\nRatios del jet:")
print(f"  Vratio = {result['Vratio']:.6f}")
print(f"  Tratio = {result['Tratio']:.6f}")
print(f"\nBalance térmico:")
print(f"  m_dot0 = {result['m_dot0']:.4f} kg/s")
print(f"  Q_total = {result['Q_total']:.2f} W")
print(f"  Q_sens = {result['Q_sens']:.2f} W")
