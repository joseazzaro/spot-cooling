# -*- coding: utf-8 -*-
"""
Busca fina: Encontrar X0 y D0 más precisos alrededor de la mejor coincidencia
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import solve_case

print("="*80)
print("BÚSQUEDA FINA: Ajustar X0 alrededor de 2.4384 m")
print("="*80)

TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
rh0 = 0.95

T0_expected = 26.96
Tj_expected = 35.25
Pj_expected = 26.87

D0 = 0.3048  # Fijo en 1 ft

print(f"\nD0 = {D0} m (1 ft) [fijo]")
print(f"Variando X0:\n")
print(f"{'X0 (m)':<12} {'X0 (ft)':<12} {'T0_calc':<12} {'Error_T0':<12} {'Tj_calc':<12} {'Error_Tj':<12} {'Pj_calc':<12} {'Error_Pj':<12}")
print("-" * 100)

best_overall = None
best_overall_error = float('inf')

# Buscar alrededor de 2.4384
for X0_ft in [6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]:
    X0 = X0_ft * 0.3048
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
    
    T0_calc = result['T0']
    Tj_calc = result['Tj']
    Pj_calc = result['Pj']
    
    error_T0 = abs(T0_calc - T0_expected)
    error_Tj = abs(Tj_calc - Tj_expected)
    error_Pj = abs(Pj_calc - Pj_expected)
    
    # Error combinado (promedio ponderado)
    error_combined = (error_T0 * 0.5 + error_Tj * 0.3 + error_Pj * 0.2)
    
    print(f"{X0:<12.4f} {X0_ft:<12.1f} {T0_calc:<12.2f} {error_T0:<12.2f} {Tj_calc:<12.2f} {error_Tj:<12.2f} {Pj_calc:<12.2f} {error_Pj:<12.2f}")
    
    if error_combined < best_overall_error:
        best_overall_error = error_combined
        best_overall = (X0, X0_ft, T0_calc, Tj_calc, Pj_calc)

print("\n" + "="*80)
print(f"MEJOR RESULTADO GENERAL:")
print("="*80)
if best_overall:
    X0, X0_ft, T0, Tj, Pj = best_overall
    print(f"\nX0 = {X0:.4f} m ({X0_ft:.1f} ft)")
    print(f"D0 = {D0:.4f} m (1.0 ft)")
    print(f"\nResultados:")
    print(f"  T0 = {T0:.2f}°C (esperado: 26.96°C, error: {abs(T0-26.96):.2f}°C)")
    print(f"  Tj = {Tj:.2f}°C (esperado: 35.25°C, error: {abs(Tj-35.25):.2f}°C)")
    print(f"  Pj = {Pj:.2f} mm Hg (esperado: 26.87 mm Hg, error: {abs(Pj-26.87):.2f} mm Hg)")
    
    # Ahora probar el mismo ratio X0/D0 con D0 = 0.6096 (si Rc es radio)
    print(f"\n" + "="*80)
    print("Al cuadrado: Si D0 = 0.6096 m (si Rc=0.3048 es radio, no diámetro):")
    print("="*80)
    D0_alt = 0.6096
    # Mantener el ratio X0/D0
    ratio = X0 / D0
    X0_alt = ratio * D0_alt
    
    result_alt = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0_alt, X0_alt, rh0=rh0)
    T0_alt = result_alt['T0']
    Tj_alt = result_alt['Ti']
    Pj_alt = result_alt['Pi']
    
    print(f"\nManteniendo ratio X0/D0 = {ratio:.2f}:")
    print(f"  X0 = {X0_alt:.4f} m ({X0_alt/0.3048:.1f} ft)")
    print(f"  D0 = {D0_alt:.4f} m (2.0 ft)")
    print(f"\n  T0 = {T0_alt:.2f}°C (esperado: 26.96°C, error: {abs(T0_alt-26.96):.2f}°C)")
    print(f"  Tj = {Tj_alt:.2f}°C (esperado: 35.25°C, error: {abs(Tj_alt-35.25):.2f}°C)")
    print(f"  Pj = {Pj_alt:.2f} mm Hg (esperado: 26.87 mm Hg, error: {abs(Pj_alt-26.87):.2f} mm Hg)")
