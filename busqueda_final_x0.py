# -*- coding: utf-8 -*-
"""
Búsqueda final: X0 es distancia vertical. Buscar entre 9-10 ft
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import solve_case

print("="*80)
print("BÚSQUEDA FINAL: X0 es la altura VERTICAL del jet")
print("="*80)

TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.3048  # 1 ft (diámetro del orificio)
rh0 = 0.95

T0_expected = 26.96
Tj_expected = 35.25
Pj_expected = 26.87

print(f"\nD0 = 0.3048 m (1 ft) - Diámetro del orificio [fijo]")
print(f"\nBuscando X0 (altura vertical del jet) entre 9.0 y 10.5 ft:\n")
print(f"{'X0 (ft)':<10} {'X0 (m)':<12} {'T0':<10} {'Err_T0':<10} {'Tj':<10} {'Err_Tj':<10} {'Pj':<10} {'Err_Pj':<10}")
print("-" * 85)

best = None
best_error = float('inf')

for X0_ft_x100 in range(900, 1051):  # 9.00 a 10.50 ft en incrementos de 0.01 ft
    X0_ft = X0_ft_x100 / 100.0
    X0 = X0_ft * 0.3048
    
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
    
    T0 = result['T0']
    Tj = result['Tj']
    Pj = result['Pj']
    
    err_T0 = abs(T0 - T0_expected)
    err_Tj = abs(Tj - Tj_expected)
    err_Pj = abs(Pj - Pj_expected)
    
    # Error ponderado
    error_combined = err_T0 * 0.5 + err_Tj * 0.3 + err_Pj * 0.2
    
    # Mostrar cada 0.1 ft
    if X0_ft_x100 % 10 == 0:
        print(f"{X0_ft:<10.2f} {X0:<12.4f} {T0:<10.2f} {err_T0:<10.2f} {Tj:<10.2f} {err_Tj:<10.2f} {Pj:<10.2f} {err_Pj:<10.2f}")
    
    if error_combined < best_error:
        best_error = error_combined
        best = (X0_ft, X0, err_T0, err_Tj, err_Pj, T0, Tj, Pj)

print("\n" + "="*80)
print("MEJOR RESULTADO:")
print("="*80)
if best:
    X0_ft, X0, err_T0, err_Tj, err_Pj, T0, Tj, Pj = best
    print(f"\n✓ X0 = {X0_ft:.2f} ft = {X0:.4f} m")
    print(f"✓ D0 = 1.00 ft = 0.3048 m")
    print(f"\nResultados obtenidos:")
    print(f"  T0 = {T0:.2f}°C (esperado 26.96°C, error ±{err_T0:.4f}°C)")
    print(f"  Tj = {Tj:.2f}°C (esperado 35.25°C, error ±{err_Tj:.4f}°C)")
    print(f"  Pj = {Pj:.2f} mmHg (esperado 26.87 mmHg, error ±{err_Pj:.4f} mmHg)")
    
    # Probar valores exactos redondos
    print(f"\n" + "="*80)
    print("PRUEBA CON VALORES REDONDOS:")
    print("="*80)
    
    for X0_test in [9.0, 9.5, 9.75, 10.0, 10.25]:
        X0_m = X0_test * 0.3048
        result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0_m, rh0=rh0)
        print(f"\nX0 = {X0_test:.2f} ft ({X0_m:.4f} m):")
        print(f"  T0 ={result['T0']:.2f}°C,  Tj = {result['Tj']:.2f}°C,  Pj = {result['Pj']:.2f} mmHg")
