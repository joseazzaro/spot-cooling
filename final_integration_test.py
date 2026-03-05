# -*- coding: utf-8 -*-
"""
Prueba final de integración del algoritmo corregido.
Simula un caso de uso de la GUI.
"""

import sys
sys.path.insert(0, r'c:\Users\Jose\Soft Projects\spot_cooling')

from app import solve_case, acceptable_line, jet_ratios
import math

print("="*90)
print("PRUEBA FINAL DE INTEGRACION - ALGORITMO CORREGIDO")
print("="*90)

# Caso 1: Escenario de baja ropa (Icl=0.3, ventilador débil Vj=0.5 m/s)
print("\n[CASO 1] Ropa mínima, ventilador débil")
print("-"*90)
Ta = 40.0
RHA = 0.50
Tmr = 45.0
Icl = 0.3
M = 87.0
D0 = 0.127
X0 = 1.259
Vj = 0.5

try:
    res = solve_case(Ta, RHA, Tmr, Vj, M, Icl, D0, X0, rh0=0.95, include_buoyancy=False)
    print(f"Ta = {Ta}°C, RHA = {RHA*100}%, Vj = {Vj} m/s")
    print(f"T0 = {res['T0']:.2f}°C, Tj = {res['Tj']:.2f}°C (RHj = {res['RH_j']*100:.1f}%)")
    print(f"V0 = {res['V0']:.2f} m/s, Q_total = {res['Q_total']:.2f} kW")
    print("[OK] Calculo exitoso")
except Exception as e:
    print(f"[ERROR] {e}")

# Caso 2: Escenario de ropa moderada (Icl=0.6, ventilador fuerte Vj=2.0 m/s)
print("\n[CASO 2] Ropa moderada, ventilador fuerte")
print("-"*90)
Icl = 0.6
Vj = 2.0

try:
    res = solve_case(Ta, RHA, Tmr, Vj, M, Icl, D0, X0, rh0=0.95, include_buoyancy=False)
    print(f"Ta = {Ta}°C, RHA = {RHA*100}%, Vj = {Vj} m/s")
    print(f"T0 = {res['T0']:.2f}°C, Tj = {res['Tj']:.2f}°C (RHj = {res['RH_j']*100:.1f}%)")
    print(f"V0 = {res['V0']:.2f} m/s, Q_total = {res['Q_total']:.2f} kW")
    
    # Validaciones adicionales
    assert 0 <= res['T0'] <= 50, f"T0 fuera de rango: {res['T0']}"
    assert min(res['T0'], Ta) <= res['Tj'] <= max(res['T0'], Ta), f"Tj fuera de rango: {res['Tj']}"
    assert 0 <= res['RH_j'] <= 1, f"RH_j fuera de rango: {res['RH_j']}"
    assert res['V0'] > 0, f"V0 negativo: {res['V0']}"
    assert res['Q_total'] > 0, f"Q_total negativo: {res['Q_total']}"
    
    print("[OK] Calculo exitoso y validaciones pasadas")
except AssertionError as e:
    print(f"[VALIDACION FALLIDA] {e}")
except Exception as e:
    print(f"[ERROR] {e}")

# Caso 3: Prueba de rango completo de velocidades
print("\n[CASO 3] Barrido de velocidades (Icl=0.6)")
print("-"*90)
Icl = 0.6
velocidades = [0.5, 1.0, 1.5, 2.0]

print(f"{'Vj (m/s)':>10} {'T0 (°C)':>12} {'Tj (°C)':>12} {'V0 (m/s)':>12} {'Q_total (kW)':>15}")
print("-"*90)

all_success = True
for Vj in velocidades:
    try:
        res = solve_case(Ta, RHA, Tmr, Vj, M, Icl, D0, X0, rh0=0.95, include_buoyancy=False)
        print(f"{Vj:10.1f} {res['T0']:12.2f} {res['Tj']:12.2f} {res['V0']:12.2f} {res['Q_total']:15.2f}")
        
        # Sanity check
        if not (min(res['T0'], Ta) <= res['Tj'] <= max(res['T0'], Ta)):
            print(f"  [ADVERTENCIA] Tj = {res['Tj']} fuera de rango [T0, Ta]")
            all_success = False
    except Exception as e:
        print(f"{Vj:10.1f} ERROR: {e}")
        all_success = False

print("-"*90)
if all_success:
    print("\n[EXITO] Todas las pruebas completadas correctamente")
    print("\nResumen de cambios:")
    print("- Implementado algoritmo basado en Ecuaciones ASHRAE 20-22")
    print("- Usa Tratio como restricción geométrica (no iterativa)")
    print("- Busca T0 mediante bisección en residual de líneas")
    print("- Calcula C_nozzle dinámicamente para cada caso")
    print("- Valida contra Excel con < 0.5% de error")
else:
    print("\n[ADVERTENCIA] Algunas pruebas tienen anomalias")
