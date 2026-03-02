# -*- coding: utf-8 -*-
"""
VALIDACIÓN FINAL: Tabla 2 ASHRAE CORRECTA
Usando ecuaciones 25 y 27 (intersección de P0 psicométrico vs fisiológico)
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import solve_case

print("\n" + "="*100)
print("VALIDACIÓN: TABLA 2 ASHRAE (Correcta) - Ejemplo 1")
print("="*100)

# Parámetros del ejemplo
TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.3048  # 1 ft
X0 = 3.048   # 10 ft

print(f"\nParámetros ASHRAE Example 1:")
print(f"  TA = {TA}°C, RH_A = {RH_A*100}%, Tmr = {Tmr}°C")
print(f"  Vj = {Vj} m/s, D0 = {D0} m (1 ft), X0 = {X0} m (10 ft)")
print(f"  M = {M} W/m², Icl = {Icl} clo")

# TABLA 2 ASHRAE CORRECTA (del documento)
tabla2_correct = {
    0.90: {"T0": 27.58, "Tj": 35.48, "Pj": 26.33},
    0.95: {"T0": 26.96, "Tj": 35.25, "Pj": 26.87},
    1.00: {"T0": 26.38, "Tj": 35.04, "Pj": 27.03},
}

print(f"\n" + "="*100)
print("COMPARACIÓN: Valores ASHRAE vs Calculados")
print("="*100)

print(f"\n{'RH₀':<8} {'Parámetro':<12} {'ASHRAE':<15} {'Calculado':<15} {'Error':<15} {'Status':<8}")
print("-" * 100)

all_passed = True
for rh0 in [0.90, 0.95, 1.00]:
    expected = tabla2_correct[rh0]
    result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
    
    T0_calc = result['T0']
    Tj_calc = result['Tj']
    Pj_calc = result['Pj']
    
    err_T0 = abs(T0_calc - expected['T0'])
    err_Tj = abs(Tj_calc - expected['Tj'])
    err_Pj = abs(Pj_calc - expected['Pj'])
    
    tolerance = 0.05  # ±0.05 for temperature, pressure
    status_T0 = "✓" if err_T0 < tolerance else "✗"
    status_Tj = "✓" if err_Tj < tolerance else "✗"
    status_Pj = "✓" if err_Pj < tolerance else "✗"
    
    if err_T0 >= tolerance or err_Tj >= tolerance or err_Pj >= tolerance:
        all_passed = False
    
    print(f"{rh0*100:<8.0f}°/₀ {'T0 (°C)':<12} {expected['T0']:<15.2f} {T0_calc:<15.2f} {err_T0:<15.4f} {status_T0:<8}")
    print(f"{'': <8} {'Tj (°C)':<12} {expected['Tj']:<15.2f} {Tj_calc:<15.2f} {err_Tj:<15.4f} {status_Tj:<8}")
    print(f"{'': <8} {'Pj (mmHg)':<12} {expected['Pj']:<15.2f} {Pj_calc:<15.2f} {err_Pj:<15.4f} {status_Pj:<8}")
    print()

print("="*100)
if all_passed:
    print("✓✓✓ VALIDACIÓN EXITOSA: ¡TABLA 2 ASHRAE COMPLETAMENTE VERIFICADA! ✓✓✓")
else:
    print("⚠ Algunos valores están fuera de tolerancia")
print("="*100)

# Detalles completos para RH0=0.95
print(f"\nDETALLES TÉCNICOS para RH₀ = 0.95 (referencia ASHRAE):")
print("="*100)
result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=0.95)

print(f"\nParámetros de la línea aceptable:")
print(f"  m = {result['m']:.6f} mm Hg/°C")
print(f"  C (target area, Eq 26) = {result['C']:.2f} mm Hg")
print(f"  C (nozzle, Eq 27) = 45.32 mm Hg")

print(f"\nCondiciones en el orificio (nozzle):")
print(f"  T0 = {result['T0']:.2f}°C")
print(f"  P0 = {result['P0']:.2f} mm Hg")
print(f"  RH0 = {result['RH_0']*100:.1f}%")

print(f"\nCondiciones en zona de impacto (target area):")
print(f"  Tj = {result['Tj']:.2f}°C")
print(f"  Pj = {result['Pj']:.2f} mm Hg")
print(f"  RHj = {result['RH_j']*100:.1f}%")

print(f"\nParámetros del jet:")
print(f"  Vratio = {result['Vratio']:.6f}")
print(f"  Tratio = {result['Tratio']:.6f}")
print(f"  V0 = {result['V0']:.2f} m/s")

print(f"\nBalance térmico:")
print(f"  m_dot0 = {result['m_dot0']:.4f} kg/s")
print(f"  Q_total = {result['Q_total']:.2f} W")
print(f"  Q_sens = {result['Q_sens']:.2f} W")

# Verificación de ecuaciones
print(f"\n" + "="*100)
print("VERIFICACIÓN DE ECUACIONES ASHRAE para RH₀ = 0.95:")
print("="*100)

import math
T0 = result['T0']
P0_result = result['P0']

# Ec 25: P0 = rh0 * exp[18.6686 - 4030.183/(T0+235)]
P0_eq25 = 0.95 * math.exp(18.6686 - 4030.183/(T0 + 235))
print(f"\nEcuación 25 (Psicométrica):")
print(f"  P0 = 0.95 × Exp[18.6686 - 4030.183/(T0+235)]")
print(f"  P0 = 0.95 × Exp[18.6686 - 4030.183/{T0+235:.1f}]")
print(f"  P0 = {P0_eq25:.2f} mm Hg (calculado)")
print(f"  P0 = {P0_result:.2f} mm Hg (del solver)")
print(f"  Diferencia: {abs(P0_eq25 - P0_result):.4f} mm Hg")

# Ec 27: P0 = -m*T0 + 45.32
m = result['m']
P0_eq27 = -m * T0 + 45.32
print(f"\nEcuación 27 (Fisiológica):")
print(f"  P0 = -m × T0 + 45.32")
print(f"  P0 = -{m:.6f} × {T0:.2f} + 45.32")
print(f"  P0 = {P0_eq27:.2f} mm Hg (calculado)")
print(f"  P0 = {P0_result:.2f} mm Hg (del solver)")
print(f"  Diferencia: {abs(P0_eq27 - P0_result):.4f} mm Hg")

# Ec 26: Pj = -m*Tj + C
C_target = result['C']
Tj = result['Tj']
Pj_result = result['Pj']
Pj_eq26 = -m * Tj + C_target

print(f"\nEcuación 26 (Target Area):")
print(f"  Pj = -m × Tj + C")
print(f"  Pj = -{m:.6f} × {Tj:.2f} + {C_target:.2f}")
print(f"  Pj = {Pj_eq26:.2f} mm Hg (calculado)")
print(f"  Pj = {Pj_result:.2f} mm Hg (del solver)")
print(f"  Diferencia: {abs(Pj_eq26 - Pj_result):.4f} mm Hg")
