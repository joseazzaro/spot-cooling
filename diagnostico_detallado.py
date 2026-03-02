# -*- coding: utf-8 -*-
"""
Diagnóstico detallado: Verificar si la bisección está encontrando el punto correcto
"""

import sys
import math
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import (solve_case, psat_mmhg, jet_ratios, acceptable_line, 
                 pv_from_rh_T, ta50_from_regression)

print("="*80)
print("DIAGNÓSTICO DETALLADO: RH0 = 0.95")
print("="*80)

# Parámetros exactos del ejemplo
TA = 40.0
RH_A = 0.50
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6
D0 = 0.3048  
X0 = 0.3048  
rh0 = 0.95

print(f"\nParámetros de entrada:")
print(f"  TA = {TA}°C, RH_A = {RH_A*100}%")
print(f"  Tmr = {Tmr}°C, Vj = {Vj} m/s")
print(f"  M = {M} W/m², Icl = {Icl} clo")
print(f"  D0 = {D0} m, X0 = {X0} m")
print(f"  RH0 = {rh0*100}%")

# Dato calculado
PA = pv_from_rh_T(RH_A, TA)
m, C, Ta50 = acceptable_line(Icl, M, Vj, Tmr)

print(f"\nParámetros calculados:")
print(f"  PA = {PA:.2f} mm Hg")
print(f"  m = {m:.6f} mm Hg/°C")
print(f"  C = {C:.2f} mm Hg")
print(f"  Ta(0.5) = {Ta50:.2f}°C")

# Línea aceptable: P = -m*T + C
print(f"\nLínea aceptable: P = -{m:.4f}*T + {C:.2f}")
print(f"             P = -m*T + C")

# Resultados esperados de TABLA 2
print(f"\n" + "="*80)
print("RESULTADOS ESPERADOS (TABLA 2):")
print("="*80)
print(f"  T0 = 26.96°C")
print(f"  Tj = 35.25°C")
print(f"  Pj = 26.87 mm Hg")
print(f"  Verificación: Pj = -m*Tj + C = -{m:.4f}*35.25 + {C:.2f} = {-m*35.25 + C:.2f} ✓")

# Ver valores de residual en el rango
print(f"\n" + "="*80)
print("ANÁLISIS DEL RESIDUAL A DIFERENTES T0:")
print("="*80)

p0_test = rh0 * psat_mmhg(26.96)  # T0 esperado
print(f"\nP0 en T0=26.96 (esperado): P0 = {rh0}*psat({26.96}) = {p0_test:.2f} mmHg")

for T0_test in [19.40, 26.96, 30.0, 40.0]:
    P0 = rh0 * psat_mmhg(T0_test)
    Vratio, Tratio = jet_ratios(X0, D0, False, TA, T0_test)
    Ti = TA - Tratio * (TA - T0_test)
    Pi = PA - Tratio * (PA - P0)
    
    # Línea aceptable: P = -m*T + C
    P_on_line = -m * Ti + C
    residual = Pi - P_on_line
    
    print(f"\nT0 = {T0_test}°C:")
    print(f"  P0 = {P0:.2f} mm Hg")
    print(f"  Vratio = {Vratio:.6f}, Tratio = {Tratio:.6f}")
    print(f"  Ti = {Ti:.2f}°C, Pi = {Pi:.2f} mm Hg")
    print(f"  P_on_line (=-m*Ti+C) = {P_on_line:.2f} mm Hg")
    print(f"  Residual = Pi - P_on_line = {residual:.2f} mm Hg")

# Ahora ejecutar el solver
print(f"\n" + "="*80)
print("RESULTADO DEL SOLVER (solve_case):")
print("="*80)

result = solve_case(TA, RH_A, Tmr, Vj, M, Icl, D0, X0, rh0=rh0)
T0 = result['T0']
Tj = result['Tj']
Pj = result['Pj']
P0 = result['P0']

print(f"\n  T0 = {T0:.2f}°C (esperado: 26.96°C, error: {abs(T0-26.96):.2f}°C)")
print(f"  Tj = {Tj:.2f}°C (esperado: 35.25°C, error: {abs(Tj-35.25):.2f}°C)")
print(f"  Pj = {Pj:.2f} mm Hg (esperado: 26.87 mm Hg, error: {abs(Pj-26.87):.2f} mm Hg)")
print(f"  P0 = {P0:.2f} mm Hg")

# Verificar si el punto está en la línea aceptable
P_on_line_at_calc = -m * Tj + C
print(f"\nVerificación de punto en línea:")
print(f"  P_on_line (=-m*Tj+C) = {P_on_line_at_calc:.2f} mm Hg")
print(f"  Pj (calculado) = {Pj:.2f} mm Hg")
print(f"  Diferencia = {abs(Pj - P_on_line_at_calc):.4f} mm Hg")

# Chequear si hay problema con X0/D0
print(f"\n" + "="*80)
print("ANÁLISIS DE X0 Y D0:")
print("="*80)
print(f"  X0/D0 = {X0/D0} (1 ft / 1 ft = 1.0)")
r = (X0/D0) + 2.572
print(f"  r = X0/D0 + 2.572 = {r}")
Vratio_example = 1.464/r
Tratio_example = 4.539/r
print(f"  Vratio = 1.464/r = {Vratio_example:.6f}")
print(f"  Tratio = 4.539/r = {Tratio_example:.6f}")

print(f"\n¿Será que X0 ≠ D0 en el ejemplo ASHRAE?")
