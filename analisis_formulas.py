# -*- coding: utf-8 -*-
"""
Analisis de formulas del proyecto Spot Cooling - ASHRAE (Azer)
Este script documenta y verifica todas las ecuaciones implementadas
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")
from coeffs import regress_Ta50_coeffs, TABLE
import math

print("="*80)
print("ANALISIS DE FORMULAS: SPOT COOLING ASHRAE (AZER)")
print("="*80)

print("\n1. TABLA DE COEFICIENTES DE REGRESION")
print("-" * 80)
print("De: Tabla 1 (Azer, ASHRAE Part 1)")
print("Formula: Ta(0.5) = a * Tmr + b")
print("\nDonde:")
print("  - Ta(0.5): Temperatura aceptable del aire al 50 porciento de personas")
print("  - Tmr: Temperatura media radiante (C)")
print("  - a, b: Coeficientes que dependen de Icl, M, V")
print("\nParametros:")
print("  - Icl: Aislamiento de ropa (clo)")
print("  - M: Metabolismo (W/m2)")
print("  - V: Velocidad del aire en la zona de ocupacion (m/s)")

print("\n\nTABLA DE COEFICIENTES:")
print("-" * 80)
for icl in [0.3, 0.6, 0.9]:
    print("\n  Icl =", icl, "clo:")
    for m in [87.0, 174.0, 232.0]:
        print("    M =", m, "W/m2:")
        for v in sorted(TABLE[icl][m].keys()):
            a, b = TABLE[icl][m][v]
            print(f"      V = {v} m/s:  a = {a:7.3f},  b = {b:6.1f}")

print("\n\n2. COEFICIENTES DE CONVECCION Y RADIACION")
print("-" * 80)
print("hc = 8.3 * V_j^0.6")
print("Donde: hc = coeficiente de conveccion (W/m2K)")
print("\nhr = 3.87 + 0.031 * Tmr")
print("Donde: hr = coeficiente de radiacion (W/m2K)")
print("\nh = hc + hr (coeficiente total)")

print("\n\n3. FACTORES DE ROPA Y EFECTIVIDAD")
print("-" * 80)
print("fcl = 1.0 + 0.2 * Icl")
print("Donde: fcl = factor de area de ropa")
print("\nFcl = 1.0 / (1.0 + 0.155 * fcl * h * Icl)")
print("Donde: Fcl = factor de efectividad de la ropa")
print("\nFpcl = 1.0 / (1.0 + 0.143 * hc * Icl)")
print("Donde: Fpcl = factor de efectividad de la ropa individual")

print("\n\n4. LINEA ACEPTABLE (ACCEPTABLE LINE)")
print("-" * 80)
print("W = 0.5  (parametro de zona de confort)")
print("m = (h * fcl * Fcl) / (2.2 * hc * Fpcl * W)")
print("Ta(0.5) = a * Tmr + b")
print("P_s(0.5) = presion de vapor saturado a Ta(0.5)")
print("C = m * Ta(0.5) + 0.5 * P_s(0.5)")
print("\nEcuacion de la linea aceptable:")
print("P_v = m * Ta + C")

print("\n\n5. RATIOS DE JET (RELACIONES DE VELOCIDAD Y TEMPERATURA)")
print("-" * 80)
print("r = (X0/D0) + 2.572")
print("\nSin buoyancy:")
print("  Vratio = 1.464 / r")
print("  Tratio = 4.539 / r")
print("\nCon buoyancy (Numero de Arquimedes):")
print("  Ar = g * beta * dT * D0 / V0^2")

print("\n\n" + "="*80)
print("EJEMPLO DE CALCULO PASO A PASO")
print("="*80)

# Ejemplo numérico
print("\nParametros de entrada:")
TA = 40.0  # C
RH_A = 0.50  # 50 porciento
Tmr = 45.0  # C
Vj = 2.0  # m/s
M = 87.0  # W/m2
Icl = 0.6  # clo
D0 = 0.127  # m
X0 = 1.259  # m

print(f"  TA = {TA} C (temperatura ambiente)")
print(f"  RH_A = {RH_A*100:.1f}% (humedad relativa ambiente)")
print(f"  Tmr = {Tmr} C (temperatura media radiante)")
print(f"  Vj = {Vj} m/s (velocidad del jet)")
print(f"  M = {M} W/m2 (metabolismo)")
print(f"  Icl = {Icl} clo (aislamiento ropa)")
print(f"  D0 = {D0} m (diametro abertura)")
print(f"  X0 = {X0} m (distancia)")

print("\nPaso 1: Coeficientes de conveccion y radiacion")
hc = 8.3 * (max(1e-6, Vj)**0.6)
hr = 3.87 + 0.031 * Tmr
h = hc + hr
print(f"  hc = 8.3 * {Vj}^0.6 = {hc:.3f} W/(m2K)")
print(f"  hr = 3.87 + 0.031 * {Tmr} = {hr:.3f} W/(m2K)")
print(f"  h = hc + hr = {h:.3f} W/(m2K)")

print("\nPaso 2: Factores de ropa")
fcl = 1.0 + 0.2 * Icl
Fcl = 1.0 / (1.0 + 0.155 * fcl * h * Icl)
Fpcl = 1.0 / (1.0 + 0.143 * hc * Icl)
print(f"  fcl = 1 + 0.2 * {Icl} = {fcl:.3f}")
print(f"  Fcl = 1 / (1 + 0.155 * {fcl:.3f} * {h:.3f} * {Icl}) = {Fcl:.4f}")
print(f"  Fpcl = 1 / (1 + 0.143 * {hc:.3f} * {Icl}) = {Fpcl:.4f}")

print("\nPaso 3: Coeficientes de regresion Ta(0.5)")
a, b = regress_Ta50_coeffs(Icl, M, Vj)
print(f"  Regresion para Icl={Icl}, M={M}, V={Vj}:")
print(f"  a = {a:.4f},  b = {b:.2f}")

print("\nPaso 4: Linea aceptable")
Ta50 = a * Tmr + b
Psat_Ta50 = math.exp(18.6686 - 4030.183/(Ta50 + 235.0))
W = 0.5
m = (h * fcl * Fcl) / (2.2 * hc * Fpcl * W)
C = m * Ta50 + 0.5 * Psat_Ta50
print(f"  Ta(0.5) = {a:.4f} * {Tmr} + {b:.2f} = {Ta50:.2f} C")
print(f"  Ps(Ta(0.5)) = {Psat_Ta50:.3f} mmHg")
print(f"  W = {W}")
print(f"  m = ({h:.3f} * {fcl:.3f} * {Fcl:.4f}) / (2.2 * {hc:.3f} * {Fpcl:.4f} * {W})")
print(f"    = {m:.6f} mmHg/C")
print(f"  C = {m:.6f} * {Ta50:.2f} + 0.5 * {Psat_Ta50:.3f} = {C:.3f} mmHg")
print(f"\n  LINEA ACEPTABLE: Pv = {m:.6f} * Ta + {C:.3f}")

print("\nPaso 5: Ratios de jet")
r = (X0 / D0) + 2.572
Vratio = 1.464 / r
Tratio = 4.539 / r
print(f"  r = ({X0}/{D0}) + 2.572 = {r:.4f}")
print(f"  Vratio = 1.464 / {r:.4f} = {Vratio:.4f}")
print(f"  Tratio = 4.539 / {r:.4f} = {Tratio:.4f}")

print("\n" + "="*80)
print("RESUMEN DE ECUACIONES IMPLEMENTADAS")
print("="*80)
print("""
Las siguientes ecuaciones estan implementadas en el codigo:

1. Psicrometria:
   - P_sat(T) = exp(18.6686 - 4030.183/(T + 235.0))  [mmHg]
   - RH = P_v / P_sat(T)
   - humidity_ratio = 0.62198 * P_v[kPa] / (P_atm[kPa] - P_v[kPa])
   - enthalpia = 1.006*T + w*(2501.0 + 1.86*T)  [kJ/kg DA]

2. Dinamica del jet:
   - Velocidad: Vratio = (1.464/r) / (1 + 0.21*Ar*r^2)^(1/3)
   - Temperatura: Tratio = (4.539/r) / (1 + 0.21*Ar*r^2)^(1/3)
   - Numero de Arquimedes: Ar = g*beta*dT*D0/V0^2

3. Termofisiologia (confort termico):
   - Conveccion: hc = 8.3 * V_j^0.6
   - Radiacion: hr = 3.87 + 0.031 * Tmr
   - Factor de ropa: fcl = 1 + 0.2*Icl
   - Efectividad: Fcl = 1 / (1 + 0.155*fcl*h*Icl)
   - Pendiente de linea aceptable: m = (h*fcl*Fcl) / (2.2*hc*Fpcl*W)

4. Regresion de temperatura aceptable:
   - Ta(0.5) = a*Tmr + b (coeficientes de Tabla 1, Azer ASHRAE)
   - Interpolacion trilineal en (Icl, M, V)
""")
