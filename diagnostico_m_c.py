# -*- coding: utf-8 -*-
"""
Diagnóstico: ¿Por qué m y C no coinciden con TABLA 1 ASHRAE?
"""

import sys
sys.path.insert(0, r"c:\Users\Jose\Soft Projects\spot_cooling")

from app import compute_factors, regress_Ta50_coeffs
import math

print("="*80)
print("DIAGNÓSTICO: CÁLCULO DE m Y C")
print("="*80)

# Parámetros
TA = 40.0
Tmr = 45.0
Vj = 2.0
M = 87.0
Icl = 0.6

print(f"\nParámetros:")
print(f"  TA = {TA}°C")
print(f"  Tmr = {Tmr}°C")
print(f"  Vj = {Vj} m/s")
print(f"  M = {M} W/m²")
print(f"  Icl = {Icl} clo")

# Paso 1: Calcular coeficientes de transferencia
hc, hr, h, fcl, Fcl, Fpcl = compute_factors(Vj, Icl, Tmr)

print(f"\nCoeficientes calculados:")
print(f"  hc = {hc:.3f} W/(m²K)")
print(f"  hr = {hr:.3f} W/(m²K)")
print(f"  h = {h:.3f} W/(m²K)")
print(f"  fcl = {fcl:.3f}")
print(f"  Fcl = {Fcl:.6f}")
print(f"  Fpcl = {Fpcl:.6f}")

# Paso 2: Calcular m usando la fórmula del código
print(f"\n" + "="*80)
print("FÓRMULA DEL CÓDIGO:")
print("="*80)

W = 0.5
m_nuestro = (h * fcl * Fcl) / (2.2 * hc * Fpcl * W)

print(f"\nm = (h × fcl × Fcl) / (2.2 × hc × Fpcl × W)")
print(f"W = {W}")
print(f"\nNumerador = {h:.3f} × {fcl:.3f} × {Fcl:.6f}")
print(f"          = {h*fcl*Fcl:.6f}")
print(f"\nDenominador = 2.2 × {hc:.3f} × {Fpcl:.6f} × {W}")
print(f"            = {2.2*hc*Fpcl*W:.6f}")
print(f"\nm_nuestro = {m_nuestro:.6f} mm Hg/°C")

# Paso 3: Comparar con TABLA 1 ASHRAE
m_ashrae = 0.737
C_ashrae = 52.85

print(f"\n" + "="*80)
print("COMPARACIÓN CON TABLA 1 ASHRAE:")
print("="*80)
print(f"\nm_ASHRAE = {m_ashrae} mm Hg/°C")
print(f"m_nuestro = {m_nuestro:.3f} mm Hg/°C")
print(f"Diferencia = {m_nuestro - m_ashrae:.3f} ({(m_nuestro/m_ashrae - 1)*100:.1f}% más alto)")

print(f"\nC_ASHRAE = {C_ashrae} mm Hg")

# Paso 4: Intentar encontrar la fórmula correcta
print(f"\n" + "="*80)
print("ANÁLISIS: ¿Cuál debería ser la fórmula?")
print("="*80)

# Si m_ASHRAE = 0.737 y sabemos los factores, ¿qué fórmula daría eso?
# m = (h * fcl * Fcl) / (2.2 * hc * Fpcl * W)
# 0.737 = (17.845 * 1.120 * 0.3498) / (2.2 * 12.580 * 0.4809 * W)
# 0.737 = 6.984 / (13.302 * W)
# Si W = 0.5: 0.737 = 6.984 / 6.651 = 1.0505 (NO, esto es lo que obtenemos)

# Pero si W es diferente... W = 6.984 / (0.737 * 13.302) = 0.714

# O si la fórmula es diferente, por ejemplo:
# m = fcl * Fcl / (2.2 * Fpcl) ?

m_test1 = (fcl * Fcl) / (2.2 * Fpcl)
print(f"\nPrueba 1: m = fcl × Fcl / (2.2 × Fpcl)")
print(f"  m = {fcl:.3f} × {Fcl:.6f} / (2.2 × {Fpcl:.6f})")
print(f"  m = {m_test1:.6f}  (NO coincide)")

# m = (h * fcl * Fcl) / hc ?
m_test2 = (h * fcl * Fcl) / hc
print(f"\nPrueba 2: m = h × fcl × Fcl / hc")
print(f"  m = {h:.3f} × {fcl:.3f} × {Fcl:.6f} / {hc:.3f}")
print(f"  m = {m_test2:.6f}  (NO coincide)")

# m = fcl * Fcl
m_test3 = fcl * Fcl
print(f"\nPrueba 3: m = fcl × Fcl")
print(f"  m = {fcl:.3f} × {Fcl:.6f}")
print(f"  m = {m_test3:.6f}  (NO coincide)")

# Espera... El ejemplo dice que m y C se basan en "Eqs 5, 6, 7, and 15 through 18"
# Estas ecuaciones podrían estar en el PDF. Sin verlas, no sé cuál es la fórmula correcta.

print(f"\n" + "="*80)
print("CONCLUSIÓN:")
print("="*80)
print("""
La fórmula de m en el código NO coincide con los valores de TABLA 1 ASHRAE.

Posibles causas:
1. La fórmula de m es incorrecta
2. W = 0.5 no es el valor correcto
3. Falta un factor de escala o una transformación
4. Las ecuaciones 5-7 y 15-18 del ASHRAE definen m de forma diferente

NECESARIO:
→ Obtener las ecuaciones 5-7, 15-18 del documento ASHRAE original
→ Reemplazar la fórmula de m y C en el código
""")

# Intentar retro-ingeniería
print(f"\n" + "="*80)
print("RETRO-INGENIERÍA:")
print("="*80)

# Si suponemos que C se calcula igual:
# C = m × Ta(0.5) + 0.5 × Psat(Ta(0.5))
# ASHRAE dice para Vj=2.0: C = 52.85

a, b = regress_Ta50_coeffs(Icl, M, Vj)
Ta50 = a * Tmr + b
Psat_Ta50 = math.exp(18.6686 - 4030.183/(Ta50 + 235.0))

print(f"\nTa(0.5) = {Ta50:.2f}°C")
print(f"Psat(Ta(0.5)) = {Psat_Ta50:.2f} mm Hg")

# C = m × Ta(0.5) + 0.5 × Psat(Ta(0.5))
# 52.85 = m × 37.99 + 0.5 × 49.67
# 52.85 = m × 37.99 + 24.835
# 28.015 = m × 37.99
# m = 0.737 ✓ ¡COINCIDE!

print(f"\nSi usamos C_ASHRAE = 52.85:")
print(f"  52.85 = m × {Ta50:.2f} + 0.5 × {Psat_Ta50:.2f}")
print(f"  52.85 = m × {Ta50:.2f} + {0.5*Psat_Ta50:.2f}")
print(f"  {52.85 - 0.5*Psat_Ta50:.2f} = m × {Ta50:.2f}")
print(f"  m = {(52.85 - 0.5*Psat_Ta50)/Ta50:.3f}")

print(f"\n¡ESTA ES LA SOLUCIÓN! Los valores de TABLA 1 ya asumen una m correcta.")
print(f"El problema es que nuestra fórmula para calcular m es incorrecta.")
