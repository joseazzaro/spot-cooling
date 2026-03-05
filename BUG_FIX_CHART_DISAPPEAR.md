# Bug Fix: Gráfico Desaparece al Ejecutar Cálculo

## Problema Identificado

El gráfico psicromético se dibujaba correctamente al ejecutar un cálculo ("Run & report"), pero luego desaparecía. El usuario tenía que volver a ejecutar el cálculo para que reapareciera.

## Causa Raíz

El problema ocurría en la interacción entre dos mecanismos:

1. **Timer de Debounce** (`_chart_update_timer`): Cuando el usuario ajusta valores en los spinboxes, hay un timer que espera 220ms después del último cambio antes de actualizar el gráfico (para evitar redibujados constantes).

2. **Lectura Inmediata de Valores**: Si el usuario estaba ajustando valores justo antes de hacer clic en "Run", el timer de debounce estaba pendiente. Después de que se dibujaba la solución, el timer se activaba y ejecutaba `_apply_chart_update()`.

3. **Limpieza Incondicional**: El método `_apply_chart_update()` siempre llamaba a `clear_operating_solution()` **después** de que la solución se dibujaba, borrándola del gráfico.

### Secuencia de Eventos (Buggy)

```
1. Usuario ajusta parámetro (ej: Vj = 1.5)
   → Timer inicia (espera 220ms)
   
2. Usuario hace clic en "Run"
   → Se calcula solución
   → Se dibuja gráfico con línea operativa
   → chart.fig.canvas.draw_idle() redibuja
   
3. Timer se dispara (220ms después)
   → _apply_chart_update() ejecuta
   → chart.clear_operating_solution() → operating_solution = None
   → chart.update_chart() redibuja
   → _plot_operating_line_and_intersections() ve que operating_solution es None
   → No dibuja la línea operativa
   → Gráfico aparece vacío/sin solución
```

## Solución Implementada

Se implementó una verificación inteligente que **solo limpia la solución si parámetros críticos cambiaron**:

### Parámetros que Invalidan la Solución

- `Icl` (Insulation de ropa): Cambia la línea aceptable
- `M` (Metabolic rate): Cambia la línea aceptable  
- `Tmr` (Mean radiant temperature): Cambia la línea aceptable
- `Vj` (Jet velocity): Parámetro de entrada del cálculo

### Parámetros que NO Invalidan la Solución

- `TA` (Temperatura ambiente): Cambia condiciones, pero solución aún es válida
- `RH_A` (Humedad relativa ambiente): Cambia condiciones, pero solución aún es válida
- `P_atm` (Presión atmosférica): Cambio menor, solución aún es aproximadamente válida
- `RH_0` (Humedad en tobera): Parámetro de entrada, no invalida solución anterior

### Código del Fix

```python
# En __init__: Track previous parameter values
self._prev_icl = self.e_ICL.value()
self._prev_m = self.e_M.value()
self._prev_tmr = self.e_TMR.value()
self._prev_vj = self.e_VJ.value()

# En _apply_chart_update(): Solo limpiar si parámetros críticos cambiaron
if icl_changed or m_changed or tmr_changed or vj_changed:
    self.chart.clear_operating_solution()
```

## Secuencia de Eventos (Fixed)

```
1. Usuario ajusta parámetro (ej: Vj = 1.5)
   → Timer inicia (espera 220ms)
   
2. Usuario hace clic en "Run"
   → Se calcula solución
   → Se dibuja gráfico con línea operativa
   → chart.fig.canvas.draw_idle() redibuja
   
3. Timer se dispara (220ms después)
   → _apply_chart_update() ejecuta
   → Verifica: ¿Icl, M, Tmr, Vj cambiaron?
   → SÍ: Vj cambió de 2.0 a 1.5
   → clear_operating_solution() se ejecuta
   → chart.update_chart() redibuja (vacío, esperando nuevo Run)
   
   ✅ PERO si usuario solo cambió TA o RH:
   → ¿Icl, M, Tmr, Vj cambiaron? NO
   → clear_operating_solution() NO se ejecuta
   → chart.update_chart() redibuja CON la solución operativa
```

## Impacto en UX

| Cambio de Parámetro | Comportamiento |
|---|---|
| Icl, M, Tmr | Solución se limpia (debe hacer Run nuevamente) |
| Vj | Solución se limpia (debe hacer Run nuevamente) |
| TA, RH_A, PATM | Solución se mantiene visible en gráfico |
| RH_0 (nozzle humidity) | Solución se mantiene - parámetro de entrada |
| Seleccionar punto en gráfico | Funciona independientemente |

## Testing

El fix se puede verificar con los siguientes pasos:

1. Abrir la aplicación
2. Ajustar algunos parámetros (ej: TA, RH)
3. Hacer clic en "Run & report"
4. El gráfico debe mostrar la línea operativa
5. Ajustar más parámetros ambientales (TA, RH)
6. El gráfico debe mantener la solución visible
7. Cambiar Icl o M
8. El gráfico debe limpiar la solución
9. Hacer clic en "Run & report" nuevamente
10. Nuevo resultado debe mostrarse

## Archivos Modificados

- `views/main_window.py`: Lógica de debounce y limpieza de solución

## Notas Técnicas

- El timer de debounce sigue siendo 220ms (no cambió)
- El único cambio es la **lógica de cuándo limpiar la solución**
- No afecta el método `update_chart()` de `psychro_chart.py`
- Completamente backward compatible
