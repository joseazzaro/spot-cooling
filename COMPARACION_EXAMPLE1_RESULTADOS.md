# Comparacion Example 1 vs Excel (Resultados.xlsx)

## Alcance
- Archivo analizado: Resultados.xlsx (Hoja1).
- Tablas encontradas en el Excel: Tabla 1 y Tabla 2 (Vj=2.0, 1.5, 1.0, 0.5).
- No se encontraron Tabla 3, Tabla 4 ni Tabla 5 en este archivo.

## Metricas de error (Tabla 2, bloque base)
- Max |dT0| = 14.2776 C
- Max |dTj| = 9.1325 C
- Max |dPj| = 3.9225 mmHg
- Mean |dT0| = 5.0609 C
- Mean |dTj| = 3.4599 C
- Mean |dPj| = 1.4654 mmHg

## Observaciones
- Para Vj=2.0 m/s, la concordancia es buena (errores pequenos).
- Para Vj=1.5, 1.0 y 0.5 m/s hay diferencias grandes en T0/Tj/Pj respecto al Excel.
- El bloque base del Excel parece no seguir el mismo modelo implementado en app.py para esos Vj.

## Archivos generados
- comparacion_tabla2_base.csv (detalle fila a fila).
- Resultados_dump.txt (mapa de celdas no vacias del Excel).