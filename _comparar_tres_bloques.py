import pandas as pd
from app import solve_case

raw = pd.read_excel('Resultados.xlsx', sheet_name='Hoja1', header=None)

# bloques: (nombre, col_inicio) con esquema [label, valor, unidad,...]
# base: cols 1..4 ; icl_mod: 6..9 ; m_mod: 11..14
blocks = [
    ('base', 1),
    ('icl_mod', 6),
    ('m_mod', 11),
]

sections = [
    (19, 2.0),
    (25, 1.5),
    (31, 1.0),
    (37, 0.5),
]

all_rows = []

for bname, c in blocks:
    TA = float(raw.iloc[2, c+1])
    RH_A = float(raw.iloc[3, c+1]) / 100.0
    Tmr = float(raw.iloc[4, c+1])
    Icl = float(raw.iloc[6, c+1])
    M = float(raw.iloc[7, c+1])
    D0 = float(raw.iloc[9, c+1])
    X0 = float(raw.iloc[10, c+1])

    for start, vj in sections:
        for r, rh0 in zip([start+2, start+3, start+4], [0.9, 0.95, 1.0]):
            # some cells are empty in right block
            t0_ref = raw.iloc[r, c+1]
            tj_ref = raw.iloc[r, c+2]
            pj_ref = raw.iloc[r, c+3]
            if pd.isna(t0_ref) or pd.isna(tj_ref) or pd.isna(pj_ref):
                continue

            out = solve_case(TA, RH_A, Tmr, vj, M, Icl, D0, X0, rh0=rh0, p_atm_kpa=101.325, include_buoyancy=False)
            t0 = float(out['T0'])
            tj = float(out['Tj'])
            pj = float(out['Pj'])

            all_rows.append({
                'bloque': bname,
                'Vj': vj,
                'rh0': rh0,
                'TA': TA,
                'RH_A': RH_A,
                'Tmr': Tmr,
                'Icl': Icl,
                'M': M,
                'D0': D0,
                'X0': X0,
                'T0_ref': float(t0_ref), 'T0_calc': t0, 'dT0': t0 - float(t0_ref),
                'Tj_ref': float(tj_ref), 'Tj_calc': tj, 'dTj': tj - float(tj_ref),
                'Pj_ref': float(pj_ref), 'Pj_calc': pj, 'dPj': pj - float(pj_ref),
            })

cmp = pd.DataFrame(all_rows)
cmp.to_csv('comparacion_resultados_tres_bloques.csv', index=False)

summary = cmp.groupby('bloque').agg(
    filas=('bloque','size'),
    max_abs_dT0=('dT0', lambda s: float(s.abs().max())),
    max_abs_dTj=('dTj', lambda s: float(s.abs().max())),
    max_abs_dPj=('dPj', lambda s: float(s.abs().max())),
    mean_abs_dT0=('dT0', lambda s: float(s.abs().mean())),
    mean_abs_dTj=('dTj', lambda s: float(s.abs().mean())),
    mean_abs_dPj=('dPj', lambda s: float(s.abs().mean())),
).reset_index()

summary.to_csv('comparacion_resultados_tres_bloques_resumen.csv', index=False)
print(summary.to_string(index=False))
print('\nSaved: comparacion_resultados_tres_bloques.csv')
print('Saved: comparacion_resultados_tres_bloques_resumen.csv')
