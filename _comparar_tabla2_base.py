import pandas as pd
from app import solve_case

raw = pd.read_excel('Resultados.xlsx', sheet_name='Hoja1', header=None)

TA = float(raw.iloc[2,2])
RH_A = float(raw.iloc[3,2]) / 100.0
Tmr = float(raw.iloc[4,2])
Icl = float(raw.iloc[6,2])
M = float(raw.iloc[7,2])
D0 = float(raw.iloc[9,2])
X0 = float(raw.iloc[10,2])

sections = [
    (19, 2.0),
    (25, 1.5),
    (31, 1.0),
    (37, 0.5),
]

rows = []
for start, vj in sections:
    for r, rh0 in zip([start+2, start+3, start+4], [0.9, 0.95, 1.0]):
        t0_ref = float(raw.iloc[r,2])
        tj_ref = float(raw.iloc[r,3])
        pj_ref = float(raw.iloc[r,4])

        out = solve_case(TA, RH_A, Tmr, vj, M, Icl, D0, X0, rh0=rh0, p_atm_kpa=101.325, include_buoyancy=False)
        t0 = float(out['T0'])
        tj = float(out['Tj'])
        pj = float(out['Pj'])

        rows.append({
            'Vj': vj, 'rh0': rh0,
            'T0_ref': t0_ref, 'T0_calc': t0, 'dT0': t0 - t0_ref,
            'Tj_ref': tj_ref, 'Tj_calc': tj, 'dTj': tj - tj_ref,
            'Pj_ref': pj_ref, 'Pj_calc': pj, 'dPj': pj - pj_ref,
        })

cmp = pd.DataFrame(rows)
cmp.to_csv('comparacion_tabla2_base.csv', index=False)
print('Base inputs:', {'TA':TA,'RH_A':RH_A,'Tmr':Tmr,'Icl':Icl,'M':M,'D0':D0,'X0':X0})
print('\nMax abs errors:')
print('dT0=', cmp['dT0'].abs().max(), 'dTj=', cmp['dTj'].abs().max(), 'dPj=', cmp['dPj'].abs().max())
print('\nMean abs errors:')
print('dT0=', cmp['dT0'].abs().mean(), 'dTj=', cmp['dTj'].abs().mean(), 'dPj=', cmp['dPj'].abs().mean())
print('\nAll rows:')
print(cmp.to_string(index=False))
print('\nSaved: comparacion_tabla2_base.csv')
