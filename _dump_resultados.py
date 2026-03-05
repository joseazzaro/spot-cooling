import pandas as pd
import numpy as np

df = pd.read_excel('Resultados.xlsx', sheet_name='Hoja1', header=None)
out_lines = []
for i,row in df.iterrows():
    cells=[]
    for j,val in enumerate(row.tolist()):
        if pd.isna(val):
            continue
        txt=str(val).strip()
        if txt=='':
            continue
        cells.append(f'C{j}:{txt}')
    if cells:
        out_lines.append(f'R{i}: ' + ' | '.join(cells))

with open('Resultados_dump.txt','w',encoding='utf-8') as f:
    f.write('\n'.join(out_lines))
print('WROTE', len(out_lines), 'non-empty rows to Resultados_dump.txt')
