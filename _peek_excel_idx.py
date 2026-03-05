import pandas as pd
raw = pd.read_excel('Resultados.xlsx', sheet_name='Hoja1', header=None)
for i in range(0,15):
    vals = [raw.iloc[i,j] for j in range(0,6)]
    print(i, vals)
