import pandas as pd
pd.set_option('display.width', 240)
pd.set_option('display.max_columns', 50)
df = pd.read_excel('Resultados.xlsx', sheet_name='Hoja1')
print('shape:', df.shape)
print('columns:', list(df.columns))
print('\nHEAD:')
print(df.head(40).to_string(index=False))
