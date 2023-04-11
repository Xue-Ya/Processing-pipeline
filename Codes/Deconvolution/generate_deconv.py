import pandas as pd

file_path = ''
output_path = ''
file_path = file_path + 'raw_meta.csv'

meta= pd.read_csv(file_path,index_col=0,encoding = "ISO-8859-1")

meta = meta[:,6:]

meta.to_csv()