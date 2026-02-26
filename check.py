import pandas as pd

# 1. Load the dataset (replace 'your_file.csv' with your actual file name)
df = pd.read_csv(r'C:\Users\LENOVO\bio sem 4\final_ppi_dataset.csv')

# 2. See headers and the first 5 rows
print(df.head())
print(list(df.columns))