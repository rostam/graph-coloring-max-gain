import glob, os, pandas as pd, numpy as np
df = pd.read_csv("results.csv")
print(df[df['MaxGain_nat']>df['ignore_nat']].sum())