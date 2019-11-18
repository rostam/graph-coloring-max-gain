import glob, os, pandas as pd, numpy as np
df = pd.read_csv("G51_res.csv")
print(df['MaxGain_lfo'])
# print(df[df['MaxGain_lfo']>df['ignore_nat']].sum())