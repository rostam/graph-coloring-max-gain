import glob, os, pandas as pd, numpy as np
os.chdir(".")
for file in glob.glob("*.csv"):
    df = pd.read_csv(file)
    print(file)
    print(sum(np.sign(df['MaxGain_nat'] - df['ignore_nat'])))