import glob, os, pandas as pd, numpy as np
os.chdir(".")
str = ""
for file in glob.glob("*.csv"):
    df = pd.read_csv(file)
    print(file)
    print(sum(np.sign(df['MaxGain_ago'] - df['ignore_ago'])))
