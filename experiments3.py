import pandas as pd, matplotlib.pyplot as plt, numpy as np, scipy.io as sio

plt.rcParams.update({'font.size': 18})
data = pd.read_csv('G51.mtxk.csv')
data.plot(x='k')
plt.show()