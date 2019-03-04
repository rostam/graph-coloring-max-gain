import pandas as pd, matplotlib.pyplot as plt, numpy as np, scipy.io as sio

# plt.rcParams.update({'font.size': 18})
data = pd.read_csv('plbuckle.mtxk.csv')
ax = data.plot(x='k')
ax.set_xlabel("$k$")
ax.set_ylabel("Number of discovered nonzeros")
# data_nos3.plot(yticks(np.arange(0, 8402, 100))
ax.legend(["NO", "AGO", "LFO", "IDO"])
plt.savefig('/home/rostam/Desktop/projs/papers/MaximizeNonzerosUncertain/PaperSources/plbucklek2.png',bbox_inches='tight')
