import pandas as pd, matplotlib.pyplot as plt, numpy as np, scipy.io as sio

plt.rcParams.update({'font.size': 18})


def vis(mat_name):
    mm = sio.mmread(mat_name)
    size = mm.shape
    nnz = mm.nnz
    data = pd.read_csv(mat_name + '0.csv')
    data['mnat'] = (data['mnat'])
    data['mnew'] = (data['mnew'])
    data['mlfo'] = (data['mlfo'])
    data['msat'] = (data['msat'])
    ax = data.plot.bar(x='numOfColor',legend=None,
                       # title='Number of discovered nonzeros based on given numbers of '
                       #       'colors \n for matrix ' + mat_name + ' with ' + str(nnz) + ' nonzeros',
                       # yticks=np.arange(0, 8402, 500),
                       # xticks=np.arange(8, 18, 1),
                       ylim=(0, nnz))
    ax.set_xlabel("Number of colors")
    ax.set_ylabel("Number of discovered nonzeros")
    # data_nos3.plot(yticks(np.arange(0, 8402, 100))
    # ax.legend(["NO", "AGO", "LFO", "IDO"])
    plt.savefig('/home/rostam/Desktop/projs/papers/MaximizeNonzerosUncertain/PaperSources/nos3_new.png',bbox_inches='tight')

mat_name = 'mats/nos3.mtx'
vis(mat_name)


