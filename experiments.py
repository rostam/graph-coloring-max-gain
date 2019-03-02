import pandas as pd, matplotlib.pyplot as plt, numpy as np
plt.rcParams.update({'font.size': 16})

# nnz = 8402
# data_nos3 = pd.read_csv('nos3.mtx.csv')
# data_nos3['mnat'] = 8402 - (data_nos3['mnat']/2)
# data_nos3['mnew'] = 8402 - (data_nos3['mnew']/2)
# data_nos3['mlfo'] = 8402 - (data_nos3['mlfo']/2)
# data_nos3['msat'] = 8402 - (data_nos3['msat']/2)
# ax = data_nos3.plot.bar(x='numOfColor',
#                         title='Number of discovered nonzeros based on given numbers of '
#                               'colors for matrix nos3.mtx with 8402 nonzeros',
#                         yticks=np.arange(0, 8402, 500),
#                         xticks=np.arange(8, 18, 1))
# ax.set_xlabel("Number of colors")
# ax.set_ylabel("Number of discovered nonzeros")
# # data_nos3.plot(yticks(np.arange(0, 8402, 100))
# ax.legend(["Natural ordering", "New ordering", "LFO ordering", "SAT ordering"])
#
#
# nnz=7017
# data_vcss08 = pd.read_csv('results_bcsstk08.mtx_329_339.csv')
# data_vcss08['mnat'] = nnz - (data_vcss08['mnat'])
# data_vcss08['mnew'] = nnz - (data_vcss08['mnew'])
# data_vcss08['mlfo'] = nnz - (data_vcss08['mlfo'])
# data_vcss08['msat'] = nnz - (data_vcss08['msat'])
# ax = data_vcss08.plot.bar(x='numOfColor',
#                       title='Number of discovered nonzeros based on given numbers '
#                             'of colors for matrix nos3.mtx with 7017 nonzeros',
#                       yticks=np.arange(0, 7017, 500),
#                       xticks=np.arange(329, 339, 1))
# ax.set_xlabel("Number of colors")
# ax.set_ylabel("Number of discovered nonzeros")
# ax.legend(["Natural ordering", "New ordering", "LFO ordering", "SAT ordering"])

nnz = 15963
data_plbuckle = pd.read_csv('results_plbuckle.mtx_22_52.csv')
data_plbuckle['mnat'] = nnz - (data_plbuckle['mnat'])
data_plbuckle['mnew'] = nnz - (data_plbuckle['mnew'])
data_plbuckle['mlfo'] = nnz - (data_plbuckle['mlfo'])
data_plbuckle['msat'] = nnz - (data_plbuckle['msat'])
ax = data_plbuckle.plot.bar(x='numOfColor',
                            title='Number of discovered nonzeros based on given numbers '
                                  'of colors for matrix plbuckle.mtx with 15963 nonzeros',
                            # yticks=np.arange(0, nnz, 500),
                            xticks=np.arange(22, 52, 1))
ax.set_xlabel("Number of colors")
ax.set_ylabel("Number of discovered nonzeros")
ax.legend(["Natural ordering", "MAT ordering", "LFO ordering", "SAT ordering"])
#
nnz = 5909
data_G51 = pd.read_csv('results_G51.mtx_147_157.csv')
data_G51['mnat'] = nnz - (data_G51['mnat'])
data_G51['mnew'] = nnz - (data_G51['mnew'])
data_G51['mlfo'] = nnz - (data_G51['mlfo'])
data_G51['msat'] = nnz - (data_G51['msat'])
ax = data_plbuckle.plot(x='numOfColor', title='The gain based on given numbers of colors for matrix G51.mtx with'
                                              '5909 nonzeros')
ax.set_xlabel("Number of colors")
ax.set_ylabel("Number of discovered nonzeros")
ax.legend(["Natural ordering", "MAT ordering", "LFO ordering", "SAT ordering"])
plt.show()

# data = data_nos3.merge(data_vcss08, on='key').merge(data_plbuckle, on='key').merge(data_G51,  on='key')
# print(data.head())

