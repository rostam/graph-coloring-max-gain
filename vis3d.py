import pandas as pd, matplotlib.pyplot as plt, numpy as np
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({'font.size': 16})

df = pd.read_csv('mats/bcsstk08.mtx' + str(0) + '.csv')
df['k'] = 0
for i in range(1, 19):
    dff = pd.read_csv('mats/bcsstk08.mtx' + str(i) + '.csv')
    dff['k'] = i
    df = df.append(dff, ignore_index=True)

# print(df)
# f1 = pd.read_csv('mats/bcsstk08.mtx0.csv')
# f2 = pd.read_csv('mats/bcsstk08.mtx1.csv')
# f3 = pd.read_csv('mats/bcsstk08.mtx2.csv')
# f4 = pd.read_csv('mats/bcsstk08.mtx3.csv')
# f5 = pd.read_csv('mats/bcsstk08.mtx4.csv')
# f6 = pd.read_csv('mats/bcsstk08.mtx5.csv')
#
# #numOfColor,mnat,mnew,mlfo,msat
# f7 = dfs[0]
# print(f7)
# for i in range(1,10):
#     f7.append(dfs[i], ignore_index=True)
# # f7 = f1.append(f2, ignore_index=True).append(f3, ignore_index=True).append(f4, ignore_index=True).append(f5, ignore_index=True).append(f6, ignore_index=True)
# print(f7)
# # threedee = plt.figure().gca(projection='3d')
# # print(color)
# #
print(df)
threedee = plt.figure().gca(projection='3d')
threedee.scatter(df.numOfColor, df.mnat, df.k)
threedee.set_xlabel('Number of Color')
threedee.set_ylabel('Discovered')
threedee.set_zlabel('k')
plt.show()
