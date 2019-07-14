import pandas as pd, matplotlib.pyplot as plt, numpy as np
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams.update({'font.size': 12})

df = pd.read_csv('bcsstk08_res.csv')
threedee = plt.figure().gca(projection='3d')
threedee.scatter(df.k, df.numOfColor, df.nat, color='red')
threedee.set_xlabel('k')
threedee.set_ylabel('Num of colors')
threedee.set_zlabel('Delivered')

threedee.scatter(df.k, df.numOfColor, df.lfo, color='blue')
threedee.set_xlabel('k')
threedee.set_ylabel('Num of colors')
threedee.set_zlabel('Delivered')

threedee.scatter(df.k, df.numOfColor, df.sat, color='green')
threedee.set_xlabel('k')
threedee.set_ylabel('Num of colors')
threedee.set_zlabel('Delivered')

threedee.scatter(df.k, df.numOfColor, df.ago, color='yellow')
threedee.set_xlabel('k')
threedee.set_ylabel('Num of colors')
threedee.set_zlabel('Delivered')

plt.show()

import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3
fig=p.figure()
ax = p3.Axes3D(fig)
ax.plot_wireframe(df.k, df.numOfColor, df.mnat)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
p.show()
plt.show()