import scipy.io as sio, numpy as np


mat = sio.mmread('nos3.mtx')
mat = np.array(mat.toarray())
mat[mat!=0] = 1
print(mat.sum())
