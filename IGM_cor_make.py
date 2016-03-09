import numpy as np
import pylab
from scipy.integrate import quad

#Ross = np.loadtxt("../Downloads/DA_DB.cat")

def D_a_eval(wav):
    return np.exp(-(3.6*10**-3)*(wav/1216.)**3.46)
    
def D_b_eval(wav):
    return np.exp(-(1.7*10**-3)*(wav/1026.)**3.46 - (1.2*10**-3)*(wav/973.)**3.46 - (9.3*10**-4)*(wav/950.)**3.46)

zarr = np.arange(0.01, 7.01, 0.01)

D = np.zeros(len(zarr)*2)
D.shape = (len(zarr), 2)

for i in range(len(zarr)):
    D[i, 0] = 1. - (1/120./(1+zarr[i]))*quad(D_a_eval, 1050*(1+zarr[i]), 1170*(1+zarr[i]))[0]
    D[i, 1] = 1. - (1/95./(1+zarr[i]))*quad(D_b_eval, 920*(1+zarr[i]), 1015*(1+zarr[i]))[0]
    
""" 
pylab.figure()
pylab.plot(zarr, D[:,0], color="blue")
pylab.plot(zarr, Ross[1:701,1], color="blue")
pylab.plot(zarr, D[:,1], color="red")
pylab.plot(zarr, Ross[1:701,2], color="red")
pylab.show()
"""

np.savetxt("IGM_Da_Db.txt", D, header="D_a D_b")
