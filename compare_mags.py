import numpy as np

data = np.loadtxt("mag_comparison.cat", usecols=(2,3,4,5,6,7,8,9,10,11,12,13,14))
mydata = np.loadtxt("1Gyrmags_70_0.3_cosmo.txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13))

print data[5,:] - mydata[2,:]
print data[5,:]
