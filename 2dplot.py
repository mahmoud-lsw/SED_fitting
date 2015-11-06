import numpy as np
import pylab
import matplotlib as mpl

chivals = np.loadtxt("synmags_T0.05/chiarrays/obj_3.txt")

x, y = np.unravel_index(chivals.argmin(), chivals.shape)

flatchivals = np.ndarray.flatten(chivals)

no = 0.

for i in range(len(flatchivals)):
    if flatchivals[i] > 9999999998.:
        no = no+1

print no, no/len(flatchivals)

pylab.figure()
cmap = mpl.colors.ListedColormap(['white','cyan','blue', "black", "red"])
pylab.plot(y, x, "x")
bounds=[10,100,1000,10000,100000, 9999999998.]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
pylab.imshow(chivals,interpolation='nearest', cmap = cmap, norm=norm, origin=np.array([80., 800.])) #"YlOrRd"
pylab.xlabel("Age'", size=16)
pylab.ylabel("Redshift", size=16)
pylab.show()

