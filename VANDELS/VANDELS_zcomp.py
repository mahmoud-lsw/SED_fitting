import numpy as np
import pylab

specz = np.loadtxt("photoz_V_2comp_LP1.txt", usecols=(1,))
photz = np.loadtxt("../VANDELS_data/list_LP_photz.txt", usecols=(4,))

pylab.figure()
pylab.scatter(specz, photz)
pylab.plot([0, 7], [0, 7], color="black")
pylab.show()
