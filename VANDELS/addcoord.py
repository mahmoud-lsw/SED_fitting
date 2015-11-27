import numpy as np
from astropy.io import fits

basetable = np.loadtxt("list.txt", dtype="str")

for i in range(len(basetable)):
    print "spectra/"+basetable[i,0]
    hdulist = fits.open("spectra/"+basetable[i,0])

