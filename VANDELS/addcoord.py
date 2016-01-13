import numpy as np
from astropy.io import fits

basetable = np.loadtxt("../../VANDELS_data/list_LP.txt", dtype="str")

finaltable = open("../../VANDELS_data/list_extra_LP.txt", "w")
finaltable.write("# fname wavzpt dwav ra dec ID \n")

for i in range(len(basetable)):
    hdulist = fits.open("../../VANDELS_data/spectra-2015-12-09-LP/"+basetable[i])

    wavzpt = hdulist[0].header["CRVAL1"]
    dwav = hdulist[0].header["CDELT1"]
    ra = hdulist[0].header["HIERARCH PND OBJRA"]
    dec = hdulist[0].header["HIERARCH PND OBJDEC"]
    ID = hdulist[0].header["HIERARCH PND OBJID"]
    finaltable.write(basetable[i] + " " + str(wavzpt) + " " +  " " + str(dwav) + " " + str(ra) + " " + str(dec) + " " + str(ID) + "\n")

finaltable.close()
