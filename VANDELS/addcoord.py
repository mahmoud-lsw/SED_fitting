import numpy as np
from astropy.io import fits

basetable = np.loadtxt("../../VANDELS_data/list.txt", dtype="str")

finaltable = open("../../VANDELS_data/list_extra.txt", "w")
finaltable.write("# fname wavzpt dwav ra dec ID \n")

for i in range(len(basetable)):
    hdulist = fits.open("../../VANDELS_data/spectra/"+basetable[i,0])

    wavzpt = hdulist[0].header["CRVAL1"]
    dwav = hdulist[0].header["CDELT1"]
    ra = hdulist[0].header["HIERARCH PND OBJRA"]
    dec = hdulist[0].header["HIERARCH PND OBJDEC"]
    ID = hdulist[0].header["HIERARCH PND OBJID"]
    finaltable.write(basetable[i,0] + " " + str(wavzpt) + " " +  " " + str(dwav) + " " + str(ra) + " " + str(dec) + " " + str(ID) + "\n")

finaltable.close()


