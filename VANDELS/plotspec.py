import numpy as np
from astropy.io import fits
import pylab

def specbin(spec, binn):
    if int(2*len(spec)/binn) != 2*int(len(spec)/binn):
        binned = np.zeros(2*len(spec)/binn -1)
    else:
        binned = np.zeros(2*len(spec)/binn)
    
    binned.shape = (len(spec)/binn, 2)
    for i in range(len(binned)):
        binned[i, 0] = np.mean(spec[binn*i:binn*i+binn, 0])
        binned[i, 1] = np.mean(spec[binn*i:binn*i+binn, 1])
    return binned

hdulist1 = fits.open("../../VANDELS_data/spectra/sc_206806_UDS_P1M1_MR_Q1_029_1.fits")
print hdulist1.info()
#sc_206024_UDS_P1M1_MR_Q1_018_3.fits

wavzpt = hdulist1[0].header["CRVAL1"]
dwav = hdulist1[0].header["CDELT1"]
print dwav




sky1 = hdulist1[2].data

fluxes1 = hdulist1[4].data#*10**19
noise1 = hdulist1[3].data#*10**19
mednoise = np.median(noise1)

for i in range(len(noise1)):
    if np.abs(noise1[i] - mednoise) > mednoise*3:
        noise1[i] = mednoise

#hdulist2 = fits.open("spectra/sc_191361_UDS_P1M1_MR_Q1_053_1.fits")
#print hdulist2.info()

#fluxes2 = hdulist2[0].data#*10**19

#sky2 = hdulist2[2].data
#noise2 = hdulist2[3].data#*10**19

maxwav = wavzpt + dwav*(len(fluxes1))
wavs = np.arange(wavzpt, maxwav, dwav)

spec = np.zeros(2*len(fluxes1), dtype="float")
spec.shape = (len(fluxes1), 2)

spec[:,1] = fluxes1
spec[:,0] = wavs

binspec = specbin(spec, 2)

pylab.figure()
pylab.plot(wavs, fluxes1, color="black", zorder=5)
#pylab.plot(wavs, sky1, color="blue")
pylab.plot(wavs, noise1, color="blue", zorder=1)
#pylab.plot(binspec[:,0], binspec[:,1], color="red", zorder=8)
pylab.xlabel("Wavelength (Angstroms)", size=16)
pylab.plot([5000, 9000], [0., 0.], color="black", zorder=0)
pylab.xlim(5000, 9000)
pylab.ylim(-0.5*10**-18, 2*10**-18)
pylab.ylabel("Flux", size=16)
pylab.show()

""" 
pylab.plot(wavs, fluxes2, color="green")
#pylab.plot(wavs, sky, color="blue")
#pylab.plot(wavs, noise2, color="green")
#print noise2[1000:1100]
pylab.show()
"""
#lines = np.array([6565, 6550, 5008, 4960, 4863, 4342, 3727, 2799, 1216])
"""
pylab.figure()
pylab.plot(wavs, (fluxes1+fluxes2)/2., color="black")
pylab.plot(wavs, noise1/np.sqrt(2), color="red")
pylab.xlabel("Wavelength (Angstroms)", size=16)
pylab.xlim(5000, 9000)
pylab.ylim(-0.5*10**-18, 2*10**-18)
pylab.ylabel("Flux", size=16)
pylab.plot([5000, 9000], [0., 0.], color="black")
for i in range(9):
    pylab.plot(np.array([lines[i], lines[i]])*2.07, np.array([0, 1*10**-18]), color="blue")
#pylab.plot(wavs, sky, color="blue")
#pylab.plot(wavs, noise2, color="green")
#print noise2[1000:1100]
pylab.show()
"""
