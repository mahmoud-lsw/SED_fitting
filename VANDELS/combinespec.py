import numpy as np
from astropy.io import fits

"""bins up spectrum consisting of a column of wavelength values, a column of fluxes and a column of flux errors by factor binn"""
def specbin(spec, binn): 
    if int(2*len(spec)/binn) != 2*int(len(spec)/binn):
        binned = np.zeros(3*len(spec)/binn -1)
    else:
        binned = np.zeros(3*len(spec)/binn)
    
    binned.shape = (len(spec)/binn, 3)
    for i in range(len(binned)):
        binned[i, 0] = np.mean(spec[binn*i:binn*i+binn, 0])
        binned[i, 1] = np.mean(spec[binn*i:binn*i+binn, 1])
        binned[i, 2] = np.sqrt(np.sum(spec[binn*i:binn*i+binn, 2]**2))
    return binned
    


hdulist = np.loadtxt("../../VANDELS_data/list_extra.txt", usecols=(0,), dtype="str")
fluxes1 = fits.open("../../VANDELS_data/spectra/sc_206806_UDS_P1M1_MR_Q1_029_1.fits")[0].data
spec = np.zeros(len(fluxes1)*3*len(hdulist), dtype="float")
spec.shape = (len(hdulist), len(fluxes1), 3)
binspec = spec[:,0:784,:]

for m in range(len(hdulist)):
    hdu1 = fits.open("../../VANDELS_data/spectra/sc_206806_UDS_P1M1_MR_Q1_029_1.fits")

    wavzpt = hdu1[0].header["CRVAL1"]
    dwav = hdu1[0].header["CDELT1"]
    fluxes1 = hdu1[4].data#*10**19
    noise1 = hdu1[3].data#*10**19
    maxwav = wavzpt + dwav*(len(fluxes1))
    objwavs_nobin = np.arange(wavzpt, maxwav, dwav)
    mednoise = np.median(noise1)

    for i in range(len(noise1)):
        if np.abs(noise1[i] - mednoise) > mednoise*3:
            noise1[i] = mednoise
            
    spec[m, :,0] = objwavs_nobin
    spec[m, :,1] = fluxes1
    spec[m, :,2] = noise1

    binspec[m,:,:] = specbin(spec[m,:,:], 2)[98:-254,:]

np.savetxt("../../VANDELS_data/combined_spec.txt", binspec[:,:,1])
np.savetxt("../../VANDELS_data/combined_spec_errs.txt", binspec[:,:,2])
np.savetxt("../../VANDELS_data/wavs.txt", binspec[0,:,0])
