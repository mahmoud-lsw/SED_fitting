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
    


hdulist = np.loadtxt("../../VANDELS_data/list_extra_LP.txt", usecols=(0,), dtype="str")
IDlist = np.loadtxt("../../VANDELS_data/list_extra_LP.txt", usecols=(5,), dtype="str")

#hdulist = hdulist[0:100]
#IDlist = IDlist[0:100]

wavzpt = fits.open("../../VANDELS_data/spectra-2015-12-09-LP/sc_UDS001001_P1M2Q4_003_1.fits")[0].header["CRVAL1"]
dwav = fits.open("../../VANDELS_data/spectra-2015-12-09-LP/sc_UDS001001_P1M2Q4_003_1.fits")[0].header["CDELT1"]
fluxes1 = fits.open("../../VANDELS_data/spectra-2015-12-09-LP/sc_UDS001001_P1M2Q4_003_1.fits")[0].data
maxwav = wavzpt + dwav*(len(fluxes1))
objwavs_nobin = np.arange(wavzpt, maxwav, dwav)
spec = np.zeros(len(fluxes1)*3*len(hdulist), dtype="float")
spec.shape = (len(hdulist), len(fluxes1), 3)
binspec = spec[:,0:784,:]

objlist = np.zeros(6*len(hdulist), dtype="float")
objlist.shape = (len(hdulist), 6)
objlist[:,0] = np.arange(1, len(hdulist)+1)


for m in range(len(hdulist)):
    hdu1 = fits.open("../../VANDELS_data/spectra-2015-12-09-LP/" + hdulist[m])
    fluxes1 = hdu1[0].data
    noise1 = hdu1[3].data
    objlist[m,1] = float(IDlist[m][3:])
    
    spec[m, :,0] = objwavs_nobin
    spec[m, :,1] = fluxes1
    spec[m, :,2] = noise1
    binspec[m,:,:] = specbin(spec[m,:,:], 2)[39:-254,:]
    objlist[m,2] = np.mean(binspec[m,0:6,1]/binspec[m,0:6,2])
    objlist[m,3] = np.mean(binspec[m,95:101,1]/binspec[m,95:101,2])
    objlist[m,4] = np.mean(binspec[m,192:198,1]/binspec[m,192:198,2])
    objlist[m,5] = np.mean(binspec[m,290:296,1]/binspec[m,290:296,2])

"""
a = []
for i in range(len(IDlist)):
    if IDlist[i] not in a:
        a.append(IDlist[i])
    elif IDlist[i] == 0:
        a.append(IDlist[i])
    elif IDlist[i] == 999999999:
        a.append(IDlist[i])
        IDlist
        
spec = np.zeros(len(fluxes1)*3*len(a), dtype="float")
spec.shape = (len(a), len(fluxes1), 3)
objlist = np.zeros(7*len(a), dtype="float")
objlist.shape = (len(a), 7)
objlist[:,0] = np.arange(1, len(a)+1)
objlist[:,1] = np.array(a)
binspec = spec[:,0:784,:]

z=0
nn=0
for m in range(len(a)):
    i=0
    ind = []
    for n in range(len(hdulist)):
        if a[m] == IDlist[n]: #and a[m] != 0 and a[m] != 999999999
            i+=1
            ind.append(n)
    if a[m] == 0:
        i=1
        ind = [ind[z]]
        z+=1

    if a[m] == 999999999:
        i=1
        ind = [ind[nn]]
        nn+=1
        
    print ind 
    if i == 2:
        objlist[m,2] = 1.
        hdu1 = fits.open("../../VANDELS_data/spectra/" + hdulist[ind[0]]) #sc_206806_UDS_P1M1_MR_Q1_029_1.fits
        hdu2 = fits.open("../../VANDELS_data/spectra/" + hdulist[ind[1]])
        fluxes1 = hdu1[4].data
        noise1 = hdu1[3].data
        fluxes2 = hdu2[4].data
        noise2 = hdu2[3].data
        
        mednoise1 = np.median(noise1)
        mednoise2 = np.median(noise2)

        for k in range(len(noise1)):
            if np.abs(noise1[k] - mednoise1) > mednoise1*3:
                noise1[k] = mednoise1
            if np.abs(noise2[k] - mednoise2) > mednoise2*3:
                noise2[k] = mednoise2
                
        spec[m, :,0] = objwavs_nobin
        spec[m, :,1] = (fluxes1 + fluxes2)/2.
        spec[m, :,2] = np.sqrt(noise1**2 + noise2**2)/2.

        binspec[m,:,:] = specbin(spec[m,:,:], 2)[98:-254,:]
        for k in range(784):
            if binspec[m,k,2] == 0.:
                binspec[m,k,2] = 99999999.

    elif i == 1:
        hdu1 = fits.open("../../VANDELS_data/spectra/" + hdulist[ind[0]])
        fluxes1 = hdu1[4].data
        noise1 = hdu1[3].data
        mednoise = np.median(noise1)

        for k in range(len(noise1)):
            if np.abs(noise1[k] - mednoise) > mednoise*3:   
                noise1[k] = mednoise
            if noise1[k] == 0.:
                noise1[k] = 99999999.
                
        spec[m, :,0] = objwavs_nobin
        spec[m, :,1] = fluxes1
        spec[m, :,2] = noise1
        binspec[m,:,:] = specbin(spec[m,:,:], 2)[98:-254,:]
        for k in range(784):
            if binspec[m,k,2] == 0.:
                binspec[m,k,2] = 99999999.
                
    else:
        print "shit..."
        raw_input()
"""      


np.savetxt("../../VANDELS_data/combined_spec_LP.txt", binspec[:,:,1])
np.savetxt("../../VANDELS_data/combined_spec_errs_LP.txt", binspec[:,:,2])
np.savetxt("../../VANDELS_data/wavs_LP.txt", binspec[0,:,0])
np.savetxt("../../VANDELS_data/ACC_obj_list_LP.txt", objlist, header="ACC_ID VANDELS_ID SNR_5000A SNR_6000A SNR_7000A SNR_8000A")
