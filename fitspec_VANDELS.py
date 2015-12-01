import numpy as np
import sys
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
import pylab

foldmin = int(sys.argv[1])
foldmax = int(sys.argv[2])

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
    
### Load up the galaxy catalogue and list of stellar population ages to be fit
hdulist1 = fits.open("../VANDELS_data/spectra/sc_206806_UDS_P1M1_MR_Q1_029_1.fits")
wavzpt = hdulist1[0].header["CRVAL1"]
dwav = hdulist1[0].header["CDELT1"]
fluxes1 = hdulist1[4].data#*10**19
noise1 = hdulist1[3].data#*10**19
maxwav = wavzpt + dwav*(len(fluxes1))
objwavs_nobin = np.arange(wavzpt, maxwav, dwav)
mednoise = np.median(noise1)

for i in range(len(noise1)):
    if np.abs(noise1[i] - mednoise) > mednoise*3:
        noise1[i] = mednoise

spec = np.zeros(len(fluxes1)*3, dtype="float")
spec.shape = (len(fluxes1), 3)
spec[:,0] = objwavs_nobin
spec[:,1] = fluxes1
spec[:,2] = noise1

binspec = specbin(spec, 2)[98:-254,:]

all_obj_fluxes = np.expand_dims(np.expand_dims(binspec[:,1], axis=0), axis=2)
all_obj_fluxerrs =  np.expand_dims(np.expand_dims(binspec[:,2], axis=0), axis=2) #erg/s/A/cm^2

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit

wavs = binspec[:,0]*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393


#ages = np.loadtxt("ages.txt")
oldages = np.array([508800000.0, 806400000.0, 1139100030.0, 1434000000.0, 1700000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0])

"""
### Apply offsets to the data as calculated by code in zeropoints/ in order to deal with various systematic calibration errors
offsets = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/ratios/mean_ratios_with_errors.txt", usecols=(0,)), axis=0), axis=2)
offset_errs = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/ratios/mean_ratios_with_errors.txt", usecols=(1,)), axis=0), axis=2)

for i in range(12): # If the offset is consistent with zero to within 1 sigma do nothing
    if np.abs(offsets[0, i, 0] - 1.) < offset_errs[0, i, 0]:
        offsets[0, i, 0] = 1.
        offset_errs[0, i, 0] = 0.

all_obj_fluxes = all_obj_fluxes/offsets
all_obj_fluxerrs = all_obj_fluxes*np.sqrt((all_obj_fluxerrs/all_obj_fluxes)**2 + (offset_errs/offsets**3)**2)


### For all objects which have not been observed in a given band, blow up the associated error so fit is not affected
for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0] == 0.:
            all_obj_fluxerrs[i, j, 0] = 9999999999.
"""
### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*8, dtype="float")
output.shape = (len(all_obj_fluxes), 8)
output[:,7] = 9999999999.

zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = cosmo.lookback_time(zarr).value
f_old_array = 1.01 - np.logspace(-2, 0, 51)

### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined



for j in range(9): #9 old ages to iterate over
    arg = 0
    while oldages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 499:
        arg = arg+1
    for l in range(3):
        th_flux_array_new_raw = np.loadtxt("models/spec/newconst/age_" + str(newages[l]) + ".txt")[:,:arg+1]
        th_flux_array_old_raw = np.loadtxt("models/spec/oldburst/age_" + str(oldages[j]) + ".txt")[:,:arg+1] #erg/s/A/cm^2
        for i in range(foldmin, foldmax+1):
            f_old_V = f_old_array[i]
            print "Initial burst age: " + str(oldages[j]) + ", fitting to redshift " + str((arg+1)*0.01) + ", with new starburst of age " + str(newages[l]) + " and f_old_V " + str(f_old_V)
            if f_old_V == 1.:
                th_flux_array_new = 0.*th_flux_array_new_raw
                th_flux_array_old = th_flux_array_old_raw
            else:
                th_flux_array_new = np.copy(th_flux_array_new_raw)
                th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw
            for k in range(61):
                EBV = 0.025*k
                th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new)*np.expand_dims((10**((-EBV*coef)/2.5)).T, axis=2), axis=0)
                const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
                chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
                for m in range(len(chivals)):
                    if np.min(chivals[m, :]) < output[m,7]:
                        zmin = np.argmin(chivals[m, :])
                        output[m,:] = np.array([m+1, 0.01*(zmin+1), oldages[j], newages[l], f_old_V, EBV, const[m, 0, zmin], chivals[m, zmin]])
    np.savetxt("photoz_V_" + str(foldmin) + "_" + str(foldmax) + ".txt", output, header="obj_no phot_z age_old age_new f_old_V  EBV norm chi")
    

"""
param = np.loadtxt("photoz_V_0_50.txt") 
print "obj_no phot_z age_old age_new f_old_V EBV norm chi"
print param
oldage = np.argmin(np.abs(oldages - param[2]))
newage = np.argmin(np.abs(newages - param[3]))
redshift_ind = (param[1]/0.01)-1
f_old_V = param[4]
old_modifier = param[5]
EBV = param[6]
const = param[7]

print "z: " + str(param[1])
print "old age: " + str(oldages[oldage]) + " new age: " + str(newages[newage])
print "E(B-V): " + str(EBV)
print "f_old: " + str(f_old_V)
print "old_modifier: " + str(old_modifier)
print "const: " + str(const)
print "Stellar Mass: " + str(np.round(const*(old_modifier + newages[newage])/10**9, 3)) + "*10^9 Solar masses"

th_flux_array_new_raw = np.loadtxt("models/spec/newconst/age_" + str(newages[newage]) + ".txt")[:,redshift_ind:redshift_ind+1]#[:,:arg+1]
th_flux_array_old_raw = np.loadtxt("models/spec/oldburst/age_" + str(oldages[oldage]) + ".txt")[:,redshift_ind:redshift_ind+1]#[:,:arg+1] #erg/s/A/cm^2
if f_old_V == 1.:
    old_modifier = np.expand_dims(np.ones(arg+1, dtype="float"), axis=0)
    th_flux_array_new = 0.*th_flux_array_new_raw
else:    
    old_modifier = (f_old_V/(1.-f_old_V))*np.expand_dims(th_flux_array_new_raw[0,:]/th_flux_array_old_raw[0,:], axis=1)
    th_flux_array_new = np.copy(th_flux_array_new_raw)
th_flux_array_old = old_modifier*th_flux_array_old_raw
th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new)*np.expand_dims((10**((-EBV*coef)/2.5)).T, axis=2), axis=0)


pylab.figure()
pylab.plot(binspec[:,0], np.squeeze(all_obj_fluxes), color="blue")
pylab.plot(binspec[:,0], np.squeeze(th_flux_array*const), color="red", zorder=10)
pylab.ylim(1.1*np.min(all_obj_fluxes), 1.1*np.max(all_obj_fluxes))
pylab.xlabel("Wavelength (A)", size=16)
pylab.ylabel("F_lambda (erg/s/A/cm^2)", size=16)
pylab.show()
"""
