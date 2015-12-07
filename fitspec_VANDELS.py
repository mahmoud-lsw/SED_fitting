import numpy as np
import sys
from astropy.cosmology import WMAP9 as cosmo
import time

foldmin = int(sys.argv[1])
foldmax = int(sys.argv[2])

    
### Load up the galaxy catalogue and list of stellar population ages to be fit

all_obj_fluxes = np.expand_dims(np.loadtxt("../VANDELS_data/combined_spec.txt"), axis=2)
all_obj_fluxerrs =  np.expand_dims(np.loadtxt("../VANDELS_data/combined_spec_errs.txt"), axis=2) #erg/s/A/cm^2


### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit

wavs = np.loadtxt("../VANDELS_data/wavs.txt")*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393


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
            t0 = time.time()
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
            t1 = time.time()
            print "time taken: " + str(t1-t0)
    np.savetxt("photoz_V_" + str(foldmin) + "_" + str(foldmax) + ".txt", output, header="obj_no phot_z age_old age_new f_old_V  EBV norm chi")
    
