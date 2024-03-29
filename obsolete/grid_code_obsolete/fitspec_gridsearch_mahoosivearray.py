import numpy as np
import sys
from astropy.cosmology import WMAP9

minage = 0#int(sys.argv[1])
maxage = 16#int(sys.argv[2])

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit
wavs = np.array([3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8])
wavs = wavs*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393

coef = np.expand_dims(np.expand_dims(np.expand_dims(coef, axis=0), axis=2), axis=3)


### Load up the galaxy catalogue and list of stellar population ages to be fit
all_obj_specz = np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(27,))
all_obj_fluxes = np.expand_dims(np.expand_dims(np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2), axis=3)
all_obj_fluxerrs = np.expand_dims(np.expand_dims(np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2), axis=3)

ages = np.loadtxt("ages.txt")


### For all objects which have not been observed in a given band, blow up the associated error so fit is not affected
for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0, 0] == 0.:
            all_obj_fluxerrs[i, j, 0, 0] = 9999999999.


### Apply offsets to the data as calculated by code in offsets/ in order to deal with various systematic calibration errors
offsets = np.expand_dims(np.expand_dims(np.expand_dims(np.loadtxt("offsets/mean_offsets_with_errors.txt", usecols=(0,)), axis=0), axis=2), axis=3)
offset_errs = np.expand_dims(np.expand_dims(np.expand_dims(np.loadtxt("offsets/mean_offsets_with_errors.txt", usecols=(1,)), axis=0), axis=2), axis=3)

for i in range(12): # If the offset is consistent with zero to within 1 sigma do nothing
    if np.abs(offsets[0, i, 0]) < offset_errs[0, i, 0]:
        offsets[0, i, 0] = 0.
        offset_errs[0, i, 0] = 0.

all_obj_fluxes = all_obj_fluxes + offsets
all_obj_fluxerrs = np.sqrt(all_obj_fluxerrs**2 + offset_errs**2)


### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*8, dtype="float")
output.shape = (len(all_obj_fluxes), 8)
output[:,7] = 9999999999.


### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined
tauvals = np.array([0.05, 1, 10], dtype="str")
zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = WMAP9.lookback_time(zarr).value
EBVarr = np.expand_dims(np.expand_dims(np.expand_dims(np.arange(0, 1.50001, 0.025), axis=0), axis=0), axis=0)
    
for j in range(minage, maxage+1): 
    arg = 0
    while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 499:
        arg = arg+1
    for l in range(3):
        print "Age: " + str(ages[j]) + ", fitting to redshift " + str((arg+1)*0.01) + ", with tau: " + tauvals[l]
        th_mag_array = np.expand_dims(np.expand_dims(np.loadtxt("models/" + tauvals[l] + "/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[0:arg, :].T, axis=0), axis=3)
        print "th_mag_array computed"
        print th_mag_array.nbytes
        th_flux_array = (10**((23.9 - EBVarr*coef - th_mag_array)/2.5)) #microjanskys
        print "th_flux_array computed"
        print th_flux_array.nbytes
        const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
        print "const array computed"
        print const.nbytes
        chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
        print chivals.nbytes
        print "chivals computed"
        for m in range(len(all_obj_specz)):
            if np.min(chivals[m, :, :]) < output[m,7]:
                zmin, EBVmin = np.unravel_index(chivals[m,:,:].argmin(), chivals[m,:,:].shape)
                output[m,:] = np.array([m+1, all_obj_specz[m], 0.01*(zmin+1), ages[j], float(tauvals[l]), 0.025*(EBVmin+1), const[m, 0, zmin, EBVmin], chivals[m, zmin, EBVmin]])
    np.savetxt("photoz_" + str(minage) + "_" + str(maxage) + ".txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")
