import numpy as np
import sys
from astropy.cosmology import WMAP9

foldmin = int(sys.argv[1])
foldmax = int(sys.argv[2])

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit
wavs = np.array([3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8])
wavs = wavs*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393


### Load up the galaxy catalogue and list of stellar population ages to be fit
all_obj_specz = np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(27,))
all_obj_fluxes = np.expand_dims(np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2)
all_obj_fluxerrs = np.expand_dims(np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)
#ages = np.loadtxt("ages.txt")
ages = np.array([508800000.0, 806400000.0, 1139100030.0, 1434000000.0, 1700000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0])

### Apply offsets to the data as calculated by code in zeropoints/ in order to deal with various systematic calibration errors
offsets = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/median_ratios_2comp.txt", usecols=(0,)), axis=0), axis=2)
offset_errs = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/median_ratios_2comp.txt", usecols=(1,)), axis=0), axis=2)

for i in range(12): # If the offset is consistent with zero to within 1 sigma do nothing
    if np.abs(offsets[0, i, 0] - 1.) > 0.:# offset_errs[0, i, 0]:
        offsets[0, i, 0] = 1.
        offset_errs[0, i, 0] = 0.

all_obj_fluxes = all_obj_fluxes/offsets
all_obj_fluxerrs = all_obj_fluxes*np.sqrt((all_obj_fluxerrs/all_obj_fluxes)**2 + (offset_errs/offsets**3)**2)


### For all objects which have not been observed in a given band, blow up the associated error so fit is not affected
for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0] == 0.:
            all_obj_fluxerrs[i, j, 0] = 9999999999.

### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*10, dtype="float")
output.shape = (len(all_obj_fluxes), 10)
output[:,9] = 9999999999.


### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined
zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = WMAP9.lookback_time(zarr).value
    
for j in range(9): #8 old ages to iterate over
    arg = 0
    while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 499:
        arg = arg+1
    for l in range(3):
        th_mag_array_old = np.loadtxt("models/burst/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[:arg, :]
        th_mag_array_new = np.loadtxt("models/const/synmags_age_" + str(newages[l]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[:arg, :]
        th_flux_array_old = np.copy(th_mag_array_old)
        th_flux_array_new = np.copy(th_mag_array_new)
        for i in range(foldmin, foldmax+1):
            f_old_V = 0.02*i
            print "Initial burst age: " + str(ages[j]) + ", fitting to redshift " + str((arg+1)*0.01) + ", with new starburst of age " + str(newages[l]) + " and f_old_V " + str(f_old_V)
            if f_old_V == 1.:
                old_modifier = np.expand_dims(np.ones(arg, dtype="float"), axis=1)
                th_flux_array_new = 0.*th_mag_array_new
            else:    
                old_modifier = (f_old_V/(1.-f_old_V))*10**((-(np.expand_dims(th_mag_array_new[:,3] - th_mag_array_old[:,3], axis=1)))/2.5)
                th_flux_array_new = 10**((23.9-th_mag_array_new)/2.5)
            th_flux_array_old = old_modifier*10**((23.9-th_mag_array_old)/2.5)
            for k in range(61):
                EBV = 0.025*k
                th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new).T*np.expand_dims((10**((-EBV*coef)/2.5)).T, axis=2), axis=0) #microjanskys
                const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
                chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
                for m in range(len(chivals)):
                    if np.min(chivals[m, :]) < output[m,9]:
                        zmin = np.argmin(chivals[m, :])
                        output[m,:] = np.array([m+1, all_obj_specz[m], 0.01*(zmin+1), ages[j], newages[l], f_old_V, old_modifier[zmin], EBV, const[m, 0, zmin], chivals[m, zmin]])
    np.savetxt("photoz_" + str(foldmin) + "_" + str(foldmax) + ".txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V old_modifier EBV norm chi")
    

