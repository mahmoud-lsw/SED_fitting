import numpy as np
import sys
from astropy.cosmology import WMAP9

foldmin = int(sys.argv[1])
foldmax = int(sys.argv[2])

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit

zarr = np.arange(0.01, 7.01, 0.01)
wavs = np.zeros(12*699)
wavs.shape = (12, 699)

for i in range(699):
    wavs[:,i] = np.array([3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8])/(1+zarr[i])

#wavs = np.array([3776.4, 4426.9, 5454.6, 6506.3, 7645.3, 9010.9, 9093.0, 9183.0, 10200.1, 12483.2, 16318.0, 22010.2]) # 3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8

wavs = wavs*10**-4
coef = np.copy(wavs)
for i in range(699):
    for j in range(12):
        if wavs[j,i] < 0.63:
            coef[j, i] = ((0.013/(wavs[j,i]**3)) - (0.232/(wavs[j,i]**2)) + (1.766/wavs[j,i]) - 0.743)
        else:
            coef[j, i] = 1.217/wavs[j,i] - 0.393

### Load up the galaxy catalogue and list of stellar population ages to be fit
all_obj_specz = np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(27,))#[0:10]   ###   galaxy_photom_specz_adam_zless5.cat  27
all_obj_fluxes = np.expand_dims(np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2)#[0:10,:,:] #../VANDELS_data/UDS_zphot_training_phot.cat
all_obj_fluxerrs = np.expand_dims(np.loadtxt("galaxy_photom_specz_adam_zless5.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)#[0:10,:,:]

#####
#all_obj_specz = np.zeros(len(all_obj_fluxes))
#####

#ages = np.loadtxt("ages.txt")
ages = np.array([508800000.0, 806400000.0, 1139100030.0, 1434000000.0, 1700000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0])


### Apply offsets to the data as calculated by code in zeropoints/ in order to deal with various systematic calibration errors
offsets = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/mean_ratios_2comp_UDS.txt", usecols=(0,)), axis=0), axis=2)
offset_errs = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/mean_ratios_2comp_UDS.txt", usecols=(1,)), axis=0), axis=2)


for i in range(12): # If the offset is consistent with zero to within 1 sigma do nothing
    if np.abs(offsets[0, i, 0] - 1.) > 0.:# offset_errs[0, i, 0]:
        offsets[0, i, 0] = 1.
        offset_errs[0, i, 0] = 0.


#all_obj_fluxes = all_obj_fluxes/offsets
#all_obj_fluxerrs = all_obj_fluxes*np.sqrt((all_obj_fluxerrs/all_obj_fluxes)**2 + (offset_errs/offsets**3)**2)


### For all objects which have not been observed in a given band, blow up the associated error so fit is not affected
for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0] == 0.:
            all_obj_fluxerrs[i, j, 0] = 9999999999.

### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*9, dtype="float")
output.shape = (len(all_obj_fluxes), 9)
output[:,8] = 9999999999.


### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined
lbtarr = WMAP9.lookback_time(zarr).value
f_old_array = 1.01 - np.logspace(-2, 0, 51)

for j in range(9): #8 old ages to iterate over
    arg = 0
    while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 699:
        arg = arg+1
    for l in range(3):
        th_mag_array_old = np.loadtxt("../models/burst/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[:arg, :] ####UDS
        th_mag_array_new = np.loadtxt("../models/const/synmags_age_" + str(newages[l]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[:arg, :]
        th_flux_array_old_raw = 10**((23.9-th_mag_array_old)/2.5)
        th_flux_array_new_raw = 10**((23.9-th_mag_array_new)/2.5)
        
        for i in range(foldmin, foldmax+1):
            f_old_V = f_old_array[i]
            print "Initial burst age: " + str(ages[j]) + ", fitting to redshift " + str((arg+1)*0.01) + ", with new starburst of age " + str(newages[l]) + " and f_old_V " + str(f_old_V)
            if f_old_V == 1.:
                th_flux_array_new = 0.*th_flux_array_new_raw
                th_flux_array_old = np.copy(th_flux_array_old_raw)
            else:
                th_flux_array_new = np.copy(th_flux_array_new_raw)
                th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw
            for k in range(61):
                EBV = 0.025*k
                #th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new).T*np.expand_dims((10**((-EBV*coef)/2.5)).T, axis=2), axis=0) #microjanskys
                th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new).T*(10**((-EBV*coef[:, :arg])/2.5)), axis=0) #microjanskys
                const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
                chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
                for m in range(len(chivals)):
                    if np.min(chivals[m, :]) < output[m,8]:
                        zmin = np.argmin(chivals[m, :])
                        output[m,:] = np.array([m+1, all_obj_specz[m], 0.01*(zmin+1), ages[j], newages[l], f_old_V, EBV, const[m, 0, zmin], chivals[m, zmin]])
    np.savetxt("photoz_" + str(foldmin) + "_" + str(foldmax) + ".txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi") #_UDS
    

