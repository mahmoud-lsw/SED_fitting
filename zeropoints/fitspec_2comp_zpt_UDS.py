import numpy as np
import sys
from astropy.cosmology import WMAP9

firstobj = int(sys.argv[1])
lastobj = int(sys.argv[2])

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit
wavs = np.array([3776.4, 4426.9, 5454.6, 6506.3, 7645.3, 9010.9, 9093.0, 9183.0, 10200.1, 12483.2, 16318.0, 22010.2]) # 3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8

wavs = wavs*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393


### Load up the galaxy catalogue and list of stellar population ages to be fit
all_obj_specz = np.loadtxt("../../VANDELS_data/UDS_zphot_training_phot_zless5.cat", usecols=(25,))
all_obj_fluxes = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_training_phot_zless5.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2)
all_obj_fluxerrs = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_training_phot_zless5.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)
#ages = np.loadtxt("ages.txt")
ages = np.array([508800000.0, 806400000.0, 1139100030.0, 1434000000.0, 1700000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0])


### For all objects which have not been observed in a given band, blow up the associated error so fit is not affected
for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0] == 0.:
            all_obj_fluxerrs[i, j, 0] = 9999999999.

### Build array for photometric redshift and parameter outputs
output = np.zeros((lastobj-firstobj+1)*9, dtype="float")
output.shape = ((lastobj-firstobj+1), 9)
output[:,8] = 9999999999.


### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined
zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = WMAP9.lookback_time(zarr).value
flux_diffs = np.copy(all_obj_fluxes[firstobj:lastobj+1, :, 0])
f_old_array = 1.01 - np.logspace(-2, 0, 51)

for m in range(firstobj, lastobj+1):
    print "object " + str(m)
    
    obj_fluxes = np.expand_dims(all_obj_fluxes[m, :, :], axis=0)
    obj_fluxerrs = np.expand_dims(all_obj_fluxerrs[m, :, :], axis=0)
    for j in range(9): #9 old ages to iterate over
        arg = 0
        z = all_obj_specz[m]
        modelno = int((np.round(z, 2)*100.) - 1.)
        while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 499:
            arg = arg+1
        for l in range(3):
            th_mag_array_old = np.expand_dims(np.loadtxt("../../models/burst/synmagsUDS_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[modelno, :], axis=0)
            th_mag_array_new = np.expand_dims(np.loadtxt("../../models/const/synmagsUDS_age_" + str(newages[l]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[modelno, :], axis=0)
        th_flux_array_old_raw = 10**((23.9-th_mag_array_old)/2.5)
        th_flux_array_new_raw = 10**((23.9-th_mag_array_new)/2.5)
        for i in range(51):
            f_old_V = f_old_array[i]
            if f_old_V == 1.:
                th_flux_array_new = 0.*th_flux_array_new_raw
                th_flux_array_old = np.copy(th_flux_array_old_raw)
            else:
                th_flux_array_new = np.copy(th_flux_array_new_raw)
                th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw
                for k in range(61):
                    EBV = 0.025*k
                    th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new).T*np.expand_dims((10**((-EBV*coef)/2.5)).T, axis=2), axis=0) #microjanskys
                    const =  np.expand_dims(np.sum(obj_fluxes*th_flux_array/obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/obj_fluxerrs**2, axis=1), axis=1)
                    chivals = np.sum((obj_fluxes/obj_fluxerrs - const*(th_flux_array/obj_fluxerrs))**2, axis=1)
                    if np.min(chivals) < output[m-firstobj,8]:
                        zmin = np.argmin(chivals)
                        flux_diffs[m-firstobj, :] = np.squeeze(obj_fluxes/(th_flux_array*const))
                        output[m-firstobj,:] = np.array([m+1, all_obj_specz[m], 0.01*(zmin+1), ages[j], newages[l], f_old_V, EBV, const, chivals])
                                   
        np.savetxt("photoz_2comp_UDS_" + str(firstobj) + "_" + str(lastobj) + ".txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi") 
        np.savetxt("flux_fracoffsets_2comp_UDS_" + str(firstobj) + "_" + str(lastobj) + ".txt", flux_diffs)
        
