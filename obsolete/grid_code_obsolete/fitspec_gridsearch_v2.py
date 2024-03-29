import numpy as np
from astropy.cosmology import WMAP9

ages = np.loadtxt("ages.txt")

wavs = np.array([3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8])
wavs = wavs*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393

all_obj_specz = np.loadtxt("galaxy_photom_specz_adam.cat", usecols=(27,))
all_obj_fluxes = np.loadtxt("galaxy_photom_specz_adam.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
all_obj_fluxerrs = np.loadtxt("galaxy_photom_specz_adam.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24))

ind = []

for i in range(len(all_obj_fluxes)):
    obj_fluxes = all_obj_fluxes[i,:]
    obj_fluxerrs = all_obj_fluxerrs[i,:]
    if np.min(obj_fluxes) > 0.:
        ind.append(i)

fulldet_obj_fluxes = np.zeros(12*len(ind))
fulldet_obj_fluxes.shape = (len(ind), 12, 1)
fulldet_obj_fluxerrs = np.copy(fulldet_obj_fluxes)
fulldet_obj_specz = np.copy(fulldet_obj_fluxes[:,0,0])

for j in range(len(fulldet_obj_fluxes)):
    fulldet_obj_fluxes[j, :, 0] = all_obj_fluxes[ind[j],:]
    fulldet_obj_fluxerrs[j, :, 0] = all_obj_fluxerrs[ind[j],:]
    fulldet_obj_specz[j] = all_obj_specz[ind[j]]

output = np.zeros(len(ind)*8, dtype="float")
output.shape = (len(ind), 8)
output[:,7] = 9999999999.

tauvals = np.array([0.05, 1, 10], dtype="str")

zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = WMAP9.lookback_time(zarr).value

#3D chi squared array over object, redshift and EBV with max age of stellar pop physically determined.

for j in range(127): 
    print "Age: " + str(ages[j]) + ", fitting to redshift " + str((arg+1)*0.01) # + ", tau: " + tauvals[l] + ", extinction: " + str(EBV)
    arg = 0
    while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 499:
        arg = arg+1
    for l in range(3):
        th_mag_array = np.loadtxt("synmags_T" + tauvals[l] + "/synmags/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[0:arg, :]
        for k in range(61):
            EBV = 0.025*k
            th_flux_array = np.expand_dims((10**((23.9 - EBV*coef - th_mag_array)/2.5)).T, axis=0) #microjanskys
            const =  np.expand_dims(np.sum(fulldet_obj_fluxes*th_flux_array/fulldet_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/fulldet_obj_fluxerrs**2, axis=1), axis=1)
            chivals = np.sum((fulldet_obj_fluxes/fulldet_obj_fluxerrs - const*(th_flux_array/fulldet_obj_fluxerrs))**2, axis=1)
            for m in range(len(chivals)):
                if np.min(chivals[m, :]) < output[m,7]:
                    zmin = np.argmin(chivals[m, :])
                    output[m,:] = np.array([m+1, fulldet_obj_specz[m], 0.01*(zmin+1), ages[j], float(tauvals[l]), EBV, const[m, 0, zmin], chivals[m, zmin]])
    np.savetxt("photoz_grid.txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")


