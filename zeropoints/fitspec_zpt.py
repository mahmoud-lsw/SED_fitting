import numpy as np
import sys
from astropy.cosmology import WMAP9

firstobj = int(sys.argv[1])
lastobj = int(sys.argv[2])

ages = np.loadtxt("../ages.txt")

wavs = np.array([3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8])
wavs = wavs*10**-4
coef = np.expand_dims(np.copy(wavs), axis=0)
coef[0, :3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[0, 3:] = 1.217/wavs[3:] - 0.393

all_obj_specz = np.loadtxt("../galaxy_photom_specz_adam_zless5.cat", usecols=(27,))
all_obj_fluxes = np.expand_dims(np.loadtxt("../galaxy_photom_specz_adam_zless5.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2)
all_obj_fluxerrs = np.expand_dims(np.loadtxt("../galaxy_photom_specz_adam_zless5.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)

for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0] == 0.:
            all_obj_fluxerrs[i, j, 0] = 9999999999.

output = np.zeros((lastobj-firstobj)*8, dtype="float")
output.shape = ((lastobj-firstobj), 8)
output[:,7] = 9999999999.

tauvals = np.array([0.05, 1, 10], dtype="str")

zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = WMAP9.lookback_time(zarr).value

flux_diffs = np.copy(all_obj_fluxes[firstobj:lastobj+1, :, 0])

EBV = np.expand_dims(np.arange(0, 1.50001, 0.025), axis=1)

for m in range(firstobj, lastobj+1):
    print "object " + str(m)
    
    obj_fluxes = all_obj_fluxes[m, :, :].T
    obj_fluxerrs = all_obj_fluxerrs[m, :, :].T
    z = all_obj_specz[m]
    modelno = int((np.round(z, 2)*100.) - 1.)
    for j in range(127): 
        arg = 0
        while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 499:
            arg = arg+1
        for l in range(3):
            th_mag_array = np.expand_dims(np.loadtxt("../models/" + tauvals[l] + "/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[modelno, :], axis=0)
            th_flux_array = (10**((23.9 - EBV*coef - th_mag_array)/2.5)) #microjanskys
            const =  np.expand_dims(np.sum(obj_fluxes*th_flux_array/obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/obj_fluxerrs**2, axis=1), axis=1)
            chivals = np.sum((obj_fluxes/obj_fluxerrs - const*(th_flux_array/obj_fluxerrs))**2, axis=1)
            if np.min(chivals) < output[m-firstobj,7]:
                minchi = np.argmin(chivals)
                output[m-firstobj,:] = np.array([m+1, all_obj_specz[m], z, ages[j], float(tauvals[l]), EBV[minchi], const[minchi, 0], chivals[minchi]])
                flux_diffs[m-firstobj, :] = obj_fluxes/(th_flux_array[minchi, :]*const[minchi, 0])
    #np.savetxt("photoz_offsetcalc.txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")
    
    np.savetxt("flux_fracoffsets_" + str(firstobj) + "_" + str(lastobj) + ".txt", flux_diffs)

