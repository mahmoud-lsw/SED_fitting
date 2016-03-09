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
fulldet_obj_fluxes.shape = (len(ind), 12)
fulldet_obj_fluxerrs = np.copy(fulldet_obj_fluxes)

for j in range(len(fulldet_obj_fluxes)):
    fulldet_obj_fluxes[j, :] = all_obj_fluxes[ind[j],:]
    fulldet_obj_fluxerrs[j, :] = all_obj_fluxerrs[ind[j],:]

output = np.zeros(len(all_obj_specz)*8, dtype="float")
output.shape = (len(all_obj_specz), 8)
output[:,7] = 9999999999.

tauvals = np.array([0.05, 1, 10], dtype="str")

#individual chi squared fits to each object (very slow)

"""
for j in range(128): 
    for l in range(3):
        th_mag_array = np.loadtxt("synmags_T" + tauvals[l] + "/synmags/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
        for i in range(1, 501): 
            z = 0.01*i
            th_mags = th_mag_array[i-1,:]
            if ages[j]*(10**-9) < 14.00 - WMAP9.lookback_time(z).value:
                for k in range(41): #more steps?
                    EBV = 0.05*k #higher density?
                    th_fluxes = 10**((23.9 - EBV*coef - th_mags)/2.5) #microjanskys
                    for m in range(len(all_obj_fluxes)):
                        obj_fluxes = all_obj_fluxes[m,:]
                        obj_fluxerrs = all_obj_fluxerrs[m,:]
                        if np.min(obj_fluxes) > 0.:
                            const =  np.sum(obj_fluxes*th_fluxes/obj_fluxerrs**2)/np.sum(th_fluxes**2/obj_fluxerrs**2)
                            chival = np.sum((obj_fluxes/obj_fluxerrs - const*th_fluxes/obj_fluxerrs)**2)
                            print chival
                            raw_input()
                            if chival < output[m,7]:
                                output[m,:] = np.array([m+1, all_obj_specz[m], z, ages[j], float(tauvals[l]), EBV, const, chival])
                print "Age: " + str(ages[j]) + ", tau: " + tauvals[l] + ", redshift: " + str(z)
                np.savetxt("photoz_grid84876.txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")

"""
#the chi squared array is now 2D over EBV and the different objects

for j in range(128): 
    for l in range(3):
        th_mag_array = np.loadtxt("synmags_T" + tauvals[l] + "/synmags/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
        for i in range(1, 501): 
            z = 0.01*i
            th_mags = th_mag_array[i-1,:]
            if ages[j]*(10**-9) < 14.00 - WMAP9.lookback_time(z).value:
                for k in range(41):
                    EBV = 0.05*k
                    th_fluxes = 10**((23.9 - EBV*coef - th_mags)/2.5) #microjanskys
                    const =  np.sum(fulldet_obj_fluxes*th_fluxes/fulldet_obj_fluxerrs**2, axis=1)/np.sum(th_fluxes**2/fulldet_obj_fluxerrs**2, axis=1)
                    chival = np.sum((fulldet_obj_fluxes/fulldet_obj_fluxerrs - (const*(th_fluxes/fulldet_obj_fluxerrs).T).T)**2, axis=1)
                    for m in range(len(fulldet_obj_fluxes)):
                        if chival[m] < output[m,7]:
                            output[m,:] = np.array([m+1, all_obj_specz[m], z, ages[j], float(tauvals[l]), EBV, const[m], chival[m]])
                print "Age: " + str(ages[j]) + ", tau: " + tauvals[l] + ", redshift: " + str(z)
                np.savetxt("photoz_grid86548.txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")


