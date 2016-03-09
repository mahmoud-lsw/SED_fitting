# l_bfgs_b fitting with approx. gradient. The process by which EBV and Norm are saved is currently wrong.

import numpy as np
import pylab
import scipy.optimize as opt
from astropy.cosmology import WMAP9
import astropy.units as u
from subprocess import call

ages = np.loadtxt("ages.txt")

wavs = np.array([3765.1, 4331.5, 5778.3, 7624.6, 7988.1, 9012.0, 10364.5, 12304.9, 15235.9, 21322.8, 35075.1, 44365.8])
wavs = wavs*10**-4
coef = np.copy(wavs)
coef[:3] = ((0.013/(wavs[:3]**3)) - (0.232/(wavs[:3]**2)) + (1.766/wavs[:3]) - 0.743)
coef[3:] = 1.217/wavs[3:] - 0.393

def chisqfunc(param, th_mags, obj_fluxes, obj_fluxerrs):    #param0 = norm, param1 = EBV 
    if param[1] < 0:
        return 999999.
    th_fluxes = param[0]*10**((23.9 - param[1]*coef - th_mags)/2.5) #microjanskys
    chisq = np.sum(((obj_fluxes - th_fluxes)/obj_fluxerrs)**2)
    return chisq

all_obj_specz = np.loadtxt("galaxy_photom_specz_adam.cat", usecols=(27,))
all_obj_fluxes = np.loadtxt("galaxy_photom_specz_adam.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
all_obj_fluxerrs = np.loadtxt("galaxy_photom_specz_adam.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24))

output = np.zeros(len(all_obj_specz)*8, dtype="float")
output.shape = (len(all_obj_specz), 8)

agesteps = 127
zsteps = 500

tauvals = np.array([0.05, 1, 10], dtype="str")

for k in range(len(all_obj_fluxes)):
    minchis = np.zeros(3, dtype="float")
    obj_fluxes = all_obj_fluxes[k,:]
    obj_fluxerrs = all_obj_fluxerrs[k,:]
    if np.min(obj_fluxes) > 0.:
        for t in range(3):
            chi_vals = np.zeros(agesteps*zsteps, dtype="float")
            chi_vals.shape = (zsteps, agesteps)
            norm_vals = np.copy(chi_vals)
            EBV_vals = np.copy(chi_vals)
            fail_flag = np.copy(chi_vals)

            for j in range(127):
                print "calculating chi squared values for object " + str(k+1)+ ", with tau = " + tauvals[t] + " and age = " + str(ages[j])

                th_mag_array = np.loadtxt("synmags_T" + tauvals[t] + "/synmags/synmags_age_" + str(ages[j]) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
                
                for i in range(1, 501):
                    z = 0.01*i
                    if ages[j]*(10**-9) < 14.00 - WMAP9.lookback_time(z).value:

                        th_mags = th_mag_array[i-1,:]
                        th_fluxes = 10**((23.9 - th_mags)/2.5) #microjanskys
                        x0 = np.mean(obj_fluxes)/np.mean(th_fluxes)
                        th_fluxes_scaled = th_fluxes*x0
                        th_mags_scaled = -2.5*np.log10(th_fluxes_scaled) + 23.9
                        #optresult = opt.minimize(chisqfunc, [1., 1.], args=(th_mags, obj_fluxes, obj_fluxerrs), method="Nelder-Mead", options={"maxfev":9999, "maxiter":9999})
                        optresult = opt.fmin_l_bfgs_b(chisqfunc, [1., 1.], args=(th_mags_scaled, obj_fluxes, obj_fluxerrs), bounds=[(-np.infty, np.infty), (0, np.infty)], approx_grad=True, pgtol=1e-10) #change to th_fluxes_scaled

                        if optresult[2]["warnflag"] != 0:
                            print "Optimisation failed"
                            #print optresult
                            fail_flag[i-1, j] = 1.
                            chi_vals[i-1, j] = 9999999999.
                            norm_vals[i-1, j] = 9999999999.
                            EBV_vals[i-1, j] = 9999999999.

                        else:
                            chi_vals[i-1, j] = optresult[1]
                            norm_vals[i-1, j] = optresult[0][0]*(1./x0)*10**5
                            EBV_vals[i-1, j] = optresult[0][1]*10**5

                    else:
                        chi_vals[i-1, j] = 9999999999.
                        norm_vals[i-1, j] = 9999999999.
                        EBV_vals[i-1, j] = 9999999999.
            minchis[t] = np.min(chi_vals)
            np.savetxt("synmags_T" + tauvals[t] + "/failflagarrays/obj_" + str(k+1) + ".txt", fail_flag)
            np.savetxt("synmags_T" + tauvals[t] + "/chiarrays/obj_" + str(k+1) + ".txt", chi_vals)
            np.savetxt("synmags_T" + tauvals[t] + "/normarrays/obj_" + str(k+1) + ".txt", norm_vals)
            np.savetxt("synmags_T" + tauvals[t] + "/EBVarrays/obj_" + str(k+1) + ".txt", EBV_vals)

        best_chi_array = np.loadtxt("synmags_T" + tauvals[np.argmin(minchis)] + "/chiarrays/obj_" + str(k+1) + ".txt")
        x, y = np.unravel_index(best_chi_array.argmin(), best_chi_array.shape)
        best_norm_array = np.loadtxt("synmags_T" + tauvals[np.argmin(minchis)] + "/normarrays/obj_" + str(k+1) + ".txt")
        best_EBV_array = np.loadtxt("synmags_T" + tauvals[np.argmin(minchis)] + "/EBVarrays/obj_" + str(k+1) + ".txt")
        output[k,:] = np.array([k+1, all_obj_specz[k], 0.01*(x+1), ages[y], float(tauvals[np.argmin(minchis)]), best_EBV_array[x, y], best_norm_array[x, y], best_chi_array[x, y]])
        np.savetxt("photoz_v3.txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")

