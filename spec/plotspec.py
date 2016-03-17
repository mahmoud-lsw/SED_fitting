import numpy as np
import sys
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
import pylab

### Load up the galaxy catalogue and list of stellar population ages to be fit

all_obj_fluxes = np.expand_dims(np.loadtxt("../../VANDELS_data/combinedspec-0008.txt"), axis=2)
all_obj_fluxerrs =  np.expand_dims(np.loadtxt("../../VANDELS_data/combinederrspec-0008.txt"), axis=2) #erg/s/A/cm^2

objdata = np.loadtxt("photoz_V_2comp_1x5.txt") #../VANDELS_data/VANDELS_fitresults_v1.txt
photz = np.loadtxt("../../VANDELS_data/catalogue_0008.txt", dtype="str", usecols=())

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit

wavs1D = np.arange(0.500166449, 0.899967, 2*0.00025529999733)


oldages = np.array([508800000.0, 718700032.0, 1015200000, 1434000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0, 101520000.0, 255000000.0, 508800000.0, 718700032.0, 1015200000.0])


### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*9, dtype="float")
output.shape = (len(all_obj_fluxes), 9)
output[:,8] = 9999999999.

### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined

# obj_no phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass

param = objdata #np.loadtxt("photoz_VANDELS_v1.txt") #obj_no phot_z age_old age_new f_old_V old_modifier EBV norm chi
for m in range(len(param)):
    oldage = np.argmin(np.abs(oldages - param[m, 2]))
    newage = np.argmin(np.abs(newages - param[m, 3]))
    redshift_ind = (param[m, 1]/0.01)-1
    f_old_V = param[m, 4]
    EBV = param[m, 5]
    const = param[m, 6]


    th_flux_array_new_raw = np.loadtxt("../../models/spec/const/age_" + str(newages[newage]) + "_EBV_" + str(EBV) + ".txt")[:,redshift_ind]
    th_flux_array_old_raw = np.loadtxt("../../models/spec/burst/age_" + str(oldages[oldage]) + "_EBV_" + str(EBV) + ".txt")[:,redshift_ind] #erg/s/A/cm^2


    if f_old_V == 1.:
        th_flux_array_new = 0.*th_flux_array_new_raw
        th_flux_array_old = np.copy(th_flux_array_old_raw)
    else:
        th_flux_array_new = np.copy(th_flux_array_new_raw)
        th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw

    th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new), axis=0) #microjanskys
    
    pylab.figure()
    pylab.plot(wavs1D*10**4, np.squeeze(all_obj_fluxes[m,:,0]), color="blue")
    pylab.plot(wavs1D*10**4, np.squeeze(th_flux_array*const), color="red", zorder=10)
    pylab.plot(wavs1D*10**4, np.squeeze(all_obj_fluxerrs[m,:,0]), color="green")
    pylab.plot(wavs1D*10**4, np.zeros(len(wavs1D)), color="grey")
    #pylab.plot(binspec[:,0], np.squeeze(th_flux_array_oldonly*const), color="green", zorder=10)
    pylab.ylim(1.1*np.min(all_obj_fluxes[m,:,0]), 1.3*np.max(all_obj_fluxes[m,:,0]))
    pylab.xlabel("Wavelength (A)", size=16)
    pylab.ylabel("F_lambda (erg/s/A/cm^2)", size=16)
    pylab.text(5100, 0.8*np.max(all_obj_fluxes[m,:,0]),  "ACC_spec_z: " + str(param[m, 1]) + "   VANDELS_phot_z: " + photz[m,3] + "\n" + "E(B-V): " + str(EBV) + "   f_old: " + str(round(f_old_V, 3)) + "\n" + "old age: " + str(oldages[oldage]) + "  new age: " + str(newages[newage]) + "\n" + "Stellar Mass: " + str(round(param[m, 9]*10**-9, 2)) + "*10^9" + "\n" + "SFR: " + str(round(param[m, 8], 3)) + "\n" + "Reduced Chi-squared value: " + str(round(param[m, 7]/784., 3)), fontsize="14")
    #pylab.show()
    pylab.savefig("../../VANDELS_data/plots-0008/" + str(m+1) + ".png")
    pylab.close()

