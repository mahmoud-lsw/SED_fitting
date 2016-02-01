import numpy as np
import sys
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
import pylab

### Load up the galaxy catalogue and list of stellar population ages to be fit

all_obj_fluxes = np.expand_dims(np.loadtxt("../VANDELS_data/combined_spec.txt"), axis=2)
all_obj_fluxerrs =  np.expand_dims(np.loadtxt("../VANDELS_data/combined_spec_errs.txt"), axis=2) #erg/s/A/cm^2

objdata = np.loadtxt("photoz_V_2comp.txt") #../VANDELS_data/VANDELS_fitresults_v1.txt
photz = np.loadtxt("../VANDELS_data/VANDELS_fitresults_v1.txt", dtype="str", usecols=())

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit

wavs1D = np.loadtxt("../VANDELS_data/wavs.txt")*10**-4
zarr = np.arange(0.01, 7.01, 0.01)
wavs = np.zeros(len(wavs1D)*700)
wavs.shape = (len(wavs1D), 700)

for i in range(700):
    wavs[:,i] = wavs1D/(1+zarr[i])

coef = np.copy(wavs)
for i in range(700):
    for j in range(784):
        if wavs[j,i] < 0.63:
            coef[j, i] = ((0.013/(wavs[j,i]**3)) - (0.232/(wavs[j,i]**2)) + (1.766/wavs[j,i]) - 0.743)
        else:
            coef[j, i] = 1.217/wavs[j,i] - 0.393


oldages = np.array([508800000.0, 1139100030.0, 1700000000.0, 2000000000.0, 4000000000.0]) #508800000.0, 806400000.0, 1139100030.0, 1434000000.0, 1700000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0########
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
output = np.zeros(len(all_obj_fluxes)*9, dtype="float")
output.shape = (len(all_obj_fluxes), 9)
output[:,8] = 9999999999.

lbtarr = cosmo.lookback_time(zarr).value
f_old_array = 1.01 - np.logspace(-2, 0, 51)

### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined

param = objdata #np.loadtxt("photoz_VANDELS_v1.txt") #obj_no phot_z age_old age_new f_old_V old_modifier EBV norm chi
for m in range(len(param)):
    print m
    oldage = np.argmin(np.abs(oldages - param[m, 2]))
    newage = np.argmin(np.abs(newages - param[m, 3]))
    redshift_ind = (param[m, 1]/0.01)-1
    f_old_V = param[m, 4]
    EBV = param[m, 5]
    const = param[m, 6]

    """
    print "z: " + str(param[m, 1])
    print "old age: " + str(oldages[oldage]) + " new age: " + str(newages[newage])
    print "E(B-V): " + str(EBV)
    print "f_old: " + str(f_old_V)
    print "const: " + str(const)
    print "Stellar Mass: " + str(round(param[m, 9]*10**-9,3)) + "*10^9 Solar masses"
    print "SFR: " + str(param[m, 8])
    print "Reduced Chi-squared value: " + str(param[m, 7]/784.)
    """

    th_flux_array_new_raw = np.loadtxt("../models/spec/newconst/age_" + str(newages[newage]) + ".txt")[:,redshift_ind:redshift_ind+1]#[:,:arg+1]
    th_flux_array_old_raw = np.loadtxt("../models/spec/oldburst/age_" + str(oldages[oldage]) + ".txt")[:,redshift_ind:redshift_ind+1]#[:,:arg+1] #erg/s/A/cm^2

    if f_old_V == 1.:
        th_flux_array_new = 0.*th_flux_array_new_raw
        th_flux_array_old = np.copy(th_flux_array_old_raw)
    else:
        th_flux_array_new = np.copy(th_flux_array_new_raw)
        th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw

    th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new)*(10**((-EBV*coef[:,redshift_ind:redshift_ind+1])/2.5)), axis=0) #microjanskys
    
    pylab.figure()
    pylab.plot(wavs1D*10**4, np.squeeze(all_obj_fluxes[m,:,0]), color="blue")
    pylab.plot(wavs1D*10**4, np.squeeze(th_flux_array*const), color="red", zorder=10)
    pylab.plot(wavs1D*10**4, np.squeeze(all_obj_fluxerrs[m,:,0]), color="green")
    pylab.plot(wavs1D*10**4, np.zeros(len(wavs1D)), color="grey")
    #pylab.plot(binspec[:,0], np.squeeze(th_flux_array_oldonly*const), color="green", zorder=10)
    pylab.ylim(1.1*np.min(all_obj_fluxes[m,:,0]), 1.3*np.max(all_obj_fluxes[m,:,0]))
    pylab.xlabel("Wavelength (A)", size=16)
    pylab.ylabel("F_lambda (erg/s/A/cm^2)", size=16)
    pylab.text(5100, 0.8*np.max(all_obj_fluxes[m,:,0]),  "ACC_spec_z: " + str(param[m, 1]) + "   VANDELS_phot_z: " + photz[m,4] + "\n" + "E(B-V): " + str(EBV) + "   f_old: " + str(round(f_old_V, 3)) + "\n" + "old age: " + str(oldages[oldage]) + "  new age: " + str(newages[newage]) + "\n" + "Stellar Mass: " + str(round(param[m, 9]*10**-9, 2)) + "*10^9" + "\n" + "SFR: " + str(round(param[m, 8], 3)) + "\n" + "Reduced Chi-squared value: " + str(round(param[m, 7]/784., 3)), fontsize="14")
    #pylab.show()
    pylab.savefig("../VANDELS_data/plots/" + str(m+1) + ".png")
    pylab.close()

