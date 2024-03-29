import numpy as np
import sys
from astropy.cosmology import WMAP9

tval = int(sys.argv[1])

zarr = np.arange(0.01, 7.01, 0.01)

### Load up the galaxy catalogue and list of stellar population ages to be fit
all_obj_specz = np.loadtxt("../../uds_validation_results_all.cat", usecols=(1,))
all_obj_fluxes = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_validation_phot.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2)
all_obj_fluxerrs = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_validation_phot.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)

SEDpath = "/disk1/adamc/progs/bc03/"
SED_file = open(SEDpath+"P94_Chab_zsun_t1.71Gyr_noabs/P94_Chab_zsun_t1.71Gyr_noabs.ised_ASCII")
ages = np.array(SED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
SED_file.close()

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
"""

### For all objects which have not been observed in a given band, blow up the associated error so fit is not affected
for i in range(len(all_obj_fluxes)):
    for j in range(12):
        if all_obj_fluxes[i, j, 0] == 0.:
            all_obj_fluxerrs[i, j, 0] = 9999999999.

### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*8, dtype="float")
output.shape = (len(all_obj_fluxes), 8)
output[:,7] = 9999999999.


### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined
tauvals = ["0.03", "0.05", "0.09", "0.16", "0.29", "0.53", "0.95", "1.71", "3.08", "5.55", "10.00", "18.01"]
lbtarr = WMAP9.lookback_time(zarr).value
    
for j in range(69, 195, 3): 
    arg = 0
    while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 699:
        arg = arg+1
    for k in range(41):
        EBV = 0.025*k
        print "Age: " + str(ages[j]) + ", fitting to redshift " + str((arg+1)*0.01) + ", EBV: " + str(EBV)
        for l in range(tval, tval+3):
            th_mag_array = np.loadtxt("../../models/tau/" + tauvals[l] + "/synmags_age_" + str(ages[j]) + "_EBV_" + str(EBV) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[0:arg+1, :]
            th_flux_array = np.expand_dims((10**((23.9 - th_mag_array)/2.5)).T, axis=0) #microjanskys
            const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
            chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
            for m in range(len(chivals)):
                if np.min(chivals[m, :]) < output[m,7]:
                    zmin = np.argmin(chivals[m, :])
                    output[m,:] = np.array([m+1, all_obj_specz[m], 0.01*(zmin+1), ages[j], float(tauvals[l]), EBV, const[m, 0, zmin], chivals[m, zmin]])
    np.savetxt("photoz_tau_" + tauvals[tval] + "_to_" + tauvals[tval+2] + ".txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")

