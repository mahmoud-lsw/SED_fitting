import numpy as np
import sys
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)

foldmin = int(sys.argv[1])
foldmax = int(sys.argv[2])

### Load up the galaxy catalogue and list of stellar population ages to be fit
all_obj_specz = np.loadtxt("../../uds_validation_results_all.cat", usecols=(1,))#np.loadtxt("UDS_zphot_validation_phot.cat", usecols=(25,))
all_obj_fluxes = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_validation_phot.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12)), axis=2)
all_obj_fluxerrs = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_validation_phot.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)

ages = np.array([508800000.0, 718700032.0, 1015200000, 1434000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0, 101520000.0, 255000000.0, 508800000.0, 718700032.0, 1015200000.0])


"""
# Apply offsets to the data as calculated by code in zeropoints/ in order to deal with various systematic calibration errors
offsets = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/mean_ratios_2comp_UDS.txt", usecols=(0,)), axis=0), axis=2)
offset_errs = np.expand_dims(np.expand_dims(np.loadtxt("zeropoints/mean_ratios_2comp_UDS.txt", usecols=(1,)), axis=0), axis=2)

for i in range(12): # If the offset is consistent with zero to within 1 sigma do nothing
    if np.abs(offsets[0, i, 0] - 1.) > 0.:# offset_errs[0, i, 0]:
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
output = np.zeros(len(all_obj_fluxes)*9, dtype="float")
output.shape = (len(all_obj_fluxes), 9)
output[:,8] = 9999999999.


### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop physically determined
zarr = np.arange(0.01, 7.01, 0.01)
lbtarr = cosmo.lookback_time(zarr).value
f_old_array = np.arange(0., 1.001, 0.02)#1.01 - np.logspace(-2, 0, 51)


for j in range(3):
    arg = 0
    while ages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 699:
        arg = arg+1
    for l in range(5):
        for k in range(41):
            EBV = 0.025*k
            print "Burst age: " + str(ages[j]) + ", fitting to redshift " + str((arg+1)*0.01) + ", const age " + str(newages[l]) + " and EBV " + str(EBV)
            th_mag_array_old = np.loadtxt("../../models/burst/synmagsUDS_age_" + str(ages[j]) + "_EBV_" + str(EBV) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[:arg+1, :]
            th_mag_array_new = np.loadtxt("../../models/const/synmagsUDS_age_" + str(newages[l]) + "_EBV_" + str(EBV) + ".txt", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))[:arg+1, :]
            th_flux_array_old_raw = 10**((23.9-th_mag_array_old)/2.5)
            th_flux_array_new_raw = 10**((23.9-th_mag_array_new)/2.5)
            
            for i in range(foldmin, foldmax+1):
                f_old_V = f_old_array[i]

                if f_old_V == 1.:
                    th_flux_array_new = 0.*th_flux_array_new_raw
                    th_flux_array_old = np.copy(th_flux_array_old_raw)
                else:
                    th_flux_array_new = np.copy(th_flux_array_new_raw)
                    th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw

                th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new).T, axis=0) #microjanskys
                const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
                chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
                for m in range(len(chivals)):
                    if np.min(chivals[m, :]) < output[m,8]:
                        zmin = np.argmin(chivals[m, :])
                        output[m,:] = np.array([m+1, all_obj_specz[m], 0.01*(zmin+1), ages[j], newages[l], f_old_V, EBV, const[m, 0, zmin], chivals[m, zmin]])
                        
    np.savetxt("photoz_UDS_" + str(foldmin) + "_" + str(foldmax) + ".txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi")
    

