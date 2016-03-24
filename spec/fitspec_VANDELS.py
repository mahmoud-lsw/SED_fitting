import numpy as np
import sys
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
import time
import pylab

foldmin = int(sys.argv[1])
foldmax = int(sys.argv[2])

### Load up the galaxy catalogue and list of stellar population ages to be fit

all_obj_fluxes = np.expand_dims(np.loadtxt("../../VANDELS_data/combinedspec-0008.txt"), axis=2)#[20:22,:,:]
all_obj_fluxerrs =  np.expand_dims(np.loadtxt("../../VANDELS_data/combinederrspec-0008.txt"), axis=2)#[20:22,:,:] #erg/s/A/cm^2

for i in range(len(all_obj_fluxes[:,0,0])):
    for j in range(len(all_obj_fluxes[0,:,0])):
        if all_obj_fluxerrs[i,j,0] == 0.:
            all_obj_fluxerrs[i,j,0] = 999999999.

### Build coefficient values to apply Calzetti et al. (2000) dust reddening law fit


"""            
for m in range(50):
    pylab.figure()
    pylab.plot(wavs*10**4, np.squeeze(all_obj_fluxes[m,:,0]), color="blue")
    pylab.plot(wavs*10**4, np.squeeze(all_obj_fluxerrs[m,:,0]), color="green")
    pylab.plot(wavs*10**4, np.zeros(len(wavs)), color="grey")
    #pylab.plot(binspec[:,0], np.squeeze(th_flux_array_oldonly*const), color="green", zorder=10)
    pylab.ylim(1.1*np.min(all_obj_fluxes[m,:,0]), 1.3*np.max(all_obj_fluxes[m,:,0]))
    pylab.xlabel("Wavelength (A)", size=16)
    pylab.ylabel("F_lambda (erg/s/A/cm^2)", size=16)
    pylab.show()
"""            
            
oldages = np.array([508800000.0, 718700032.0, 1015200000, 1434000000.0, 2000000000.0, 2500000000.0, 3000000000.0, 4000000000.0])
newages = np.array([10000000.0, 25119000.0, 50000000.0, 101520000.0, 255000000.0, 508800000.0, 718700032.0, 1015200000.0])

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

### Build array for photometric redshift and parameter outputs
output = np.zeros(len(all_obj_fluxes)*8, dtype="float")
output.shape = (len(all_obj_fluxes), 8)
output[:,7] = 9999999999.


### output = np.loadtxt("photoz_V_" + str(foldmin) + "_" + str(foldmax) + ".txt") ###########################################################################

zarr = np.arange(0.01, 7.01, 0.01)
lbtarr = cosmo.lookback_time(zarr).value

### Perform model fitting using a 2D chi squared array over object and redshift with max age of stellar pop phy sically determined

for j in range(1):
    arg = 0
    while oldages[j]*(10**-9) < 14.0 - lbtarr[arg] and arg < 699:
        arg = arg+1
    for k in range(41):
        EBV = 0.025*k
        print "Initial burst age: " + str(oldages[j]) + ", EBV: " + str(EBV)
        for l in range(5):
            th_flux_array_new_raw = np.loadtxt("../../models/spec/const/age_" + str(newages[l]) + "_EBV_" + str(EBV) + ".txt")[:,:arg+1]
            th_flux_array_old_raw = np.loadtxt("../../models/spec/burst/age_" + str(oldages[j]) + "_EBV_" + str(EBV) + ".txt")[:,:arg+1] #erg/s/A/cm^2
            for i in range(foldmin, foldmax+1):
                f_old_V = 0.02*i
                if f_old_V == 1.:
                    th_flux_array_new = 0.*th_flux_array_new_raw
                    th_flux_array_old = np.copy(th_flux_array_old_raw)
                else:
                    th_flux_array_new = np.copy(th_flux_array_new_raw)
                    th_flux_array_old = (f_old_V/(1.-f_old_V))*th_flux_array_old_raw

                th_flux_array = np.expand_dims((th_flux_array_old + th_flux_array_new), axis=0) #microjanskys
                const =  np.expand_dims(np.sum(all_obj_fluxes*th_flux_array/all_obj_fluxerrs**2, axis=1)/np.sum(th_flux_array**2/all_obj_fluxerrs**2, axis=1), axis=1)
                chivals = np.sum((all_obj_fluxes/all_obj_fluxerrs - const*(th_flux_array/all_obj_fluxerrs))**2, axis=1)
                for m in range(len(chivals)):
                    if np.min(chivals[m, :]) < output[m,7]:
                        zmin = np.argmin(chivals[m, :])
                        output[m,:] = np.array([m+1, 0.01*(zmin+1), oldages[j], newages[l], f_old_V, EBV, const[m, 0, zmin], chivals[m, zmin]])
                        """
                        if m == 0:
                            print output[m,:]
                            
                            pylab.figure()
                            pylab.plot(wavs1D*10**4, np.squeeze(all_obj_fluxes[m,:,0]), color="blue")
                            pylab.plot(wavs1D*10**4, np.squeeze(th_flux_array[:, :, zmin]*const[m, 0, zmin]), color="red", zorder=10)
                            pylab.plot(wavs1D*10**4, np.squeeze(all_obj_fluxerrs[m,:,0]), color="green")
                            pylab.plot(wavs1D*10**4, np.zeros(len(wavs1D)), color="grey")

                            pylab.plot(wavs1D*10**4, 0.1*np.squeeze(np.max(all_obj_fluxes[m,:,0])*(np.abs(((all_obj_fluxes[m,:,0] - th_flux_array[:, :, zmin]*const[m, 0, zmin])/all_obj_fluxerrs[m,:,0])))), color="black")
                            pylab.ylim(1.1*np.min(all_obj_fluxes[m,:,0]), 1.3*np.max(all_obj_fluxes[m,:,0]))
                            pylab.xlabel("Wavelength (A)", size=16)
                            pylab.ylabel("F_lambda (erg/s/A/cm^2)", size=16)
                            #pylab.text(5100, 0.8*np.max(all_obj_fluxes[m,:,0]),  "ACC_spec_z: " + str(param[m, 1]) + "   VANDELS_phot_z: " + photz[m,4] + "\n" + "E(B-V): " + str(EBV) + "   f_old: " + str(round(f_old_V, 3)) + "\n" + "old age: " + str(oldages[oldage]) + "  new age: " + str(newages[newage]) + "\n" + "Stellar Mass: " + str(round(param[m, 9]*10**-9, 2)) + "*10^9" + "\n" + "SFR: " + str(round(param[m, 8], 3)) + "\n" + "Reduced Chi-squared value: " + str(round(param[m, 7]/784., 3)), fontsize="14")
                            
                            chisq = np.sum(((all_obj_fluxes[m,:,0] - th_flux_array[:, :, zmin]*const[m, 0, zmin])/all_obj_fluxerrs[m,:,0])**2)
                            print chisq
                            
                            pylab.show()
                        """
    np.savetxt("photoz_V_" + str(foldmin) + "_" + str(foldmax) + ".txt", output, header="obj_no phot_z age_old age_new f_old_V  EBV norm chi")
    
