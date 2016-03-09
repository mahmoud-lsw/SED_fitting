import numpy as np
import pylab
from astropy.cosmology import FlatLambdaCDM
SEDpath = "/disk1/adamc/progs/bc03/"
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)

# obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass

D = np.loadtxt("../IGM_Da_Db.txt")
model = np.loadtxt("photoz_2comp_UDS_8x8.txt")
model1 = np.loadtxt("photoz_2comp_UDS_1x5.txt")
fluxes = np.loadtxt("../../VANDELS_data/UDS_zphot_validation_phot.cat", usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
fluxerrs = np.expand_dims(np.loadtxt("../../VANDELS_data/UDS_zphot_validation_phot.cat", usecols=(13,14,15,16,17,18,19,20,21,22,23,24)), axis=2)
bandwavs = np.array([3776.4, 4426.9, 5454.6, 6506.3, 7645.3, 9010.9, 9093.0, 9183.0, 10200.1, 12483.2, 16318.0, 22010.2])


newSED_file = open(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII")
newages = np.array(newSED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
newagevals = np.array([69, 89, 107, 115, 123, 129, 132, 135])
newSED_file.close()
newspecdata = np.genfromtxt(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr
newwavs = np.copy(newspecdata[0, :])
newwavs = newwavs*10**-4
newcoef = np.copy(newwavs)
for j in range(6917):
    if newwavs[j] < 0.12:
        newcoef[j] = 2.659*(-2.156 + (1.509/newwavs[j]) - (0.198/(newwavs[j]**2)) + (0.011/(newwavs[j]**3))) + 4.05 - 69.44*(0.12-newwavs[j])
    elif newwavs[j] < 0.63:
        newcoef[j] = 2.659*(-2.156 + (1.509/newwavs[j]) - (0.198/(newwavs[j]**2)) + (0.011/(newwavs[j]**3))) + 4.05
    else:
        newcoef[j] = 2.659*(-1.857 + (1.040/newwavs[j])) + 4.05



oldSED_file = open(SEDpath+"models/Padova1994/chabrier/bc2003_hr_stelib_m62_chab_ssp.ised_ASCII")
oldages = np.array(oldSED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
oldagevals = np.array([129,132,135,138,144,149,152,156])
oldSED_file.close()
oldspecdata = np.genfromtxt(SEDpath+"models/Padova1994/chabrier/bc2003_hr_stelib_m62_chab_ssp.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr
oldwavs = np.copy(oldspecdata[0, :])
oldwavs = oldwavs*10**-4
oldcoef = np.copy(oldwavs)
for j in range(6917):
    if oldwavs[j] < 0.12:
        oldcoef[j] = 2.659*(-2.156 + (1.509/oldwavs[j]) - (0.198/(oldwavs[j]**2)) + (0.011/(oldwavs[j]**3))) + 4.05 - 69.44*(0.12-oldwavs[j])
    elif oldwavs[j] < 0.63:
        oldcoef[j] = 2.659*(-2.156 + (1.509/oldwavs[j]) - (0.198/(oldwavs[j]**2)) + (0.011/(oldwavs[j]**3))) + 4.05
    else:
        oldcoef[j] = 2.659*(-1.857 + (1.040/oldwavs[j])) + 4.05
        

for i in range(len(model)):

    EBV = model[i,6]
    norm = model[i,7]
    f_old_V = model[i,5]
    z = model[i,2]
    zind = int(z/0.01 - 1)
    oldageind = np.argmin(np.abs(model[i,3] - oldages))
    newageind = np.argmin(np.abs(model[i,4] - newages))
    print oldages[oldageind]
    print newages[newageind]

    newagenorms = np.loadtxt("../../models/const/agenormsUDS_" + str(newages[newageind]) + "_EBV_" + str(EBV) + ".txt")
    oldagenorms = np.loadtxt("../../models/burst/agenormsUDS_" + str(oldages[oldageind]) + "_EBV_" + str(EBV) + ".txt")

    newspectrum = np.array([newspecdata[0, :], newspecdata[newageind+1, :]]).T

    newzspectrum = np.copy(newspectrum)
    newzspectrum[:,1] = newzspectrum[:,1]*10**(-EBV*newcoef/2.5)

    for k in range(len(newzspectrum)):
        if newzspectrum[k,0] < 912.:
            newzspectrum[k,1] = 0.
        elif newzspectrum[k,0] > 912. and newzspectrum[k,0] < 1026.:
            newzspectrum[k,1] = newzspectrum[k,1]*(1-D[int(z/0.01 - 1),1])#*(1-D[int(z/0.01 - 1),0])
        elif newzspectrum[k,0] > 1026. and newzspectrum[k,0] < 1216.:
            newzspectrum[k,1] = newzspectrum[k,1]*(1-D[int(z/0.01 - 1),0])
                            
    newzspectrum[:,1] = newzspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
    newzspectrum[:,1] = newzspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z).value*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
    newzspectrum[:,1] = newzspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
    newzspectrum[:,1] = newzspectrum[:,1]/newagenorms[zind]
    newzspectrum[:,0] = newzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting

    oldspectrum = np.array([oldspecdata[0, :], oldspecdata[oldageind+1, :]]).T    
    oldzspectrum = np.copy(oldspectrum)
    oldzspectrum[:,1] = oldzspectrum[:,1]*10**(-EBV*oldcoef/2.5)

    for k in range(len(oldzspectrum)):
        if oldzspectrum[k,0] < 912.:
            oldzspectrum[k,1] = 0.
        elif oldzspectrum[k,0] > 912. and oldzspectrum[k,0] < 1026.:
            oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[int(z/0.01 - 1),1])#*(1-D[int(z/0.01 - 1),0])
        elif oldzspectrum[k,0] > 1026. and oldzspectrum[k,0] < 1216.:
            oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[int(z/0.01 - 1),0])
                            
    oldzspectrum[:,1] = oldzspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
    oldzspectrum[:,1] = oldzspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z).value*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
    oldzspectrum[:,1] = oldzspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
    oldzspectrum[:,1] = oldzspectrum[:,1]/oldagenorms[zind]
    oldzspectrum[:,0] = oldzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting

    if f_old_V == 1.:
        th_flux_array_new = 0.*newzspectrum[:,1]*norm
        th_flux_array_old = np.copy(oldzspectrum[:,1])*norm
    else:
        th_flux_array_new = np.copy(newzspectrum[:,1])*norm
        th_flux_array_old = (f_old_V/(1.-f_old_V))*oldzspectrum[:,1]*norm

    
    th_flux_array = (th_flux_array_new + th_flux_array_old)

    obs_flux_array = (10**-29)*fluxes[i,:]*(3*10**18/bandwavs/bandwavs)
    obs_fluxerr_array = (10**-29)*np.squeeze(fluxerrs[i,:])*(3*10**18/bandwavs/bandwavs)
    print "spec_z = " + str(model[i,1]) + ", phot_z = " + str(z) + ", alt phot_z = " + str(model1[i,2])
    pylab.figure()
    pylab.plot(oldzspectrum[:,0], th_flux_array_old, color="red", lw=0.5)
    pylab.plot(newzspectrum[:,0], th_flux_array_new, color="blue", lw=0.5)
    pylab.plot(newzspectrum[:,0], th_flux_array, color="black", lw=0.5)
    pylab.errorbar(bandwavs, obs_flux_array, yerr=obs_fluxerr_array, color="green", ls="none", zorder=10, lw=3, capsize=4)
    pylab.xlim(0, 22500)
    pylab.xlabel("Wavelength (A)", size=16)
    pylab.ylabel("f_lambda (erg/s/cm^2/A)", size=16)
    pylab.show()    

