#call(["csp", "$bc03/models/Padova1994/chabrier/bc2003_lr_BaSeL_m62_chab_ssp.ised", )
#stellar_mass = 7.6128*10**9 at 13Gyr (ages[-53])

import sys
sys.path.append("/disk1/adamc/my_code")

SEDpath = "/disk1/adamc/progs/bc03/"

import numpy as np
from subprocess import call
import synphot as s
from astropy.cosmology import WMAP9 as cosmo

tauvals = np.array([0.05, 1, 10], dtype="str")


for t in range(3):
    SED_file = open(SEDpath+"P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs/out.ised_ASCII")
    ages = np.array(SED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
    SED_file.close()
    np.savetxt("ages.txt", ages[69:197])
    specdata = np.genfromtxt(SEDpath+"P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs/out.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr

    for j in range(181, 197):
 
        spectrum = np.array([specdata[0, :], specdata[j+1, :]]).T

        synmags = np.zeros(13*500, dtype="float")
        synmags.shape = (500, 13)
        synmags[:,0] = np.arange(0.01, 5.01, 0.01)

        for i in range(1, 501):
            z = 0.01*i
            if ages[j]*(10**-9) < 14.00 - cosmo.lookback_time(z):
                print "Calculating colours for age: " + str(ages[j]) + " and redshift: " + str(z) + ", total age: " + str(ages[j]*(10**-9) + cosmo.lookback_time(z))
                zspectrum = np.copy(spectrum)
                zspectrum[:,1] = zspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
                zspectrum[:,1] = zspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z)*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
                zspectrum[:,1] = zspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
                zspectrum[:,0] = zspectrum[:,0]*(1+z) #change wavelength values to account for redshifting

                th_mags = s.AB_mags_Guo(zspectrum)
                synmags[i-1, 1:] = th_mags

        np.savetxt("synmags_T" + tauvals[t] + "/synmags/synmags_age_" + str(ages[j]) + ".txt", synmags)

