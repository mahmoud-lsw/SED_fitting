import numpy as np
from astropy.cosmology import FlatLambdaCDM
import sys
sys.path.append("/disk1/adamc/my_code")
import synphot as s
tval = int(sys.argv[1])
SEDpath = "/disk1/adamc/progs/bc03/"
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)


tauvals = ["0.03", "0.05", "0.09", "0.16", "0.29", "0.53", "0.95", "1.71", "3.08", "5.55", "10.00", "18.01"]
D = np.loadtxt("../IGM_Da_Db.txt")


for t in range(tval, tval+1):
    SED_file = open(SEDpath+"P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs/P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs.ised_ASCII")
    ages = np.array(SED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
    SED_file.close()

    specdata = np.genfromtxt(SEDpath+"P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs/P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr

    wavs = np.copy(specdata[0, :])
    wavs = wavs*10**-4
    coef = np.copy(wavs)
    for j in range(6917):
        if wavs[j] < 0.12:
            coef[j] = 2.659*(-2.156 + (1.509/wavs[j]) - (0.198/(wavs[j]**2)) + (0.011/(wavs[j]**3))) + 4.05 - 69.44*(0.12-wavs[j])
        elif wavs[j] < 0.63:
            coef[j] = 2.659*(-2.156 + (1.509/wavs[j]) - (0.198/(wavs[j]**2)) + (0.011/(wavs[j]**3))) + 4.05
        else:
            coef[j] = 2.659*(-1.857 + (1.040/wavs[j])) + 4.05
    
    for e in range(41):
        EBV = 0.025*e
        for j in range(69, 195, 3):
     
            spectrum = np.array([specdata[0, :], specdata[j+1, :]]).T
            spectrum[:,1] = spectrum[:,1]*10**(-EBV*coef/2.5)
            synmags = np.zeros(13*700, dtype="float")
            synmags.shape = (700, 13)
            synmags[:,0] = np.arange(0.01, 7.01, 0.01)

            for i in range(1, 701):
                z = 0.01*i
                if ages[j]*(10**-9) < 14.00 - cosmo.lookback_time(z).value:
                    print "Calculating colours for age: " + str(ages[j]) + " and redshift: " + str(z) + ", total age: " + str(ages[j]*(10**-9) + cosmo.lookback_time(z).value) + ", tau = " + tauvals[t]
                    zspectrum = np.copy(spectrum)
                    for k in range(len(zspectrum)):
                        if zspectrum[k,0] < 912.:
                            zspectrum[k,1] = 0.
                        elif zspectrum[k,0] > 912. and zspectrum[k,0] < 1026.:
                            zspectrum[k,1] = zspectrum[k,1]*(1-D[i-1,1])
                        elif zspectrum[k,0] > 1026. and zspectrum[k,0] < 1216.:
                            zspectrum[k,1] = zspectrum[k,1]*(1-D[i-1,0])
                    zspectrum[:,1] = zspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
                    zspectrum[:,1] = zspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z)*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
                    zspectrum[:,1] = zspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
                    zspectrum[:,0] = zspectrum[:,0]*(1+z) #change wavelength values to account for redshifting

                    th_mags = s.AB_mags_UDS(zspectrum)
                    synmags[i-1, 1:] = th_mags

            np.savetxt("../../models/tau/" + tauvals[t] + "/synmags_age_" + str(ages[j]) + "_EBV_" + str(EBV) + ".txt", synmags)
