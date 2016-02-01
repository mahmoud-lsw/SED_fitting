import sys
sys.path.append("../my_code")
SEDpath = "../progs/bc03/"

import numpy as np
from subprocess import call
import synphot as s
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
import pylab

oldSED_file = open(SEDpath+"P94_Chab_zsun_const/P94_Chab_zsun_const.ised_ASCII")
oldages = np.array(oldSED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
oldSED_file.close()

D = np.loadtxt("IGM_Da_Db.txt")

oldspecdata = np.genfromtxt(SEDpath+"P94_Chab_zsun_const/P94_Chab_zsun_const.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr
oldspectrum = np.array([oldspecdata[0, :], oldspecdata[125, :]]).T
oldsynmags = np.zeros(14*6, dtype="float")
oldsynmags.shape = (6, 14)
oldsynmags[:,0] = np.array([2., 4., 6., 2., 4., 6.])

for j in range(2):
    for i in range(1, 4):
        z = 2*i
        oldzspectrum = np.copy(oldspectrum)
        oldzspectrum[:,1]*(10**10)/oldages[124]
        if j == 1:
            EBV = 0.247
            wavs1D = np.array([3722.9, 4299.0, 5844.1, 7667.7, 8003.1, 9001.8, 10483.0, 12425.7, 13830.1, 15324.7, 21425.2, 35378.4, 44780.0])*10**-4
            zarr = oldsynmags[0:3,0]
            wavs = np.zeros(len(wavs1D)*3)
            wavs.shape = (len(wavs1D), 3)

            for k in range(3):
                wavs[:,k] = wavs1D/(1+zarr[k])

            coef = np.copy(wavs)
            for l in range(3):
                for m in range(len(wavs1D)):
                    if wavs[m,l] < 0.63:
                        coef[m, l] = 2.659*(-2.156 + (1.509/wavs[m, l]) - (0.198/(wavs[m, l]**2)) + (0.011/(wavs[m, l]**3))) + 4.05 
                    else:
                        coef[m, l] = 2.659*(-1.857 + (1.040/wavs[m, l])) + 4.05 

        for k in range(len(oldzspectrum)):
            if oldzspectrum[k,0] < 912.:
                oldzspectrum[k,1] = 0.
            elif oldzspectrum[k,0] > 912. and oldzspectrum[k,0] < 1026.:
                oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[int(z/0.01 - 1),1])*(1-D[int(z/0.01 - 1),0])
            elif oldzspectrum[k,0] > 1026. and oldzspectrum[k,0] < 1216.:
                oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[int(z/0.01 - 1),0])
        oldzspectrum[:,1] = oldzspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
        oldzspectrum[:,1] = oldzspectrum[:,1]*10**10/(4*np.pi*(cosmo.luminosity_distance(z).value*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
        oldzspectrum[:,1] = oldzspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
        oldzspectrum[:,0] = oldzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting
        oldth_mags = s.AB_mags_Ross(oldzspectrum)
        if j == 1:
            oldth_mags = oldth_mags + EBV*coef[:,i-1] #microjanskys
        oldsynmags[3*j + i-1, 1:] = oldth_mags
                

    np.savetxt("300Myrburst.txt", oldsynmags)



