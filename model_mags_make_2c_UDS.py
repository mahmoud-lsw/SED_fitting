#call(["csp", "$bc03/models/Padova1994/chabrier/bc2003_lr_BaSeL_m62_chab_ssp.ised", )
#stellar_mass = 7.6128*10**9 at 13Gyr (ages[-53])

import sys
sys.path.append("../my_code")
arg = int(sys.argv[1])
SEDpath = "../progs/bc03/"

import numpy as np
from subprocess import call
import synphot_UDS as s
from astropy.cosmology import WMAP9 as cosmo
import pylab



"""bins up spectrum consisting of a column of wavelength values, a column of fluxes and a column of flux errors by factor binn"""
def specbin(spec, binn): 
    if int(2*len(spec)/binn) != 2*int(len(spec)/binn):
        binned = np.zeros(len(spec)/binn -1)
    else:
        binned = np.zeros(len(spec)/binn)
    
    for i in range(len(binned)):
        binned[i] = np.mean(spec[binn*i:binn*i+binn])
    return binned



"""Takes bin midpoints and creates bin LHS and bin width arrays, first bin starts at the first midpoint minus half the distance between the first two midpoints. If make_rhs is "True", the RHS of the final bin is included as the last entry in the bin_lhs array."""
def make_bins(wavelengths, make_rhs="False"):
    bin_widths = np.zeros(len(wavelengths), dtype="float")
    if make_rhs == "True":
        bin_lhs = np.zeros(len(wavelengths)+1, dtype="float")
    else:
        bin_lhs = np.zeros(len(wavelengths), dtype="float")

    bin_lhs[0] = wavelengths[0] - (wavelengths[1]-wavelengths[0])/2
    if make_rhs == "True":
        bin_lhs[-1] = wavelengths[-1] + (wavelengths[-1]-wavelengths[-2])/2

    for i in range(1, len(wavelengths)):
        bin_lhs[i] = (wavelengths[i] + wavelengths[i-1])/2
        bin_widths[i-1] = bin_lhs[i]-bin_lhs[i-1]
    bin_widths[-1] = (wavelengths[-1] - wavelengths[-2])
    return bin_lhs, bin_widths



def rebin_to_filter(spec, filter):
    spec_lhs, spec_widths = make_bins(spec[:,0])
    filter_lhs, filter_widths = make_bins(filter[:,0], make_rhs="True")
    i=0
    i_bins_in_f_bin = []
    rebinned = []

    for j in range(len(filter_lhs)-1):
        i_bins_in_f_bin = []
        while spec_lhs[i] < filter_lhs[j]:
            i = i+1

        if spec_lhs[i] != filter_lhs[j]:
            i_bins_in_f_bin.append(spec_lhs[i-1])
            bin_start = i-1
        else:
            bin_start = i

        while spec_lhs[i] < filter_lhs[j+1]:
            i_bins_in_f_bin.append(spec_lhs[i])
            i=i+1

        i = i-1
        bin_end = i
        perc_under_f_bin = np.ones(len(i_bins_in_f_bin), dtype="float")
    
        if len(i_bins_in_f_bin) == 1:
            perc_under_f_bin[0] = (filter_lhs[j+1] - filter_lhs[j])/(spec_lhs[bin_start+1] - spec_lhs[bin_start])
        else:
            perc_under_f_bin[0] = (spec_lhs[bin_start+1] - filter_lhs[j])/(spec_lhs[bin_start+1] - spec_lhs[bin_start])
            perc_under_f_bin[-1] = (filter_lhs[j+1] - spec_lhs[bin_end])/(spec_lhs[bin_end+1] - spec_lhs[bin_end])

        spec_flux_d = spec[bin_start:bin_end+1, 1]
        widths = spec_widths[bin_start:bin_end+1]
        rebinned.append(np.sum(perc_under_f_bin*spec_flux_d)/np.sum(perc_under_f_bin))
    return np.array(rebinned)

normfilt = np.zeros(4*2, dtype="float")
normfilt.shape = (4, 2)
normfilt[:,0] = 4962.5 + np.arange(4)*25.0
normfilt[:,1] = np.ones(4)

newagenorms = np.zeros(3*700, dtype="float")
newagenorms.shape = (700, 3)

newSED_file = open(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII")
newages = np.array(newSED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
newSED_file.close()

newspecdata = np.genfromtxt(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr
newagevals = np.array([69, 89, 107])

oldagenorms = np.zeros(9*700, dtype="float")
oldagenorms.shape = (700, 9)

oldSED_file = open(SEDpath+"models/Padova1994/chabrier/bc2003_hr_stelib_m62_chab_ssp.ised_ASCII")
oldages = np.array(oldSED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
oldSED_file.close()
#print ages[129], ages[133], ages[136], ages[138], ages[141], ages[144], ages[149], ages[152], ages[156]

D = np.loadtxt("IGM_Da_Db.txt")

if arg < 3:

    for j in range(arg, arg+1):
        newspectrum = np.array([newspecdata[0, :], newspecdata[newagevals[j]+1, :]]).T
        newsynmags = np.zeros(13*700, dtype="float")
        newsynmags.shape = (700, 13)
        newsynmags[:,0] = np.arange(0.01, 7.01, 0.01)

        for i in range(1, 701):
            z = 0.01*i
            if newages[j]*(10**-9) < 14.00 - cosmo.lookback_time(z).value:
                print "Calculating colours for age: " + str(newages[newagevals[j]]) + " and redshift: " + str(z) + ", total age: " + str(newages[newagevals[j]]*(10**-9) + cosmo.lookback_time(z).value )
                newzspectrum = np.copy(newspectrum)
                for k in range(len(newzspectrum)):
                    if newzspectrum[k,0] < 912.:
                        newzspectrum[k,1] = 0.
                    elif newzspectrum[k,0] > 912. and newzspectrum[k,0] < 1026.:
                        newzspectrum[k,1] = newzspectrum[k,1]*(1-D[i-1,1])*(1-D[i-1,0])
                    elif newzspectrum[k,0] > 1026. and newzspectrum[k,0] < 1216.:
                        newzspectrum[k,1] = newzspectrum[k,1]*(1-D[i-1,0])
                        
                """      
                pylab.figure()
                pylab.plot(newspectrum[:,0], newspectrum[:,1])
                pylab.plot(newzspectrum[:,0], newzspectrum[:,1])
                pylab.xlim(0, 3000)
                pylab.show()
                """
                newzspectrum[:,1] = newzspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
                newzspectrum[:,1] = newzspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z).value*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
                newzspectrum[:,1] = newzspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
                newagenorms[i-1, j] = 25.*np.sum(rebin_to_filter(newzspectrum, normfilt))
                newzspectrum[:,1] = newzspectrum[:,1]/newagenorms[i-1, j] #reduce flux by a factor of 1/(1+z) to account for redshifting
                newzspectrum[:,0] = newzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting
                newth_mags = s.AB_mags_UDS(newzspectrum) #np.ones(12)#
                newsynmags[i-1, 1:] = newth_mags
                
        np.savetxt("../models/const/synmagsUDS_age_" + str(newages[newagevals[j]]) + ".txt", newsynfluxes)
        np.savetxt("../models/const/agenormsUDS_" + str(newages[newagevals[j]]) + ".txt", newagenorms[:,j])

else:

    oldspecdata = np.genfromtxt(SEDpath+"models/Padova1994/chabrier/bc2003_hr_stelib_m62_chab_ssp.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr
    oldagevals = np.array([129,133,136,138,141,144,149,152,156])

    for j in range(arg-3, arg-2):
        oldspectrum = np.array([oldspecdata[0, :], oldspecdata[oldagevals[j]+1, :]]).T
        oldsynmags = np.zeros(13*700, dtype="float")
        oldsynmags.shape = (700, 13)
        oldsynmags[:,0] = np.arange(0.01, 7.01, 0.01)

        for i in range(1, 701):
            z = 0.01*i
            if oldages[oldagevals[j]]*(10**-9) < 14.00 - cosmo.lookback_time(z).value:
                print "Calculating colours for age: " + str(oldages[oldagevals[j]]) + " and redshift: " + str(z) + ", total age: " + str(oldages[oldagevals[j]]*(10**-9) + cosmo.lookback_time(z).value)
                oldzspectrum = np.copy(oldspectrum)
                for k in range(len(oldzspectrum)):
                    if oldzspectrum[k,0] < 912.:
                        oldzspectrum[k,1] = 0.
                    elif oldzspectrum[k,0] > 912. and oldzspectrum[k,0] < 1026.:
                        oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[int(z/0.01 - 1),1])*(1-D[int(z/0.01 - 1),0])
                    elif oldzspectrum[k,0] > 1026. and oldzspectrum[k,0] < 1216.:
                        oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[int(z/0.01 - 1),0])
                oldzspectrum[:,1] = oldzspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
                oldzspectrum[:,1] = oldzspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z).value*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
                oldzspectrum[:,1] = oldzspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
                oldagenorms[i-1, j] = 25.*np.sum(rebin_to_filter(oldzspectrum, normfilt))
                oldzspectrum[:,1] = oldzspectrum[:,1]/oldagenorms[i-1, j]
                oldzspectrum[:,0] = oldzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting
                oldth_mags = s.AB_mags_UDS(oldzspectrum)                   
                oldsynmags[i-1, 1:] = oldth_mags

        np.savetxt("../models/burst/synmagsUDS_age_" + str(oldages[oldagevals[j]]) + ".txt", oldsynfluxes)
        np.savetxt("../models/burst/agenormsUDS_" + str(oldages[oldagevals[j]]) + ".txt", oldagenorms[:,j])

