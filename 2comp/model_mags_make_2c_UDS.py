import numpy as np
from astropy.cosmology import FlatLambdaCDM
import sys
sys.path.append("/disk1/adamc/my_code")
import synphot as s
SEDpath = "/disk1/adamc/progs/bc03/"
cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3)
sys_arg = int(sys.argv[1])



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



zarr = np.arange(0.01, 7.01, 0.01)
D = np.loadtxt("../IGM_Da_Db.txt")

normfilt = np.zeros(4*2, dtype="float")
normfilt.shape = (4, 2)
normfilt[:,0] = 4962.5 + np.arange(4)*25.0
normfilt[:,1] = np.ones(4)


newSED_file = open(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII")
newages = np.array(newSED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
newagevals = np.array([69, 89, 107, 115, 123, 129, 132, 135])
newSED_file.close()
newspecdata = np.genfromtxt(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr
newsynmags_raw = np.zeros(13*700, dtype="float")
newsynmags_raw.shape = (700, 13)
newsynmags_raw[:,0] = zarr
newagenorms_raw = np.zeros(700)
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
oldsynmags_raw = np.zeros(13*700, dtype="float")
oldsynmags_raw.shape = (700, 13)
oldsynmags_raw[:,0] = zarr
oldagenorms_raw = np.zeros(700)
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


if sys_arg < 8:
    for e in range(41):
        EBV = 0.025*e
        for j in range(sys_arg, sys_arg+1):
            newspectrum = np.array([newspecdata[0, :], newspecdata[newagevals[j]+1, :]]).T
            newsynmags = np.copy(newsynmags_raw)
            newagenorms = np.copy(newagenorms_raw)
            for i in range(0, 700):
                z = zarr[i]
                if newages[newagevals[j]]*(10**-9) < 14.00 - cosmo.lookback_time(z).value:
                    print "Calculating colours for EBV: " + str(EBV) + ", age: " + str(newages[newagevals[j]]) + " and redshift: " + str(z)
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
                    newagenorms[i] = np.mean(rebin_to_filter(newzspectrum, normfilt))
                    newzspectrum[:,1] = newzspectrum[:,1]/newagenorms[i]
                    newzspectrum[:,0] = newzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting
                    newth_mags = s.AB_mags_UDS(newzspectrum)
                    
                    newsynmags[i, 1:] = newth_mags
            
            np.savetxt("../../models/const/synmagsUDS_age_" + str(newages[newagevals[j]]) + "_EBV_" + str(EBV) + ".txt", newsynmags)
            np.savetxt("../../models/const/agenormsUDS_" + str(newages[newagevals[j]]) + "_EBV_" + str(EBV) + ".txt", newagenorms)
        
        
else:
    for e in range(41):
        EBV = 0.025*e
        for j in range(sys_arg-8, sys_arg-7):
            oldspectrum = np.array([oldspecdata[0, :], oldspecdata[oldagevals[j]+1, :]]).T
            oldsynmags = np.copy(oldsynmags_raw)
            oldagenorms = np.copy(oldagenorms_raw)
            for i in range(0, 700):
                z = zarr[i]
                if oldages[oldagevals[j]]*(10**-9) < 14.00 - cosmo.lookback_time(z).value:
                    print "Calculating colours for EBV: " + str(EBV) + ", age: " + str(oldages[oldagevals[j]]) + " and redshift: " + str(z)
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
                    oldagenorms[i] = np.mean(rebin_to_filter(oldzspectrum, normfilt))
                    oldzspectrum[:,1] = oldzspectrum[:,1]/oldagenorms[i]
                    oldzspectrum[:,0] = oldzspectrum[:,0]*(1+z) #change wavelength values to account for redshifting
                    oldth_mags = s.AB_mags_UDS(oldzspectrum)                
                    oldsynmags[i, 1:] = oldth_mags

            np.savetxt("../../models/burst/synmagsUDS_age_" + str(oldages[oldagevals[j]]) + "_EBV_" + str(EBV) + ".txt", oldsynmags)
            np.savetxt("../../models/burst/agenormsUDS_" + str(oldages[oldagevals[j]]) + "_EBV_" + str(EBV) + ".txt", oldagenorms)
        
