import numpy as np
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
import pylab
import sys
SEDpath = "/disk1/adamc/progs/bc03/"
arg = int(sys.argv[1])
sys.path.append("/disk1/adamc/my_code")
import synphot as s

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



"""calculates the total flux observed when a given spectrum (spec) is observed through a filter (filter).  \\
The spectrum input should contain a column of wavelength values and a column of associated flux per unit wavelength values. \\
The filter input should contain a column of wavelength values and a column of associated transmission coefficients. \\
Make sure the spectrum has wavelength coverage over the whole range of any filters used and leave a margin either side."""
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



def band_flux_dev(spec, filter):
    spec_bounds, spec_widths = make_bins(spec[:,0], make_rhs="True")
    filter_bounds, filter_widths = make_bins(filter[:,0], make_rhs="True")
    top = np.zeros(len(spec_widths)*len(filter_widths))
    bottom = np.copy(top)

    top = (np.expand_dims(filter_bounds[1:], axis=1) - np.expand_dims(spec_bounds[:-1], axis=0))/np.expand_dims(filter_widths, axis=1)
    bottom = (np.expand_dims(spec_bounds[1:], axis=0) - np.expand_dims(filter_bounds[:-1], axis=1))/np.expand_dims(filter_widths, axis=1)
    part1 = np.logical_and(top>=1, bottom>=1).astype("float")
    part2 = np.logical_and(np.logical_and(0<top, top<1), bottom>=1).astype("float")
    part3 = np.logical_and(np.logical_and(0<bottom, bottom<1), top>=1).astype("float")
    part4 = np.logical_and(np.logical_and(0<top, top<1), np.logical_and(0<bottom, bottom<1)).astype("float")
    crossover = part1 + top*part2 + bottom*part3 + (top+bottom-1)*part4 
    flux_total = np.sum(np.sum(crossover*np.expand_dims(filter[:,1], axis=1)*np.expand_dims(spec[:,1], axis=0)*np.expand_dims(filter_widths, axis=1), axis=1), axis=0)
    return flux_total
    
zarr = np.arange(0.01, 7.001, 0.01)
lbtarr = cosmo.lookback_time(zarr).value
D = np.loadtxt("../IGM_Da_Db.txt")

### Load up the data necessary to build the "filter" to rebin the model spectra onto

hdulist1 = fits.open("../../VANDELS_data/spectra/sc_206806_UDS_P1M1_MR_Q1_029_1.fits")
wavzpt = hdulist1[0].header["CRVAL1"]
dwav = hdulist1[0].header["CDELT1"]
fluxes1 = hdulist1[4].data

maxwav = wavzpt + dwav*(len(fluxes1))
objwavs_nobin = np.arange(wavzpt, maxwav, dwav)

objwavs = specbin(objwavs_nobin, 2)

objfilter = np.ones(len(objwavs[98:-254])*2, dtype="float")
objfilter.shape = (len(objwavs[98:-254]), 2)
objfilter[:,0] = objwavs[98:-254]

### Load up the spectral models to fit

old_SED_file = open(SEDpath+"models/Padova1994/chabrier/bc2003_hr_stelib_m62_chab_ssp.ised_ASCII")
oldages = np.array(old_SED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
old_SED_file.close()
oldageind = np.array([129,132,135,138,144,149,152,156])
oldspecdata = np.genfromtxt(SEDpath+"models/Padova1994/chabrier/bc2003_hr_stelib_m62_chab_ssp.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #oldspecdata[1, :] is for 0 Gyr

oldwavspec = oldspecdata[0,:]

oldwavs = oldwavspec*10**-4
oldcoef = np.copy(oldwavs)
for j in range(6917):
    if oldwavs[j] < 0.12:
        oldcoef[j] = 2.659*(-2.156 + (1.509/oldwavs[j]) - (0.198/(oldwavs[j]**2)) + (0.011/(oldwavs[j]**3))) + 4.05 - 69.44*(0.12-oldwavs[j])
    elif oldwavs[j] < 0.63:
        oldcoef[j] = 2.659*(-2.156 + (1.509/oldwavs[j]) - (0.198/(oldwavs[j]**2)) + (0.011/(oldwavs[j]**3))) + 4.05
    else:
        oldcoef[j] = 2.659*(-1.857 + (1.040/oldwavs[j])) + 4.05


oldoutput = np.zeros(len(objfilter)*700, dtype="float")
oldoutput.shape = (len(objfilter), 700)


new_SED_file = open(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII")
newages = np.array(new_SED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
new_SED_file.close()
newageind = np.array([69, 89, 107, 115, 123, 129, 132, 135])
newspecdata = np.genfromtxt(SEDpath+"P94_chab_zsun_const/P94_chab_zsun_const.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #newspecdata[1, :] is for 0 Gyr
newwavspec = newspecdata[0,:]

newwavs = newwavspec*10**-4
newcoef = np.copy(newwavs)
for j in range(6917):
    if newwavs[j] < 0.12:
        newcoef[j] = 2.659*(-2.156 + (1.509/newwavs[j]) - (0.198/(newwavs[j]**2)) + (0.011/(newwavs[j]**3))) + 4.05 - 69.44*(0.12-newwavs[j])
    elif newwavs[j] < 0.63:
        newcoef[j] = 2.659*(-2.156 + (1.509/newwavs[j]) - (0.198/(newwavs[j]**2)) + (0.011/(newwavs[j]**3))) + 4.05
    else:
        newcoef[j] = 2.659*(-1.857 + (1.040/newwavs[j])) + 4.05


newoutput = np.zeros(len(objfilter)*700, dtype="float")
newoutput.shape = (len(objfilter), 700)

normfilt = np.zeros(4*2, dtype="float")
normfilt.shape = (4, 2)
normfilt[:,0] = 4962.5 + np.arange(4)*25.0
normfilt[:,1] = np.ones(4)

oldagenorms = np.zeros(8*700, dtype="float")
oldagenorms.shape = (8, 700)

newagenorms = np.zeros(8*700, dtype="float")
newagenorms.shape = (8, 700)

if arg < 8:
    for e in range(41):
        EBV = 0.025*e
        print "new age: " + str(newages[newageind[arg]]) + ", EBV: " + str(EBV)
        for i in range(arg, arg+1):
            for j in range(len(zarr)):
                newzspectrum = np.zeros(len(newwavspec)*2, dtype="float")
                newzspectrum.shape = (len(newwavspec), 2)
                newzspectrum[:,0] = newwavspec
                newzspectrum[:,1] = newspecdata[newageind[i]+1,:]*3.826*10**33 #luminosity in erg/s/A
                newzspectrum[:,1] = newzspectrum[:,1]*10**(-EBV*newcoef/2.5)
                for k in range(len(newzspectrum)):
                    if newzspectrum[k,0] < 912.:
                        newzspectrum[k,1] = 0.
                    elif newzspectrum[k,0] > 912. and newzspectrum[k,0] < 1026.:
                        newzspectrum[k,1] = newzspectrum[k,1]*(1-D[j,1])
                    elif newzspectrum[k,0] > 1026. and newzspectrum[k,0] < 1216.:
                        newzspectrum[k,1] = newzspectrum[k,1]*(1-D[j,0])
                newzspectrum[:,1] = newzspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(zarr[j])*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
                newzspectrum[:,1] = newzspectrum[:,1]/(1+zarr[j]) #reduce flux by a factor of 1/(1+z) to account for redshifting
                newagenorms[i, j] = 25.*np.sum(rebin_to_filter(newzspectrum, normfilt))
                newzspectrum[:,1] = newzspectrum[:,1]/newagenorms[i, j] #reduce flux by a factor of 1/(1+z) to account for redshifting
                newzspectrum[:,0] = newzspectrum[:,0]*(1+zarr[j]) #change wavelength values to account for redshifting
                newoutput[:,j] = rebin_to_filter(newzspectrum, objfilter)
            
            np.savetxt("../../models/spec/const/age_" + str(newages[newageind[i]]) + "_EBV_" + str(EBV) + ".txt", newoutput)
            np.savetxt("../../models/spec/const/agenorms_" + str(newages[newageind[i]]) + "_EBV_" + str(EBV) + ".txt", newagenorms[i,:])
         
else:
    for e in range(41):
        EBV = 0.025*e
        print "old age: " + str(oldages[oldageind[arg-8]]) + ", EBV: " + str(EBV)
        for i in range(arg-8, arg-7):
            for j in range(len(zarr)):
                oldzspectrum = np.zeros(len(oldwavspec)*2, dtype="float")
                oldzspectrum.shape = (len(oldwavspec), 2)
                oldzspectrum[:,0] = oldwavspec
                oldzspectrum[:,1] = oldspecdata[oldageind[i]+1,:]*3.826*10**33 #luminosity in erg/s/A
                oldzspectrum[:,1] = oldzspectrum[:,1]*10**(-EBV*oldcoef/2.5)
                for k in range(len(oldzspectrum)):
                    if oldzspectrum[k,0] < 912.:
                        oldzspectrum[k,1] = 0.
                    elif oldzspectrum[k,0] > 912. and oldzspectrum[k,0] < 1026.:
                        oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[j,1])
                    elif oldzspectrum[k,0] > 1026. and oldzspectrum[k,0] < 1216.:
                        oldzspectrum[k,1] = oldzspectrum[k,1]*(1-D[j,0])

                oldzspectrum[:,1] = oldzspectrum[:,1]/(4.*np.pi*(cosmo.luminosity_distance(zarr[j])*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2        
                oldzspectrum[:,1] = oldzspectrum[:,1]/(1.+zarr[j]) #reduce flux by a factor of 1/(1+z) to account for redshifting
                oldagenorms[i, j] = 25.*np.sum(rebin_to_filter(oldzspectrum, normfilt))
                oldzspectrum[:,1] = oldzspectrum[:,1]/oldagenorms[i, j] #reduce flux by a factor of 1/(1+z) to account for redshifting
                oldzspectrum[:,0] = oldzspectrum[:,0]*(1+zarr[j]) #change wavelength values to account for redshifting
                oldoutput[:,j] = rebin_to_filter(oldzspectrum, objfilter)

            np.savetxt("../../models/spec/burst/age_" + str(oldages[oldageind[i]]) + "_EBV_" + str(EBV) + ".txt", oldoutput)
            np.savetxt("../../models/spec/burst/agenorms_" + str(oldages[oldageind[i]]) + "_EBV_" + str(EBV) + ".txt", oldagenorms[i,:])


