import numpy as np
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
import pylab
import sys
SEDpath = "../progs/bc03/"

tau = int(sys.argv[1])

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
    
zarr = np.arange(0.01, 5.001, 0.01)
lbtarr = cosmo.lookback_time(zarr).value


### Load up the data necessary to build the "filter" to rebin the model spectra onto

hdulist1 = fits.open("../VANDELS_data/spectra/sc_206806_UDS_P1M1_MR_Q1_029_1.fits")
wavzpt = hdulist1[0].header["CRVAL1"]
dwav = hdulist1[0].header["CDELT1"]
fluxes1 = hdulist1[4].data#*10**19

maxwav = wavzpt + dwav*(len(fluxes1))
objwavs_nobin = np.arange(wavzpt, maxwav, dwav)
tauvals = np.array([0.05, 1, 10], dtype="str")
objwavs = specbin(objwavs_nobin, 2)

normfilt = np.zeros(4*2, dtype="float")
normfilt.shape = (4, 2)
normfilt[:,0] = 4962.5 + np.arange(4)*25.0
normfilt[:,1] = np.ones(4)

objfilter = np.ones(len(objwavs[98:-254])*2, dtype="float")
objfilter.shape = (len(objwavs[98:-254]), 2)
objfilter[:,0] = objwavs[98:-254]

output = np.zeros(len(objfilter)*500, dtype="float")
output.shape = (len(objfilter), 500)

agenorms = np.zeros(128*500, dtype="float")
agenorms.shape = (128, 500)

### Load up the spectral models to fit

for t in range(tau, tau+1):
    SED_file = open(SEDpath+"P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs/out.ised_ASCII")
    ages = np.array(SED_file.readline().split(), dtype="float")[1:] #ages[0] = 0.
    SED_file.close()
    np.savetxt("ages.txt", ages[69:197])
    specdata = np.genfromtxt(SEDpath+"P94_Chab_zsun_t" + tauvals[t] + "Gyr_noabs/out.ised_ASCII", skip_header = 6, skip_footer=12, usecols=np.arange(1, 6918), dtype="float") #specdata[1, :] is for 0 Gyr

    for j in range(69, 197):
        print "Calculating colours for age: " + str(ages[j]) + " and Tau: " + tauvals[t]
        spectrum = np.array([specdata[0, :], specdata[j+1, :]]).T

        synmags = np.zeros(13*500, dtype="float")
        synmags.shape = (500, 13)
        synmags[:,0] = np.arange(0.01, 5.01, 0.01)

        for i in range(1, 501):
            z = 0.01*i
            if ages[j]*(10**-9) < 14.00 - cosmo.lookback_time(z).value:

                zspectrum = np.copy(spectrum)
                zspectrum[:,1] = zspectrum[:,1]*3.826*10**33 #luminosity in erg/s/A
                zspectrum[:,1] = zspectrum[:,1]/(4*np.pi*(cosmo.luminosity_distance(z).value*3.086*10**24)**2) #convert to observed flux at given redshift in erg/s/A/cm^2
                zspectrum[:,1] = zspectrum[:,1]/(1+z) #reduce flux by a factor of 1/(1+z) to account for redshifting
                
                agenorms[j, i-1] = 25.*np.sum(rebin_to_filter(zspectrum, normfilt))
                zspectrum[:,1] = zspectrum[:,1]/agenorms[j,i-1] 
                zspectrum[:,0] = zspectrum[:,0]*(1+z) #change wavelength values to account for redshifting

                output[:,j] = rebin_to_filter(zspectrum, objfilter)
            else:
                output[:,j] = np.zeros(len(objfilter))

        np.savetxt("models/spec/" + tauvals[t] + "/age_" + str(ages[j]) + ".txt", output)
        np.savetxt("models/spec/" + tauvals[t] + "/agenorms_" + str(ages[j]) + ".txt", agenorms[j,:])



