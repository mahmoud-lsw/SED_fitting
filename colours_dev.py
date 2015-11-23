# !!! ALL MAGNITUDES OUTPUTTED BY THIS FILE ARE ON THE VEGA SYSTEM !!!

import numpy as np

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
def band_flux(spec, filter):
    flux_total = 0.
    spec_lhs, spec_widths = make_bins(spec[:,0])
    filter_lhs, filter_widths = make_bins(filter[:,0], make_rhs="True")

    i=0
    i_bins_in_f_bin = []

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
        flux_total = flux_total+np.sum(perc_under_f_bin*spec_flux_d*widths)*filter[j,1]
    return flux_total

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
    
if __name__ == "__main__":
    spec = np.array([[-3, 1.6264], [1, 3.2514], [5, 3.7654], [9, 2.1617], [15, 0.9525]])
    filter = np.ones(20, dtype="float")
    filter.shape = (10, 2)
    filter[:,0] = np.arange(10)
    filter[:,0] = filter[:,0] + 1
    filter[:,1] +=1
    print band_flux(spec, filter)
    print band_flux_dev(spec, filter)
    """
    #Calculate QSO fluxes in the three bands, then magnitudes and colours for a range of redshifts.
   
    redshift = []
    r_mag = []
    i_mag = []
    z_mag = []
    Y_mag = []
    J_mag = []
    H_mag = []
    K_mag = []
    W1_mag = []
    W2_mag = []
    W3_mag = []
   
    for i in range(181):
        redshift_value = 5.7+0.01*i
        print "Calculating colour for redshift " + str(redshift_value)
        zspec = np.loadtxt("high_z_templates/SDSS_mean_comp-z"+str(("%.2f" % redshift_value)))
        #QSO_r_flux= band_flux(zspec, r_filter)
        #QSO_i_flux= band_flux(zspec, i_filter) #these are commented out because the highest redshift templates the code is currently set up to generate data for don't have the right wavelength coverage, all you need to do to get r and i band magnitudes is to uncomment the relevant lines and reduce the maximum redshift the file iterates to.
        QSO_z_flux = band_flux(zspec, z_filter)
        QSO_Y_flux = band_flux(zspec, Y_filter)
        QSO_J_flux = band_flux(zspec, J_filter)
        QSO_H_flux = band_flux(zspec, H_filter)
        QSO_K_flux = band_flux(zspec, K_filter)
        QSO_W1_flux = band_flux(zspec, W1)
        QSO_W2_flux = band_flux(zspec, W2)
    
        redshift.append(redshift_value)
        #r_mag.append(-2.5*np.log10(QSO_r_flux/vega_r_flux))
        #i_mag.append(-2.5*np.log10(QSO_i_flux/vega_i_flux))
        z_mag.append(-2.5*np.log10(QSO_z_flux/vega_z_flux))
        Y_mag.append(-2.5*np.log10(QSO_Y_flux/vega_Y_flux))
        J_mag.append(-2.5*np.log10(QSO_J_flux/vega_J_flux))
        H_mag.append(-2.5*np.log10(QSO_H_flux/vega_H_flux))
        K_mag.append(-2.5*np.log10(QSO_K_flux/vega_K_flux))
        W1_mag.append(-2.5*np.log10(QSO_W1_flux/vega_W1_flux))
        W2_mag.append(-2.5*np.log10(QSO_W2_flux/vega_W2_flux))
        
    #Save trace to file.
    tosave = np.zeros((len(z_mag), 10), dtype="float")
    tosave[:,0] = redshift
    #tosave[:,1] = r_mag
    #tosave[:,2] = i_mag
    tosave[:,3] = z_mag
    tosave[:,4] = Y_mag
    tosave[:,5] = J_mag
    tosave[:,6] = H_mag
    tosave[:,7] = K_mag
    tosave[:,8] = W1_mag
    tosave[:,9] = W2_mag
    np.savetxt("NEW_trace_allbands.txt", tosave, header="redshift rmag_Vega imag_Vega zmag_Vega Ymag_Vega Jmag_Vega Hmag_Vega Kmag_Vega W1mag W2mag")
    """

