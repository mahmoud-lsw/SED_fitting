import numpy as np
import pylab

# obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass
data = np.loadtxt("photoz_feb16/photoz_3tau_42age.txt")

no = len(data)
cat = 0.
dz = (data[:,1] - data[:,2])/(1+data[:,1])
sig_dz = np.std(dz)
MAD = np.median(np.abs(dz - np.median(dz)))

gooddata = []

for i in range(len(data)):
    if np.abs(dz[i]) > 0.15:
        cat = cat+1.
    else:
        gooddata.append(i)

data_nocat = np.zeros(len(gooddata)*8)
data_nocat.shape = (len(gooddata), 8)

for i in range(len(gooddata)):
    data_nocat[i, :] = data[gooddata[i], :]

sigma_dz_clipped = np.std((data_nocat[:,1] - data_nocat[:,2])/(1+data_nocat[:,1]))

print "Including catastrophic outliers" 
print "No of data points: " + str(no) 
print "Catastrophic outliers: " + str(cat)
print "sigma_dz: " + str(sig_dz)
print "sigma_dz_NMAD: " + str(MAD*1.483)
print "sigma_dz_clipped: " + str(sigma_dz_clipped)


#Plots spec_z vs phot_z with a colorbar of stellar pop age
maxxy = np.max([np.max(data[:,1]), np.max(data[:,2])]) + 0.05
pylab.figure()
pylab.scatter(data[:,1], data[:,2])#, c = data[:,5])#, norm=mpl.colors.LogNorm()) #3
pylab.plot([0., maxxy], [0., maxxy], color="black", lw=2)
pylab.plot([0., maxxy], [0. - 0.15, 0.85*maxxy -0.15], color="grey")
pylab.plot([0., maxxy], [0. + 0.15, 1.15*maxxy + 0.15], color="grey")
pylab.xlim(0, maxxy)
pylab.ylim(0, maxxy)
pylab.ylabel("phot_z_ACC")
pylab.xlabel("spec_z")
#cbar = pylab.colorbar()
#cbar.set_label("Age of Stellar Pop. (yrs)", size=16)
pylab.show()



results = np.loadtxt("../uds_validation_results_all.cat")

pylab.figure()
pylab.scatter(data[:,2], results[:,2])
pylab.plot([0., 7.], [0., 7.], color="black", lw=2)
pylab.xlim(0, 7)
pylab.ylim(0, 7)
pylab.xlabel("phot_z_ACC")
pylab.ylabel("phot_z_Alice_bc03")
pylab.show()

names = ["z_alice_bc03",   "z_alice_cos",   "z_barros",   "z_ciras",   "z_emq_peg",   "z_emq_pca",   "z_fink",   "z_font",   "z_micol",   "z_shegy",   "z_will_eazy",   "z_will_zebra",   "z_rjm",   "z_jarvis"]   

for i in range(14):

    no = len(results)
    cat = 0.
    dz = (results[:,1] - results[:,2+i])/(1+results[:,1])
    sig_dz = np.std(dz)
    MAD = np.median(np.abs(dz - np.median(dz)))

    goodresults = []

    for k in range(len(results)):
        if np.abs(dz[k]) > 0.15:
            cat = cat+1.
        else:
            goodresults.append(k)

    results_nocat = np.zeros(len(goodresults)*16)
    results_nocat.shape = (len(goodresults), 16)

    for k in range(len(goodresults)):
        results_nocat[k, :] = results[goodresults[k], :]

    sigma_dz_clipped = np.std((results_nocat[:,1] - results_nocat[:,2+i])/(1+results_nocat[:,1]))

    print " "
    print names[i]
    print "Including catastrophic outliers" 
    print "No of data points: " + str(no) 
    print "Catastrophic outliers: " + str(cat)
    print "sigma_dz: " + str(sig_dz)
    print "sigma_dz_NMAD: " + str(MAD*1.483)
    print "sigma_dz_clipped: " + str(sigma_dz_clipped)



    #Plots spec_z vs phot_z with a colorbar of stellar pop age
    maxxy = np.max([np.max(results[:,1]), np.max(results[:,2+i])]) + 0.05
    pylab.figure()
    pylab.scatter(results[:,1], results[:,2+i])#, c = results[:,5])#, norm=mpl.colors.LogNorm()) #3
    pylab.plot([0., maxxy], [0., maxxy], color="black", lw=2)
    pylab.plot([0., maxxy], [0. - 0.15, 0.85*maxxy -0.15], color="grey")
    pylab.plot([0., maxxy], [0. + 0.15, 1.15*maxxy + 0.15], color="grey")
    pylab.xlim(0, maxxy)
    pylab.ylim(0, maxxy)
    pylab.ylabel(names[i])
    pylab.xlabel("spec_z")
    #cbar = pylab.colorbar()
    #cbar.set_label("Age of Stellar Pop. (yrs)", size=16)
    pylab.show()
