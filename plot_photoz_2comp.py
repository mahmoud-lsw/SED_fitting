import numpy as np
import pylab
import matplotlib as mpl

#data = np.loadtxt("photoz/photoz_EBV1.5.txt")
#data = np.loadtxt("photoz/photoz_allbandobj_v3.txt")
#data = np.loadtxt("photoz/photoz_ratios_EBV2.5.txt")
#data = np.loadtxt("photoz/photoz_medianratios_EBV0.5.txt")
data = np.loadtxt("photoz/photoz_2comp_medratios_EBV1.5.txt")

# obj_no spec_z phot_z age_old age_new f_old_V old_modifier EBV norm chi

### Find the number of data points, catastrophic outliers and sigma_dz, sigma_dz_NMAD and sigma_dz_clipped
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

data_nocat = np.zeros(len(gooddata)*10)
data_nocat.shape = (len(gooddata), 10)

for i in range(len(gooddata)):
    data_nocat[i, :] = data[gooddata[i], :]

sigma_dz_clipped = np.std((data_nocat[:,1] - data_nocat[:,2])/(1+data_nocat[:,1]))

print "Including catastrophic outliers" 
print "No of data points: " + str(no) 
print "Catastrophic outliers: " + str(cat)
print "sigma_dz: " + str(sig_dz)
print "sigma_dz_NMAD: " + str(MAD*1.483)
print "sigma_dz_clipped: " + str(sigma_dz_clipped)


### Calculate star formation rates and stellar masses for the full sample and for the well fit sample
SFR = np.copy(data[:,0])
mass = np.copy(data[:,0])
good_SFR = np.copy(data_nocat[:,0])
good_mass = np.copy(data_nocat[:,0])

for i in range(len(data)):
    if data[i,5] == 1.:
        SFR[i] = 0.00001 #actually zero but choose arbitrary small value to plot on log scale.
        mass[i] = data[i,8]*data[i,6]
    else:
        SFR[i] = data[i,8]
        mass[i] = data[i,8]*(data[i,6] + data[i,4])

for j in range(len(data_nocat)):
    if data_nocat[j,5] == 1.:
        good_SFR[j] = 0.00001 #actually zero but choose arbitrary small value to plot on log scale.
        good_mass[j] = data_nocat[j,8]*data_nocat[j,6]
    else:
        good_SFR[j] = data_nocat[j,8]
        good_mass[j] = data_nocat[j,8]*(data_nocat[j,6] + data_nocat[j,4]) #4 is age new, 6 is old_modifier


### Various plotting codes

#Plots spec_z vs phot_z with a colorbar of stellar pop age
maxxy = np.max([np.max(data[:,1]), np.max(data[:,2])]) + 0.05
pylab.figure()
pylab.scatter(data[:,1], data[:,2], c = data[:,3], norm=mpl.colors.LogNorm())
pylab.plot([0., maxxy], [0., maxxy], color="black", lw=2)
pylab.plot([0., maxxy], [0. - 0.15, 0.85*maxxy -0.15], color="grey")
pylab.plot([0., maxxy], [0. + 0.15, 1.15*maxxy + 0.15], color="grey")
pylab.xlim(0, maxxy)
pylab.ylim(0, maxxy)
pylab.ylabel("phot_z")
pylab.xlabel("spec_z")
cbar = pylab.colorbar()
cbar.set_label("Age of Stellar Pop. (yrs)", size=16)
pylab.show()


#Plots age of stellar population vs extinction EBV value
pylab.figure()
pylab.scatter(data[:, 3]*10**-9, data[:,7], color="red")
pylab.scatter(data_nocat[:, 3]*10**-9, data_nocat[:,7], color="blue")
pylab.xlabel("Age of Stellar Pop. (Gyr)", size="16")
pylab.xlim(5*10**-3, 15)
pylab.xscale("log")
pylab.ylabel("E(B - V) (mag)", size="16")
pylab.show()


#Plots stellar mass vs star formation rate
pylab.figure()
pylab.scatter(mass, SFR, color="red")
pylab.scatter(good_mass, good_SFR, color="blue")
pylab.xlabel("Stellar Mass (Solar Masses)", size="16")
pylab.yscale("log")
pylab.xscale("log")
pylab.ylabel("SFR (Solar Masses per year)", size="16")
pylab.ylim(10**-6, 5*10**3)
pylab.show()





#The following code investigates the effects of changing the range of allowed reddening values.
"""
highred = []

for i in range(len(data)):
    if data[i,5] > 0.5:
        highred.append(i)

print len(highred)
diffred = np.zeros(8*len(highred), dtype="float")
print len(highred), len(highred)*8
diffred.shape = (len(highred), 8)
oldred = np.copy(diffred)

for i in range(len(highred)):
    diffred[i,:] = data[highred[i],:]
    oldred[i,:] = datalored[highred[i],:]


maxage = np.loadtxt("max_redshift_for_age.txt")
pylab.figure()
#pylab.scatter(datalored[:,3]*10**-9, datalored[:,2], color="blue")
for i in range(len(highred)):
    pylab.plot(np.array([diffred[i,3], datalored[highred[i],3]])*10**-9, np.array([diffred[i,2], datalored[highred[i],2]]), color="black", zorder=0)
    pylab.scatter(datalored[highred[i],3]*10**-9, datalored[highred[i],2], color="blue")
pylab.scatter(diffred[:,3]*10**-9, diffred[:,2], color="red")
pylab.plot(maxage[:,0]*10**-9, maxage[:,1], lw=2, color="black")
pylab.ylabel("phot_z", size=16)
pylab.xlabel("Age (Gyr)", size=16)
pylab.xscale("log")
pylab.xlim(0.009, 15)
pylab.show()


maxxy = 5.
pylab.figure()
#pylab.scatter(datalored[:,3]*10**-9, datalored[:,2], color="blue")
for i in range(len(highred)):
    pylab.plot(np.array([diffred[i,1], datalored[highred[i],1]]), np.array([diffred[i,2], datalored[highred[i],2]]), color="black", zorder=0)
    pylab.scatter(datalored[highred[i],1], datalored[highred[i],2], color="blue")
pylab.scatter(diffred[:,1], diffred[:,2], color="red")
pylab.ylabel("phot_z", size=16)
pylab.xlabel("spec_z", size=16)
pylab.plot([0, 5], [0, 5], color="black")
pylab.plot([0., maxxy], [0., maxxy], color="black", lw=2)
pylab.plot([0., maxxy], [0. - 0.15, 0.85*maxxy -0.15], color="grey")
pylab.plot([0., maxxy], [0. + 0.15, 1.15*maxxy + 0.15], color="grey")
pylab.xlim(0, maxxy)
pylab.ylim(0, maxxy)
pylab.show()


pylab.figure()
for i in range(len(highred)):
    pylab.plot(np.array([diffred[i,3]*10**-9, datalored[highred[i],3]*10**-9]), np.array([diffred[i,5], datalored[highred[i],5]]), color="black", zorder=0)
    pylab.scatter(datalored[highred[i],3]*10**-9, datalored[highred[i],5], color="blue")
pylab.scatter(diffred[:,3]*10**-9, diffred[:,5], color="red")
pylab.ylabel("EBV", size=16)
pylab.xlabel("Age (Gyr)", size=16)
pylab.xscale("log")
pylab.xlim(0.009, 15)
pylab.show()


pylab.figure()
for i in range(len(highred)):
    pylab.plot(np.array([diffred[i,2], datalored[highred[i],2]]), np.array([diffred[i,5], datalored[highred[i],5]]), color="black", zorder=0)
    pylab.scatter(datalored[highred[i],2], datalored[highred[i],5], color="blue")
pylab.scatter(diffred[:,2], diffred[:,5], color="red")
pylab.ylabel("EBV", size=16)
pylab.xlabel("phot_z", size=16)
pylab.show()


olddz = (oldred[:,1] - oldred[:,2])/(1+oldred[:,1])
old_sig_dz = np.std(olddz)
newdz = (diffred[:,1] - diffred[:,2])/(1+diffred[:,1])
new_sig_dz = np.std(newdz)

print "old: " + str(old_sig_dz) + " new: " + str(new_sig_dz)
"""
