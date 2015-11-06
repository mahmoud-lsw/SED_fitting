import numpy as np
import pylab
import matplotlib as mpl

data = np.loadtxt("photoz.txt")

no = 0.
cat = 0.
dz = (data[:,1] - data[:,2])/(1+data[:,1])
sig_dz = np.std(dz)

for i in range(len(data)):
    if data[i,1] > 0.:
        no = no+1.

gooddata = []

for i in range(len(data)):
    if np.abs(dz[i]) > 0.15:
        cat = cat+1.
    else:
        gooddata.append(i)


print no, cat, sig_dz


wellfit = np.zeros(len(gooddata)*8)
wellfit.shape = (len(gooddata), 8)

for i in range(len(gooddata)):
    wellfit[i, :] = data[gooddata[i], :]


maxxy = np.max([np.max(data[:,1]), np.max(data[:,2])]) + 0.05

SFR = (1./(data[:,4]*10**9))*np.exp(-data[:,3]/(data[:,4]*10**9))*data[:,6]
mass = (1. - np.exp(-data[:,3]/(data[:,4]*10**9)))*data[:,6]
"""
for i in range(len(SFR)):
    if SFR[i] < 10**-5:
        print SFR[i], (1./(data[i,4]*10**9)), np.exp(-data[i,3]/(data[i,4]*10**9)), data[i,6]
"""



good_SFR = (1./(wellfit[:,4]*10**9))*np.exp(-wellfit[:,3]/(wellfit[:,4]*10**9))*wellfit[:,6]
good_mass = (1. - np.exp(-wellfit[:,3]/(wellfit[:,4]*10**9)))*wellfit[:,6]

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

pylab.figure()
pylab.scatter(data[:, 3]*10**-9, data[:,5], color="red")
pylab.scatter(wellfit[:, 3]*10**-9, wellfit[:,5], color="blue")
pylab.xlabel("Age of Stellar Pop. (Gyr)", size="16")
#pylab.yscale("log")
pylab.xlim(5*10**-3, 15)
pylab.xscale("log")
pylab.ylabel("E(B - V) (mag)", size="16")
pylab.show()

pylab.figure()
pylab.scatter(mass, SFR, color="red")
pylab.scatter(good_mass, good_SFR, color="blue")
pylab.xlabel("Stellar Mass (Solar Masses)", size="16")
pylab.yscale("log")
pylab.xscale("log")
pylab.ylabel("SFR (Solar Masses per year)", size="16")
pylab.ylim(10**-11, 5*10**3)
pylab.show()


