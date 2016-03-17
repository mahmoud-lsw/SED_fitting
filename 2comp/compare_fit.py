import numpy as np
import pylab

# obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass

model8 = np.loadtxt("photoz/photoz_2comp_UDS_8x8.txt")
model1 = np.loadtxt("photoz/photoz_2comp_UDS_1x5.txt")

spec_z = model8[:,1]

phot_z_1 = model1[:,2]
phot_z_8 = model8[:,2]

badobj = []

for i in range(len(spec_z)):
    if np.abs(spec_z[i] - phot_z_1[i]) < np.abs(spec_z[i] - phot_z_8[i]) and phot_z_8[i] - phot_z_1[i]:
        badobj.append(i)
        
badobj = np.array(badobj)

model1badobj = np.zeros(11*len(badobj))
model1badobj.shape = (len(badobj), 11)

for a in range(len(badobj)):
    model1badobj[a,:] = model1[badobj[a],:]
    

#Plots stellar mass vs star formation rate
pylab.figure()

pylab.scatter(model1[:,4], model1[:,5], color="blue")
pylab.scatter(model1badobj[:,4], model1badobj[:,5], color="red")
pylab.xlabel("age_new", size="16")
#pylab.yscale("log")
#pylab.xscale("log")
pylab.ylabel("f_old", size="16")
#pylab.ylim(10**-5, 5*10**3)
#pylab.xlim(10**5, 5*10**12)
pylab.show()

print np.mean(model1[:,5]), np.std(model1[:,5])
print np.mean(model1badobj[:,5]), np.std(model1badobj[:,5])

"""
better = 0
same = 0
worse = 0
better_diffs = []
worse_diffs = []
better_fold = []
worse_fold = []
same_fold = []
spec_z_worse = []
phot_z_worse1 = []
spec_z_better = []

for i in range(len(spec_z)):
    #print spec_z[i], phot_z_1[i], phot_z_8[i]

    if np.abs(spec_z[i] - phot_z_1[i]) > np.abs(spec_z[i] - phot_z_8[i]):
        better = better + 1
        better_diffs.append(phot_z_1[i] - phot_z_8[i])
        spec_z_better.append(spec_z[i])
        better_fold.append(model1[i,5])
        #print "better"
    elif np.abs(spec_z[i] - phot_z_1[i]) == np.abs(spec_z[i] - phot_z_8[i]):    
        same = same + 1
        #print "same"
        same_fold.append(model1[i,5])
    elif spec_z[i] - phot_z_1[i] < 0.2:
        #print "worse"
        spec_z_worse.append(spec_z[i])
        worse_diffs.append(phot_z_8[i] - phot_z_1[i])
        phot_z_worse1.append(phot_z_1[i])
        worse = worse + 1
        worse_fold.append(model1[i,5])
    #raw_input()
   
print worse
pylab.figure()
pylab.scatter(np.array(spec_z_worse), np.array(worse_diffs), c = np.array(spec_z_worse) - np.array(phot_z_worse1))
pylab.plot([0, 7], [0, 0], color="grey")
pylab.plot([0, 7], [-0.05, -0.05], color="grey")
pylab.xlabel("spec_z", size=16)
pylab.ylabel("phot_z_8 - phot_z_1", size=16)

cbar = pylab.colorbar()
cbar.set_label("spec_z - phot_z_1", size=16)#Age of Stellar Pop. (yrs)

pylab.show()


down = 0

for i in range(len(spec_z)):
    if phot_z_8[i] - phot_z_1[i] < 0. and np.abs(spec_z[i] - phot_z_1[i]) > np.abs(spec_z[i] - phot_z_8[i]):
        down = down + 1
        
print down

   
    
print "Better: " + str(better)
print "Same: " + str(same)
print "Worse: " + str(worse)
print "Better diffs mean: " + str(np.mean(np.array(better_diffs)))
print "Worse diffs mean: " + str(np.mean(np.array(worse_diffs)))
print "Better mean fold: " + str(np.mean(np.array(better_fold))) + " +/- " + str(np.std(np.array(better_fold)))
print "Same mean fold: " + str(np.mean(np.array(same_fold))) + " +/- " + str(np.std(np.array(same_fold)))
print "Worse mean fold: " + str(np.mean(np.array(worse_fold))) + " +/- " + str(np.std(np.array(worse_fold)))

print np.mean(np.array(worse_diffs)), np.std(np.array(worse_diffs))
print np.mean(np.array(better_diffs)), np.std(np.array(better_diffs))

pylab.figure()
pylab.scatter(spec_z, phot_z_1, color="blue")
pylab.scatter(spec_z, phot_z_8, color="red")
pylab.xlabel("spec_z", size=16)
pylab.ylabel("phot_z_1", size=16)
pylab.plot([0, 7], [0, 7], color="grey")
pylab.show()

pylab.figure()
pylab.scatter(spec_z - phot_z_1, phot_z_1 - phot_z_8)
pylab.xlabel("spec_z - phot_z_1", size=16)
pylab.ylabel("phot_z_1 - phot_z_8", size=16)
pylab.plot([1, -1], [-1, 1], color="grey")
pylab.plot([0, 0], [-1, 1], color="grey")
pylab.plot([-1, 1], [0, 0], color="grey")
pylab.show()

spec_z_rest = []
phot_z_rest1 = []
phot_z_rest8 = []

for i in range(len(spec_z)):
    if spec_z[i] > 0.5 and spec_z[i] < 1.5 and np.abs(spec_z[i] - phot_z_1[i]) < 0.2:
        spec_z_rest.append(spec_z[i])
        phot_z_rest1.append(phot_z_1[i])
        phot_z_rest8.append(phot_z_8[i])

#print np.mean(spec_z - phot_z_1), np.std(spec_z - phot_z_1)
#print np.mean(np.array(spec_z_rest) - np.array(phot_z_rest1)), np.std(np.array(spec_z_rest) - np.array(phot_z_rest1))

"""
