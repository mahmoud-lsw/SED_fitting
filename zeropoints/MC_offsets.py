import numpy as np
import pylab
from scipy.optimize import minimize
import random as r

def chisq(param, binmidp, weights):
    errs = np.sqrt(weights)
    for i in range(len(errs)):
        if errs[i] == 0.:
            errs[i] = 1.
    return np.sum(((param[0]*(1/param[2]/np.sqrt(2*np.pi))*np.exp(-((binmidp - param[1])**2)/2./param[2]**2) - weights)/errs)**2)

data = np.loadtxt("ratios/flux_ratios_offsetcalc.txt")

for j in range(12):
    weights, bins = np.histogram(data[:,j], bins=40, range=(0.8, 1.2))#MC
    weights = np.array(weights, dtype="float")
    binmidp = np.copy(bins[0:-1])
    for i in range(len(binmidp)):
        binmidp[i] = np.mean(bins[i:i+2])
               
    optresult = minimize(chisq, [5., 1., 0.01], args=(binmidp, weights), method="Nelder-Mead", options={"maxfev":9999, "maxiter":9999})
    param = optresult["x"]
    #print param, optresult["success"], optresult["fun"]
    model = param[0]*(1/param[2]/np.sqrt(2*np.pi))*np.exp(-((binmidp - param[1])**2)/2./param[2]**2)
        
    pylab.figure()
    pylab.plot(binmidp, weights, color="blue")
    pylab.plot(binmidp, model, color="red")
    pylab.plot([1., 1.], [0., np.max(weights)*1.2], color="black", lw=2)
    pylab.ylim(0, np.max(weights)*1.2)
    pylab.show()


meanvals = np.zeros(12000, dtype="float")
meanvals.shape = (1000, 12)
"""
for k in range(1000):
    print k
    dataMC = np.copy(data[0:1000, :])
    for i in range(1000):
        dataMC[i,:] = data[int(r.random()*len(data)), :]

    for j in range(12):
        weights, bins = np.histogram(dataMC[:,j], bins=40, range=(0.8, 1.2))#MC
        weights = np.array(weights, dtype="float")
        binmidp = np.copy(bins[0:-1])

        for i in range(len(binmidp)):
            binmidp[i] = np.mean(bins[i:i+2])
               
        optresult = minimize(chisq, [5., 1., 0.01], args=(binmidp, weights), method="Nelder-Mead", options={"maxfev":9999, "maxiter":9999})
        param = optresult["x"]
        #print param, optresult["success"], optresult["fun"]
        #model = param[0]*(1/param[2]/np.sqrt(2*np.pi))*np.exp(-((binmidp - param[1])**2)/2./param[2]**2)
        
        pylab.figure()
        pylab.plot(binmidp, weights, color="blue")
        pylab.plot(binmidp, model, color="red")
        pylab.plot([1., 1.], [0., np.max(weights)*1.2], color="black", lw=2)
        pylab.ylim(0, np.max(weights)*1.2)
        pylab.show()
        
        meanvals[k,j] = param[1]

np.savetxt("meanvals.txt", meanvals)
"""
meanvals = np.loadtxt("ratios/meanvals.txt")

offsets = np.zeros(24, dtype="float")
offsets.shape = (12, 2)

for i in range(12):
    
    pylab.figure()
    pylab.hist(meanvals[:, i])
    pylab.show()
    
    offsets[i,:] = np.mean(meanvals[:, i]), np.std(meanvals[:, i])
#print offsets
"""
for i in range(12):
    if (offsets[i,0] - 1) < offsets[i,1]:
        print 1., 0.
    else:
        print offsets[i,0], offsets[i,1]
"""
for i in range(12): #use medians instead of all that complicated stuff
    offsets[i,0] = np.median(data[i,:])
    offsets[i,1] = 0.
print offsets
np.savetxt("mean_ratios_with_errors.txt", offsets)
