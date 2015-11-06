import numpy as np
import pylab
from scipy.optimize import minimize
import random as r

f = open("offsets.txt", "a")

def chisq(param, binmidp, weights):
    errs = np.sqrt(weights)
    for i in range(len(errs)):
        if errs[i] == 0.:
            errs[i] = 1.
    return np.sum(((param[0]*(1/param[2]/np.sqrt(2*np.pi))*np.exp(-((binmidp - param[1])**2)/2./param[2]**2) - weights)/errs)**2)

data = np.loadtxt("flux_diffs_offsetcalc.txt")


"""
for k in range(1000):

    print k
    dataMC = np.copy(data[0:1000, :])
    for i in range(1000):
        dataMC[i,:] = data[int(r.random()*len(data)), :]

    for j in range(12):
        weights, bins = np.histogram(dataMC[:,j], bins=200, range=(-1, 1))
        weights = np.array(weights, dtype="float")
        binmidp = np.copy(bins[0:-1])

        for i in range(len(binmidp)):
            binmidp[i] = np.mean(bins[i:i+2])
               
        optresult = minimize(chisq, [5., 0., 0.1], args=(binmidp, weights), method="Nelder-Mead", options={"maxfev":9999, "maxiter":9999})
        param = optresult["x"]
        #print param, optresult["success"], optresult["fun"]
        model = param[0]*(1/param[2]/np.sqrt(2*np.pi))*np.exp(-((binmidp - param[1])**2)/2./param[2]**2)
        
        pylab.figure()
        pylab.plot(binmidp, weights, color="blue")
        pylab.plot(binmidp, model, color="red")
        pylab.plot([0., 0.], [0., np.max(weights)*1.2], color="black", lw=2)
        pylab.ylim(0, np.max(weights)*1.2)
        #pylab.show()
        
        f.write(str(param[1]) + " ")

    f.write(" \n")
f.close()
"""
meanvals = np.loadtxt("offsets.txt")

for i in range(12):
    """
    pylab.figure()
    pylab.hist(meanvals[:, i])
    pylab.show()
    """
    print np.mean(meanvals[:, i]), np.std(meanvals[:, i])
