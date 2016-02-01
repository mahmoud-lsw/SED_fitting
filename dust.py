import numpy as np
import pylab



wavs = np.arange(0.1, 5.01, 0.01)

coef1 = np.copy(wavs)
coef2 = np.copy(wavs)
for i in range(len(wavs)):
    if wavs[i] < 0.63:
        coef1[i] = 2.659*(-2.156 + (1.509/wavs[i]) - (0.198/(wavs[i]**2)) + (0.011/(wavs[i]**3))) + 4.05     
        coef2[i] = ((0.013/(wavs[i]**3)) - (0.232/(wavs[i]**2)) + (1.766/wavs[i]) - 0.743)
    else:
        coef1[i] = 2.659*(-1.857 + (1.040/wavs[i])) + 4.05
        coef2[i] = 1.217/wavs[i] - 0.393
            
            
#coef2 = coef2*coef1[-1]/coef2[-1]
pylab.figure()
pylab.plot(wavs, coef1, color="blue")
pylab.plot(wavs, coef2, color="red")
pylab.plot(wavs, np.zeros(len(wavs)), color="black")
pylab.show()
