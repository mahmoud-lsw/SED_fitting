import numpy as np
from subprocess import call

file1 = np.loadtxt("photoz_V_0_13.txt")
file2 = np.loadtxt("photoz_V_14_25.txt")
file3 = np.loadtxt("photoz_V_26_37.txt")
file4 = np.loadtxt("photoz_V_38_50.txt")
"""

file1 = np.loadtxt("photoz_V_0_25.txt")
file2 = np.loadtxt("photoz_V_26_50.txt")
file3 = np.loadtxt("photoz_V_51_75.txt")
file4 = np.loadtxt("photoz_V_76_100.txt")
"""
output = np.zeros(10*len(file1))
output.shape = (len(file1), 10)

for m in range(len(file1)):
    
    minarr = np.argmin([file1[m,7], file2[m,7], file3[m,7], file4[m,7]])
    if minarr == 0:
        output[m,:8] = file1[m,:]
    elif minarr == 1:
        output[m,:8] = file2[m,:]
    elif minarr == 2:
        output[m,:8] = file3[m,:]
    elif minarr == 3:
        output[m,:8] = file4[m,:]
    else:
        print "well shit..."
    
    #output = np.loadtxt("photoz_V_obj1_full.txt")
    old_norm = np.loadtxt("../../models/spec/burst/agenorms_" + str(output[m,2]) + "_EBV_" + str(output[m,5]) + ".txt")[int(output[m,1]/0.01) - 1]
    new_norm = np.loadtxt("../../models/spec/const/agenorms_" + str(output[m,3]) + "_EBV_" + str(output[m,5]) + ".txt")[int(output[m,1]/0.01) - 1]
    
    if output[m,4] == 1.:
        output[m,8] = 0.
        output[m,9] = output[m,6]/old_norm #replace "data" with "output"
    else:
        output[m,8] = output[m,6]/new_norm
        output[m,9] = output[m,6]*((output[m,4]/(1.-output[m,4]))/old_norm + output[m,3]/new_norm)
        
np.savetxt("photoz_V_2comp.txt", output, header="obj_no phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass")

call(["rm", "photoz_V_0_13.txt", "photoz_V_14_25.txt", "photoz_V_26_37.txt", "photoz_V_38_50.txt"])
