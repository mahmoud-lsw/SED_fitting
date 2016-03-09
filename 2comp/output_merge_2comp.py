import numpy as np
from subprocess import call

file1 = np.loadtxt("photoz_UDS_0_13.txt") #_UDS
file2 = np.loadtxt("photoz_UDS_14_25.txt")
file3 = np.loadtxt("photoz_UDS_26_37.txt")
file4 = np.loadtxt("photoz_UDS_38_50.txt")

output = np.zeros(11*len(file1))
output.shape = (len(file1), 11)

for m in range(len(file1)):
    minarr = np.argmin([file1[m,8], file2[m,8], file3[m,8], file4[m,8]])
    if minarr == 0:
        output[m,:9] = file1[m,:]
    elif minarr == 1:
        output[m,:9] = file2[m,:]
    elif minarr == 2:
        output[m,:9] = file3[m,:]
    elif minarr == 3:
        output[m,:9] = file4[m,:]
    else:
        print "well shit..."
        
    old_norm = np.loadtxt("../../models/burst/agenormsUDS_" + str(output[m, 3]) + "_EBV_" + str(output[m,6]) + ".txt")[int(output[m, 2]/0.01) - 1] ####UDS
    new_norm = np.loadtxt("../../models/const/agenormsUDS_" + str(output[m, 4]) + "_EBV_" + str(output[m,6]) + ".txt")[int(output[m, 2]/0.01) - 1]
    if output[m, 5] == 1.:
        output[m, 9] = 0.
        output[m, 10] = output[m, 6]/old_norm
    else:
        output[m, 9] = output[m, 7]/new_norm
        output[m, 10] = output[m, 7]*((output[m, 5]/(1.-output[m, 5]))/old_norm + output[m, 4]/new_norm)  #4 is age new, 6 is old_modifier
        

np.savetxt("photoz_2comp_UDS.txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass") #_UDS
call(["rm", "photoz_UDS_0_13.txt", "photoz_UDS_14_25.txt", "photoz_UDS_26_37.txt", "photoz_UDS_38_50.txt"])

##"obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi"
