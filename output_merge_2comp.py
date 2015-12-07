import numpy as np

file1 = np.loadtxt("photoz_0_13.txt")
file2 = np.loadtxt("photoz_14_25.txt")
file3 = np.loadtxt("photoz_26_37.txt")
file4 = np.loadtxt("photoz_38_50.txt")

output = np.zeros(11*len(file1))
output.shape = (11, len(file1))

for m in range(len(file1)):
    minarr = np.argmin([file1[m,8], file2[m,8], file3[m,8], file4[m,8]])
    if minarr == 0:
        output[m,:] = file1[m,:9]
    elif minarr == 1:
        output[m,:] = file2[m,:9]
    elif minarr == 2:
        output[m,:] = file3[m,:9]
    elif minarr == 3:
        output[m,:] = file4[m,:9]
    else:
        print "well shit..."
        
    old_norm = np.loadtxt("models/burst/agenorms_" + str(output[m, 3]) + ".txt")[int(output[m, 2]/0.01) - 1]
    new_norm = np.loadtxt("models/const/agenorms_" + str(output[m, 4]) + ".txt")[int(output[m, 2]/0.01) - 1]
    
    if output[m, 4] == 1.:
        output[m, 9] = 0.
        output[m, 10] = output[m, 6]/old_norm
    else:
        output[m, 9] = output[m, 7]/new_norm
        output[m, 10] = output[m, 7]*((output[m, 5]/(1.-output[m, 5]))/old_norm + output[m, 4]/new_norm)  #4 is age new, 6 is old_modifier
        

np.savetxt("photoz_2comp.txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass")


##"obj_no spec_z phot_z age_old age_new f_old_V EBV norm chi"
