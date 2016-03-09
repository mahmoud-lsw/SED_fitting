import numpy as np

file1 = np.loadtxt("photoz_0_467.txt")
file2 = np.loadtxt("photoz_468_935.txt")
file3 = np.loadtxt("photoz_936_1402.txt")
file4 = np.loadtxt("photoz_1403_1870.txt")

output = np.copy(file1)

for m in range(len(file1)):
    minarr = np.argmin([file1[m,9], file2[m,9], file3[m,9], file4[m,9]])
    if minarr == 0:
        output[m,:] = file1[m,:]
    elif minarr == 1:
        output[m,:] = file2[m,:]
    elif minarr == 2:
        output[m,:] = file3[m,:]
    elif minarr == 3:
        output[m,:] = file4[m,:]
    else:
        print "well shit..."

np.savetxt("photoz_2comp_zpt.txt", output, header="obj_no spec_z phot_z age_old age_new f_old_V old_modifier EBV norm chi")
