import numpy as np

file1 = np.loadtxt("photoz_tau_0.03_to_0.09.txt")
file2 = np.loadtxt("photoz_tau_0.16_to_0.53.txt")
file3 = np.loadtxt("photoz_tau_0.95_to_3.08.txt")
file4 = np.loadtxt("photoz_tau_5.55_to_18.01.txt")


output = np.copy(file1)

for m in range(len(file1)):
    minarr = np.argmin([file1[m,7], file2[m,7], file3[m,7], file4[m,7]])
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

np.savetxt("photoz.txt", output, header="obj_no spec_z phot_z age tau EBV norm chi")
