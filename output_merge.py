import numpy as np

file1 = np.loadtxt("photoz_0_20.txt")
file2 = np.loadtxt("photoz_21_42.txt")
file3 = np.loadtxt("photoz_43_64.txt")
file4 = np.loadtxt("photoz_65_126.txt")

output = np.copy(file1)

for m in range(len(file1)):
    minarr = np.argmin([file1[:,7], file2[:,7], file3[:,7], file4[:,7]])
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

np.savetxt("photoz.txt", output)
