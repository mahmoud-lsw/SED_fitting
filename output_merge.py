import numpy as np

file1 = np.loadtxt("photoz_0_20.txt")
file2 = np.loadtxt("photoz_21_42.txt")
file3 = np.loadtxt("photoz_43_64.txt")
file4 = np.loadtxt("photoz_65_126.txt")

print file1, file2, file3, file4
raw_input()

output = np.copy(file1)

for m in range(len(file1)):
    minarr = np.argmin([file1[m,7], file2[m,7], file3[m,7], file4[m,7]])
    print minarr
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
