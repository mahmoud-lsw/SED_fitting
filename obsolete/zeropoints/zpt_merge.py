import numpy as np

file1 = np.loadtxt("flux_fracoffsets_0_467.txt")
file2 = np.loadtxt("flux_fracoffsets_468_935.txt")
file3 = np.loadtxt("flux_fracoffsets_936_1402.txt")
file4 = np.loadtxt("flux_fracoffsets_1403_1870.txt")

print file1.shape, file2.shape

output = np.zeros(12*1871, dtype="float")
output.shape = (1871, 12)
print output.shape, output[0:468,:].shape

output[0:468,:] = file1
output[468:936,:] = file2
output[936:1403.:] = file3
output[1403:1871.:] = file4

np.savetxt("flux_fracoffsets.txt", output)
