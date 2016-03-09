import numpy as np

file1 = np.loadtxt("flux_fracoffsets_2comp_UDS_0_255.txt")
file2 = np.loadtxt("flux_fracoffsets_2comp_UDS_256_510.txt")
file3 = np.loadtxt("flux_fracoffsets_2comp_UDS_511_765.txt")
file4 = np.loadtxt("flux_fracoffsets_2comp_UDS_766_1019.txt")

print file1.shape, file2.shape

output = np.zeros(12*1020, dtype="float")
output.shape = (1020, 12)

output[0:256,:] = file1
output[256:511,:] = file2
output[511:766,:] = file3
output[766:1020,:] = file4

np.savetxt("flux_ratios_2comp_UDS.txt", output)
