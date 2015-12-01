import numpy as np

file1 = np.loadtxt("photoz_V_0_13.txt")
file2 = np.loadtxt("photoz_V_14_25.txt")
file3 = np.loadtxt("photoz_V_26_37.txt")
file4 = np.loadtxt("photoz_V_38_50.txt")

output = np.zeros(10) #len(file1)*
output.shape = (10) #len(file1), 

"""
for m in range(len(file1)):
    minarr = np.argmin([file1[m,8], file2[m,8], file3[m,8], file4[m,8]])
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
"""

for m in range(1):#len(file1)):
    
    minarr = np.argmin([file1[7], file2[7], file3[7], file4[7]])
    if minarr == 0:
        output[:8] = file1[:]
    elif minarr == 1:
        output[:8] = file2[:]
    elif minarr == 2:
        output[:8] = file3[:]
    elif minarr == 3:
        output[:8] = file4[:]
    else:
        print "well shit..."
    
    #output = np.loadtxt("photoz_V_obj1_full.txt")
    old_norm = np.loadtxt("models/spec/oldburst/agenorms_" + str(output[2]) + ".txt")[int(output[1]/0.01) - 1]
    new_norm = np.loadtxt("models/spec/newconst/agenorms_" + str(output[3]) + ".txt")[int(output[1]/0.01) - 1]
    
    if output[4] == 1.:
        output[8] = 0.
        output[9] = output[6]/old_norm #replace "data" with "output"
    else:
        output[8] = output[6]/new_norm
        output[9] = output[6]*((output[4]/(1.-output[4]))/old_norm + output[3]/new_norm)  #4 is age new, 6 is old_modifier
        
 
np.savetxt("photoz_V_2comp.txt", output, header="obj_no phot_z age_old age_new f_old_V EBV norm chi SFR stellar_mass")
