import numpy as np

data = np.loadtxt("max_redshift_for_age.txt")

total = np.sum(data[:,1])
eighth = total/4.
print eighth

print len(data)
runtot = 0.
for i in range(len(data)):
    runtot += data[i,1]
    if runtot > eighth:
        print i
        runtot = 0.

print runtot
