import numpy as np

# load data
time, ld1 = np.loadtxt("ld-data.txt", usecols=(0,1), unpack=True)
time *= 50.0 # because each frame is written after stride=50ps

print(ld1.shape, time.shape)
print(time[-1])

# load bias 
opes_b, uwall_b, lwall_b = np.loadtxt("../COLVAR", usecols=(2,8,9), unpack=True, max_rows=int(time[-1])+50)

opes_b = opes_b[::50]
uwall_b = uwall_b[::50]
lwall_b = lwall_b[::50]

print(opes_b.shape)
print(uwall_b.shape)
print(lwall_b.shape)

f_o = open("colvar-out.txt", "w")
f_o.write("#! FIELDS time ld1 opes.bias u_wall.bias l_wall.bias\n")
for i in range(len(time)):
    f_o.write(str(time[i])+" "+str(ld1[i])+" "+str(opes_b[i])+" "+str(uwall_b[i])+" "+str(lwall_b[i])+"\n")
f_o.close()


#! FIELDS time ld1 opes.bias opes.rct opes.zed opes.neff opes.nker opes.work u_wall.bias l_wall.bias
