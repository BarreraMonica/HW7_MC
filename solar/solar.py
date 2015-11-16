import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt("solar.dat")
a=data[:,0]
b=data[:,1]
c=data[:,2]
d=data[:,3]
l=data[:,4]

plt.figure()
f,fig = plt.subplots(3,3,figsize=(12,10))

fig[0,0].scatter(a,b)
fig[0,0].set_ylabel('b')
fig[0,0].set_xlabel('a')

fig[1,0].scatter(a,c)
fig[1,0].set_ylabel('c')
fig[1,0].set_xlabel('a')

fig[2,0].scatter(a,d)
fig[2,0].set_ylabel('d')
fig[2,0].set_xlabel('a')

fig[0,1].axis('off')

fig[1,1].scatter(b,c)
fig[1,1].set_ylabel('c')
fig[1,1].set_xlabel('b')

fig[2,1].scatter(b,d)
fig[2,1].set_ylabel('d')
fig[2,1].set_xlabel('b')

fig[0,2].axis('off')

fig[1,2].axis('off')

fig[2,2].scatter(c,d)
fig[2,2].set_ylabel('d')
fig[2,2].set_xlabel('c')

plt.savefig("Probabilidades_solar.pdf",format="pdf")
plt.close()
