import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt("poblaciones.dat")
a=data[:,0]
b=data[:,1]
c=data[:,2]
d=data[:,3]
l=data[:,4]

plt.figure()
f,fig = plt.subplots(3,3,figsize=(12,10))

fig[0,0].scatter(a,b)
fig[0,0].set_ylabel(r'$\beta$')
fig[0,0].set_xlabel(r'$\alpha$')

fig[1,0].scatter(a,c)
fig[1,0].set_ylabel(r'$\gamma$')
fig[1,0].set_xlabel(r'$\alpha$')

fig[2,0].scatter(a,d)
fig[2,0].set_ylabel(r'$\delta$')
fig[2,0].set_xlabel(r'$\alpha$')

fig[0,1].axis('off')

fig[1,1].scatter(b,c)
fig[1,1].set_ylabel(r'$\gamma$')
fig[1,1].set_xlabel(r'$\beta$')

fig[2,1].scatter(b,d)
fig[2,1].set_ylabel(r'$\delta$')
fig[2,1].set_xlabel(r'$\beta$')

fig[0,2].axis('off')

fig[1,2].axis('off')

fig[2,2].scatter(c,d)
fig[2,2].set_ylabel(r'$\delta$')
fig[2,2].set_xlabel(r'$\gamma$')

plt.savefig("Probabilidades_poblaciones.png",format="png")
plt.close()
