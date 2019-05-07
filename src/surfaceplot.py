import numpy as np
import matplotlib.pyplot as plt
x=[0.0,2.0, 2.0, 2.0, 2.0, 4.0, 4.0, 4.0, 4.0, 7.5, 7.5, 7.5, 7.5]
y=[0.0,2.0, 4.0, 8.0, 16.0, 2.0, 4.0, 8.0, 16.0, 2.0, 4.0, 8.0, 16.0]
z=[0.0,0.228346, 0.330935, 0.503597, 0.597122, 0.255639, 0.478261, 0.635714, 0.714286, 0.27907, 0.618705, 0.733813, 0.878571]
x.append(100)
y.append(100)
z.append(1.0)
#f, ax = plt.subplots(1,2, sharex=True, sharey=True)
plt.figure(figsize=(7,6))
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
plt.clf()
bounds = [0.0,1.0]
ax = plt.tricontourf(x,y,z, 20, vmin=0.0, vmax=1.0) # choose 20 contour levels, just to show how good its interpolation is

plt.xlabel("Rate of trait evolution ($\\tau$)", fontsize=20)
plt.ylabel("Association strength ($\\lambda$)", fontsize=20)
plt.xlim(2.0,7.5)
plt.ylim(2.0,16.0)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#cb = plt.colorbar(shrink=1.0, boundaries=bounds)
cb = plt.colorbar()
plt.clim(0.0,1.0)
#cb.set_ticklabels(["0.0","0.25","0.5","0.75","1.0"])
#cb.set_ticks([0.0,0.25,0.5,0.75,1.0])
#cb.set_ticklabels()
cb.set_label(label="Recall",size=20)
cb.ax.tick_params(labelsize=20)
#plt.title("Sensitivity",size=24)
plt.savefig('benchmarks_contour.pdf')