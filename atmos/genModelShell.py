import numpy as np
import matplotlib.pyplot as pl

yearInSec = 365.0*24.0*3600.0
solarMassPerYear = 1.99e33 / yearInSec
RStar = 4e13
TStar = 2330.0
MStar = 0.8 * 1.99e33
R0 = 1.2 * RStar
Rc = 5 * RStar
Rw = 20.0 * RStar
vexp = 14.5 * 1e5
vturb = 1.0
MLoss = 2e-5 * solarMassPerYear
G = 6.67259e-8
k = 1.381e-16
mg = 2.3 * 1.6605402e-24
alpha = 0.55
nStar = 1.8e14
gamma = 0.89
dStar = 130  # pc
dShellIn = 15.0 * dStar / 206265.0 * 3.08e18
dShellOut = 19.0 * dStar / 206265.0 * 3.08e18

n = [10,30,30,30,20]

B = np.ones(sum(n)) * 0.1

vexp1 = np.ones(n[0]) * 0.0
vexp2 = np.ones(n[1]) * 5.0 * 1e5
vexp3 = np.ones(n[2]) * 11.0 * 1e5
vexp4 = np.ones(n[3]) * 14.5 * 1e5
vexp5 = np.ones(n[4]) * 14.5 * 1e5

r1 = np.linspace(RStar,R0,n[0])
r2 = np.linspace(R0,Rc,n[1])
r3 = np.linspace(Rc,Rw,n[2])
r4 = np.linspace(Rw,dShellIn,n[3])
r5 = np.linspace(dShellIn,dShellOut,n[4])

n3 = MLoss / (4.0*np.pi*r3**2 * mg * vexp3)
n4 = MLoss / (4.0*np.pi*r4**2 * mg * vexp4)
n5 = np.ones(n[4]) * 5e4
TDust3 = 800.0 * (r3/Rc)**(-0.375)
TDust4 = np.zeros(n[3])
TDust5 = np.zeros(n[4])

t = G * MStar * mg / (k * TStar * RStar * (1.0-alpha)) * (1.0 - (RStar/r1)**(1-alpha))
n1 = nStar * (r1 / RStar)**alpha * np.exp(-t)
TDust1 = np.zeros(n[0])

t = G * MStar * mg * (1.0-gamma**2) / (k * TStar * RStar**alpha * R0**(1.0-alpha) * (1.0-alpha)) * (1.0 - (r2/R0)**(alpha-1.0))
n2 = n1[-1] * np.exp(-t)
TDust2 = np.zeros(n[1])

r = np.concatenate((r1,r2,r3,r4,r5))
n = np.concatenate((n1,n2,n3,n4,n5))
TDust = np.concatenate((TDust1,TDust2,TDust3,TDust4,TDust5))
v = np.concatenate((vexp1,vexp2,vexp3,vexp4,vexp5))

u, order = np.unique(r, return_index=True)

CNAbundance = np.zeros(r.size)
CNAbundance[(r < 2.0*RStar)] = 8e-6
CNAbundance[( (r > 2.0*RStar) & (r < 2.5*RStar) )] = 1e-6
CNAbundance[( (r > 2.5*RStar) & (r < dShellIn))] = 0.0
CNAbundance[( (r > dShellIn) )] = 1e-5

pl.close('all')
fig, ax1 = pl.subplots()
ax1.loglog(r, n, 'b')
# ax1.loglog(r, CNAbundance * n, 'g')
ax1.set_xlabel('r [cm]')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel(r'n [H$_2$]')
#for tl in ax1.get_yticklabels():
    #tl.set_color('b')
ax1.set_ylim([1e4,1e15])

Tk = TStar * (r/RStar)**(-0.55)
Tk[(r > dShellIn)] = 20.0

ax2 = ax1.twinx()
ax2.plot(r, Tk, 'r')
ax2.set_ylabel('T [K]')
#for tl in ax2.get_yticklabels():
    #tl.set_color('r')
ax1.set_xlim([4e13,4e16])
ax2.set_ylim([-200,2400])
ax1.axvline(R0, color='k', linestyle='--')
ax1.axvline(Rc, color='k', linestyle='--')
ax1.axvline(Rw, color='k', linestyle='--')
ax1.annotate(r'$v_\mathrm{exp}$=5 km s$^{-1}$', (6e13,1e14), size='large')
ax1.annotate(r'$v_\mathrm{exp}$=11 km s$^{-1}$', (2.5e14,5e13), size='large')
ax1.annotate(r'$v_\mathrm{exp}$=14.5 km s$^{-1}$', (1e15,1e13), size='large')
ax1.annotate(r'$R_0$', (5e13,5e5), size='x-large')
ax1.annotate(r'$R_c$', (2.2e14,5e5), size='x-large')
ax1.annotate(r'$R_w$', (9e14,5e5), size='x-large')
pl.tight_layout()
pl.savefig('model.pdf')

f = open('modelWithShell100mG1e-5.atmos', 'w')
f.write("r [cm]      n[cm^-3]     n(CN) [cm^-3]  Tk [K]    Tdust[K]     v[km/s]    B[G]\n")
f.write(str(r.size)+'\n')
for j in range(order.size):
	i = order[j]
	f.write("{0:10.3e}  {1:10.3e}  {2:10.3e}  {3:10.3f}  {4:10.3f}  {5:10.3f}  {6:10.3f}\n".format(r[i], n[i], CNAbundance[i], Tk[i], TDust[i], v[i]*1e-5, B[i]))
f.close()
