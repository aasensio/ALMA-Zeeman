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
pc = 3.0857e18

n = 150

radius = 414.0 * pc * (10.0 / 206265.0)

r = np.linspace(0.001*radius, radius, n)
nH2 = np.ones(n) * 1.5e7
SOAbundance = np.ones(n) * 20e15 / 4.2e23
Tk = np.ones(n) * 220.0
TDust = np.zeros(n)
v = np.zeros(n)

B = np.ones(n) * 10.0 * 1e-3

f = open('model10mG_hotcore.atmos', 'w')
f.write("r [cm]      n[cm^-3]     n(CN) [cm^-3]  Tk [K]    Tdust[K]     v[km/s]    B[G]\n")
f.write("{0}\n".format(n))
for i in range(n):
	f.write("{0:10.3e}  {1:10.3e}  {2:10.3e}  {3:10.3f}  {4:10.3f}  {5:10.3f}  {6:10.3f}\n".format(r[i], nH2[i], SOAbundance[i], Tk[i], TDust[i], v[i]*1e-5, B[i]))
f.close()


radius = 414.0 * pc * (20.0 / 206265.0)

r = np.linspace(0.001*radius, radius, n)
nH2 = np.ones(n) * 5e6
SOAbundance = np.ones(n) * 69e15 / 2.1e23
Tk = np.ones(n) * 150.0
TDust = np.zeros(n)
v = np.zeros(n)

B = np.ones(n) * 10.0 * 1e-3

f = open('model10mG_plateau.atmos', 'w')
f.write("r [cm]      n[cm^-3]     n(CN) [cm^-3]  Tk [K]    Tdust[K]     v[km/s]    B[G]\n")
f.write("{0}\n".format(n))
for i in range(n):
    f.write("{0:10.3e}  {1:10.3e}  {2:10.3e}  {3:10.3f}  {4:10.3f}  {5:10.3f}  {6:10.3f}\n".format(r[i], nH2[i], SOAbundance[i], Tk[i], TDust[i], v[i]*1e-5, B[i]))
f.close()

radius = 414.0 * pc * (3.0 / 206265.0)

r = np.linspace(0.001*radius, radius, n)
nH2 = np.ones(n) * 1.5e7
SOAbundance = np.ones(n) * 20e15 / 4.2e23
Tk = np.ones(n) * 220.0
TDust = np.zeros(n)
v = np.zeros(n)

B = np.ones(n) * 10.0 * 1e-3

f = open('model10mG_hotcorecentral.atmos', 'w')
f.write("r [cm]      n[cm^-3]     n(CN) [cm^-3]  Tk [K]    Tdust[K]     v[km/s]    B[G]\n")
f.write("{0}\n".format(n))
for i in range(n):
    f.write("{0:10.3e}  {1:10.3e}  {2:10.3e}  {3:10.3f}  {4:10.3f}  {5:10.3f}  {6:10.3f}\n".format(r[i], nH2[i], SOAbundance[i], Tk[i], TDust[i], v[i]*1e-5, B[i]))
f.close()