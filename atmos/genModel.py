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

inputModel = np.loadtxt('rpfit_iktau.dat', skiprows=3)
n = inputModel.shape[0]

B = np.ones(n) * 1.0
r = inputModel[:,0]
nH2 = inputModel[:,1]
SOAbundance = inputModel[:,11]
Tk = inputModel[:,2]
TDust = inputModel[:,5]
v = inputModel[:,4]

f = open('model1G.atmos', 'w')
f.write("r [cm]      n[cm^-3]     A(mol)     Tk [K]    Tdust[K]     v[km/s]    B[G]\n")
f.write("{0}\n".format(n))
for i in range(n):
	f.write("{0:10.3e}  {1:10.3e}  {2:10.3e}  {3:10.3f}  {4:10.3f}  {5:10.3f}  {6:10.3f}\n".format(r[i], nH2[i], SOAbundance[i], Tk[i], TDust[i], v[i], B[i]))
f.close()

v = inputModel[:,4] * 0.0

f = open('model1G_rest.atmos', 'w')
f.write("r [cm]      n[cm^-3]     A(mol)       Tk [K]    Tdust[K]     v[km/s]    B[G]\n")
f.write("{0}\n".format(n))
for i in range(n):
    f.write("{0:10.3e}  {1:10.3e}  {2:10.3e}  {3:10.3f}  {4:10.3f}  {5:10.3f}  {6:10.3f}\n".format(r[i], nH2[i], SOAbundance[i], Tk[i], TDust[i], v[i], B[i]))
f.close()
