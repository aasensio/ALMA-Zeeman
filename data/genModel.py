import numpy as np
from ipdb import set_trace as stop

"""
Generate a model molecule from the linelist
"""

def zeemanLande(S, N, J):
	"""
	Return the Zeeman splitting factor in Hz/microG
	"""	
	if (J == 0):
		gJ = 0.0
	else:
		gJ = (J*(J+1)+S*(S+1)-N*(N+1)) / (J*(J+1))
		
	return gJ

dat = np.loadtxt('so.linelist',skiprows=2)
# datv1 = np.loadtxt('sov1.linelist',skiprows=2)

# dat = np.vstack([dat,datv1])

nTransitions = dat.shape[0]
nLevels = np.max(dat[:,0])

energy = np.loadtxt('so.energy', skiprows=2)
# energyv1 = np.loadtxt('sov1.energy', skiprows=2)
# energyv1[:,3] += 1637.41
# energy = np.vstack([energy, energyv1])
nMax = int(np.max(energy[:,1:3]))
lev = np.zeros((nMax+1,nMax+1))
E = np.zeros((nMax+1,nMax+1))
for i in range(energy.shape[0]):
	lev[int(energy[i,1]),int(energy[i,2])] = int(energy[i,0])
	E[int(energy[i,1]),int(energy[i,2])] = energy[i,3] / 6.62606876e-27 * 1.3806503e-16

print("up   low     nu [Hz]          Elow [Hz]  Aul [Hz]      gu        gl      Ju  Jl     geff")
with open('soFinal.linelist','w') as f:
	f.write("up   low     nu [Hz]          Elow [Hz]  Aul [Hz]      gu        gl      Ju  Jl\n")
	for i in range(nTransitions):
		levUp = lev[int(dat[i,1]),int(dat[i,2])]
		levLow = lev[int(dat[i,3]),int(dat[i,4])]
		nu = dat[i,5] * 1e6
		Elow = E[int(dat[i,3]),int(dat[i,4])]
		Aul = dat[i,8]
		
		Nu, Ju, Nl, Jl = dat[i,1:5]
				
		gu = zeemanLande(1.5,Nu,Ju)
		gl = zeemanLande(1.5,Nl,Jl)

		geff = 0.5*(gu+gl + 0.25 *(gu-gl*(Ju*(Ju+1)-Jl*(Jl+1))))

		print("{0:3d} {1:3d} {2:11.3f} {3:14.5f} {4:10.3e} {5:10.5f} {6:10.5f} {7:3.1f} {8:3.1f} {9:10.5f}".format(int(levUp),int(levLow),nu,Elow,Aul,gu,gl,Ju,Jl,geff))

		f.write("{0:3d} {1:3d} {2:11.3f} {3:14.5f} {4:10.3e} {5:10.5f} {6:10.5f} {7:3.1f} {8:3.1f}\n".format(int(levUp),int(levLow),nu,Elow,Aul,gu,gl,Ju,Jl))
	f.close()

with open('soFinal.energy','w') as f:
	f.write("ind      E [Hz]               E [K]      g\n")
	for i in range(int(energy.shape[0])):
		f.write("{0:3d}   {1:13.6e}    {2:14.5f}     {3:3d}\n".format(i, energy[i,3] / 6.62606876e-27 * 1.3806503e-16, energy[i,3], int(energy[i,4])))
	f.close()