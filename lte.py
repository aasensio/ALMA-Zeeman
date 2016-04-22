import numpy as np
import sys
from lowerToSep import lowerToSep
from configobj import ConfigObj
import wigner
import pydust
from profileVoigtFaraday import profileVoigtFaraday
import pydelopar
import matplotlib.pyplot as pl

class lteSynthesizer(object):
	
#--------------------------------		
# Constructor
#--------------------------------
	def __init__(self, config):
		
		self.eMass = 9.10938188e-28
		self.eCharge = 4.8032e-10
		self.lightSpeed = 2.99792458e10
		self.planck = 6.62606876e-27
		self.boltzmann = 1.3806503e-16
		self.microturbulence = float(config['general']['microturbulence']) * 1e5
		self.amu = 1.66053873e-24
		self.sqrtpi = np.sqrt(np.pi)
		self.hc2 = 2.0 * self.planck / self.lightSpeed**2
		self.hk = self.planck / self.boltzmann
		self.dustToGas = 1.0 / 500.0
		
		self.fileLinelist = config['general']['file with linelist']
		self.fileEnergy = config['general']['file with energies']
		self.readLineList()
		
		self.fileModelAtm = config['general']['file with model atmosphere']
		self.readAtmosphere()
		
		self.nRegions = int(config['wavelength regions']['number of regions'])		
		self.ranges = np.zeros((3,self.nRegions))

		self.fieldTopology = config['general']['field topology']

		self.typeSource = config['general']['type']
		self.distance = float(config['general']['distance'])
		self.impactParameters = np.asarray(config['general']['impact parameters'], dtype='float')

		self.impactParameters *= self.distance / 206265.0 * 3.08e18

		self.CNMass = self.amu * float(config['general']['mass'])

		print("Molecular mass -> {0} AMU".format(float(config['general']['mass'])))
		
# Read the regions
		for i in range(self.nRegions):
			s = "region {0}".format(i)
			info = config['wavelength regions'][s]
			print("Region {0} -> [{1},{2}] MHz with {3} points".format(i, info[0], info[1], info[2]))
			self.ranges[:,i] = np.asarray(config['wavelength regions'][s])
			self.ranges[0:2,i] *= 1e9   # to Hz
			
		self.generateLineInfo()
				
		self.larmor = self.eCharge / (4.0*np.pi*self.eMass*self.lightSpeed)

		self.config = config
		
#--------------------------------		
# Read the atmosphere
#--------------------------------
	def readAtmosphere(self):
		self.atmosphere = np.loadtxt(self.fileModelAtm, skiprows=2)
		
#--------------------------------		
# Read the linelist
#--------------------------------
	def readLineList(self):
		self.linelist = np.loadtxt(self.fileLinelist, skiprows=1)
		energies = np.loadtxt(self.fileEnergy, skiprows=1)
		self.energy = energies[:,2][:,np.newaxis]
		self.gLevel = energies[:,3][:,np.newaxis]
		
#--------------------------------
# Extract the atmosphere for a given impact parameter
#--------------------------------
	def extractAtmosphere(self, impactParameter):
		
# Compute x position along the ray with the given impact parameter
		x = np.sqrt(self.atmosphere[:,0]**2 - impactParameter**2)
		ind = ~np.isnan(x)
				
		self.x = x[ind]
		self.nPoints = self.x.size
		
		self.nH2 = self.atmosphere[ind,1]
		self.nCN = self.atmosphere[ind,2] * self.nH2
		self.Tk = self.atmosphere[ind,3]
		self.TDust = self.atmosphere[ind,4]
		self.vBulk = self.atmosphere[ind,5] * self.x / self.atmosphere[ind,0] * 1e5
		self.BField = self.atmosphere[ind,6]
		if (self.fieldTopology == 'radial'):
			self.cosGamma = (self.x / self.atmosphere[ind,0])
			self.sinGamma = np.sqrt(1.0 - self.cosGamma**2)
		if (self.fieldTopology == 'perpendicular'):
			self.cosGamma = np.zeros(self.nPoints)
			self.sinGamma = np.ones(self.nPoints)
		self.cos2Chi = np.ones(self.nPoints)
		self.sin2Chi = np.zeros(self.nPoints)
		self.aDamp = np.zeros(self.nPoints)
				
		
# Test whether the ray impacts the stellar surface
# If not, then use zero as boundary condition and replicate the physical conditions
		if (self.typeSource == 'cloud'):
			self.boundaryCondition = 'cmb'
			self.x = np.concatenate((-self.x[::-1], self.x[1:]))
			self.nH2 = np.concatenate((self.nH2[::-1], self.nH2[1:]))
			self.nCN = np.concatenate((self.nCN[::-1], self.nCN[1:]))
			self.Tk = np.concatenate((self.Tk[::-1], self.Tk[1:]))
			self.TDust = np.concatenate((self.TDust[::-1], self.TDust[1:]))
			self.vBulk = np.concatenate((-self.vBulk[::-1], self.vBulk[1:]))
			self.BField = np.concatenate((self.BField[::-1], self.BField[1:]))
			self.cosGamma = np.concatenate((-self.cosGamma[::-1], self.cosGamma[1:]))
			self.sinGamma = np.sqrt(1.0 - self.cosGamma**2)
			self.nPoints = self.x.size
		
			self.cos2Chi = np.ones(self.nPoints)
			self.sin2Chi = np.zeros(self.nPoints)
			self.aDamp = np.zeros(self.nPoints)

		if (self.typeSource == 'star'):
			self.boundaryCondition = 'star'
			if (impactParameter > self.atmosphere[0,0]):
				self.boundaryCondition = 'cmb'			
				self.x = np.concatenate((-self.x[::-1], self.x[1:]))
				self.nH2 = np.concatenate((self.nH2[::-1], self.nH2[1:]))
				self.nCN = np.concatenate((self.nCN[::-1], self.nCN[1:]))
				self.Tk = np.concatenate((self.Tk[::-1], self.Tk[1:]))
				self.TDust = np.concatenate((self.TDust[::-1], self.TDust[1:]))
				self.vBulk = np.concatenate((-self.vBulk[::-1], self.vBulk[1:]))
				self.BField = np.concatenate((self.BField[::-1], self.BField[1:]))
				self.cosGamma = np.concatenate((-self.cosGamma[::-1], self.cosGamma[1:]))
				self.sinGamma = np.sqrt(1.0 - self.cosGamma**2)
				self.nPoints = self.x.size
			
				self.cos2Chi = np.ones(self.nPoints)
				self.sin2Chi = np.zeros(self.nPoints)
				self.aDamp = np.zeros(self.nPoints)
			
		self.sinGamma = self.sinGamma[:,np.newaxis]
		self.cosGamma = self.cosGamma[:,np.newaxis]
		self.cos2Chi = self.cos2Chi[:,np.newaxis]
		self.sin2Chi = self.sin2Chi[:,np.newaxis]			

#--------------------------------		
# Compute splitting and strength of each transition. This works for fine and hyperfine structure in the Zeeman regime
#--------------------------------
	def computeSplittingStrength(self, JUp, JLow, gu, gl):
		nUp = 2*JUp+1
		nLow = 2*JLow+1
		
		iPi = 0
		iBlue = 0
		iRed = 0
		
		nComponents = np.zeros(3)
		for iUp in range(1,int(nUp)+1):
			MUp = JUp + 1 - iUp
			for iLow in range(1,4):
				MLow = MUp - 2 + iLow
				if (abs(MLow) <= JLow):
					nComponents[iLow-1] += 1
		
		splitting = np.zeros((3,int(np.max(nComponents))))
		strength = np.zeros((3,int(np.max(nComponents))))
		
		for iUp in range(1,int(nUp)+1):
			MUp = JUp + 1 - iUp
			for iLow in range(1,4):
				MLow = MUp - 2 + iLow
				if (abs(MLow) <= JLow):
					if (iLow == 1):
						which = iBlue
						iBlue += 1						
						
					if (iLow == 2):
						which = iPi
						iPi += 1						
						
					if (iLow == 3):
						which = iRed
						iRed += 1						
					
					strength[iLow-1,which] = 3.0 * wigner.Wigner3j(JUp,JLow,1,MUp,-MLow,MLow-MUp)**2
					splitting[iLow-1,which] = gu*MUp - gl*MLow
					
		return nComponents, splitting, strength
		
#--------------------------------		
# Filter the linelist
# Go through each line and accept it if is inside a region and compute their main properties
#--------------------------------
	def generateLineInfo(self):
		self.lineInfo = []		
		self.linesInRegion = [[] for x in range(self.nRegions)]
		loop = 0
		for i in range(self.linelist.shape[0]):
			info = {}
			accept = False			
			for j in range(self.nRegions):				
				if ((self.linelist[i,2] > self.ranges[0,j]) and (self.linelist[i,2] < self.ranges[1,j])):
					accept = True
					info['region'] = j
			if (accept):
				info['up'] = self.linelist[i,0]
				info['low'] = self.linelist[i,1]
				info['freq'] = self.linelist[i,2]
				info['elow'] = self.linelist[i,3] * self.planck / self.boltzmann
				info['aul'] = self.linelist[i,4]				
				info['landeu'] = self.linelist[i,5]
				info['landel'] = self.linelist[i,6]
				info['ju'] = self.linelist[i,7]
				info['jl'] = self.linelist[i,8]
				info['hnu_k'] = self.linelist[i,2] * self.planck / self.boltzmann
				
				gu = 2.0*info['ju']+1.0
				info['gl_blu'] = self.lightSpeed**2 / (2.0 * self.planck * info['freq']**3) * gu * info['aul']

				info['geff'] = 0.5*(info['landeu'] + info['landel'])+0.25*(info['landeu'] - info['landel']) * (info['ju']*(info['ju']+1.0)-info['jl']*(info['jl']+1.0))
				
				info['ncomponents'], info['splitting'], info['strength'] = self.computeSplittingStrength(info['ju'],info['jl'],info['landeu'],info['landel'])
				
				self.lineInfo.append(info)
				self.linesInRegion[info['region']].append(loop)
				
				loop += 1
		print("{0} lines included".format(len(self.lineInfo)))
		
#--------------------------------		
# Compute the Zeeman profiles taking into account the splitting of the lines
#--------------------------------
	def zeemanProfile(self, whichLine, BField, vBulk, deltaNu, aDamp, frequency):
		iPi = 0
		iBlue = 0
		iRed = 0
		
		JUp = self.lineInfo[whichLine]['ju']
		JLow = self.lineInfo[whichLine]['jl']
		
		nUp = 2*JUp+1
		nLow = 2*JLow+1
		
		nz = BField.size
		nf = frequency.size
		zeemanVoigt = np.zeros((3,nz,nf))
		zeemanFaraday = np.zeros((3,nz,nf))
				
		for iUp in range(1,int(nUp)+1):
			MUp = JUp + 1 - iUp
			for iLow in range(1,4):
				MLow = MUp - 2 + iLow
				if (abs(MLow) <= JLow):
					if (iLow == 1):
						which = iBlue
						iBlue += 1						
						
					if (iLow == 2):
						which = iPi
						iPi += 1						
						
					if (iLow == 3):
						which = iRed
						iRed += 1						
						
					strength = self.lineInfo[whichLine]['strength'][iLow-1,which]
					splitting = self.larmor * BField * self.lineInfo[whichLine]['splitting'][iLow-1,which]
										
					lineFreq = self.lineInfo[whichLine]['freq']
																									
					voigt, faraday = profileVoigtFaraday(((lineFreq - frequency[np.newaxis,:] - lineFreq * vBulk[:,np.newaxis] / self.lightSpeed + splitting[:,np.newaxis]) / deltaNu[:,np.newaxis]), aDamp[:,np.newaxis])
					
					zeemanVoigt[iLow-1,:,:] += strength * voigt
					zeemanFaraday[iLow-1,:,:] += strength * faraday		
					
		
		return zeemanVoigt, zeemanFaraday

#--------------------------------		
# Compute the eta_i and rho_i terms from the Voigt and Faraday profiles and the geometry of the magnetic field
#--------------------------------
	def zeemanOpacity(self, zeemanVoigt, zeemanFaraday, sinGamma, cosGamma, sin2Chi, cos2Chi):
				
		ni, nz, nf = zeemanVoigt.shape
		
		coefficients = np.zeros((nz,7,nf))
# Dichroic terms
		coefficients[:,0,:] = 0.5 * (zeemanVoigt[1,:,:]*sinGamma**2 + 0.5*(zeemanVoigt[0,:,:]+zeemanVoigt[2,:,:])*(1.0+cosGamma**2))  # eta_I
		coefficients[:,1,:] = 0.5 * (zeemanVoigt[1,:,:] - 0.5*(zeemanVoigt[0,:,:]+zeemanVoigt[2,:,:])) * sinGamma**2*cos2Chi  # eta_Q
		coefficients[:,2,:] = 0.5 * (zeemanVoigt[1,:,:] - 0.5*(zeemanVoigt[0,:,:]+zeemanVoigt[2,:,:])) * sinGamma**2*sin2Chi  # eta_U
		coefficients[:,3,:] = 0.5 * (zeemanVoigt[2,:,:]-zeemanVoigt[0,:]) * cosGamma  # eta_V

# Magneto-optical coefficients
		coefficients[:,4,:] = 0.5 * (zeemanFaraday[1,:,:] - 0.5*(zeemanFaraday[0,:,:]+zeemanFaraday[2,:,:])) * sinGamma**2*cos2Chi  # rho_Q
		coefficients[:,5,:] = 0.5 * (zeemanFaraday[1,:,:] - 0.5*(zeemanFaraday[0,:,:]+zeemanFaraday[2,:,:])) * sinGamma**2*sin2Chi  # rho_U
		coefficients[:,6,:] = 0.5 * (zeemanFaraday[2,:,:]-zeemanFaraday[0,:,:]) * cosGamma  # rho_V
		
		return coefficients
	
#--------------------------------		
# Compute partition function
#--------------------------------
	def partitionFunction(self, T):
		out = np.zeros_like(T)
		cut = 5000.0
		out[T < cut] = np.sum(self.gLevel * np.exp(-self.energy / T[T < cut]),axis=0)
		out[T > cut] = 10.0**(4.7963 - 2.1308 * np.log10(5040./T[T>cut]) + 0.5224 * np.log10(5040./T[T>cut])**2)
		return out
				
#--------------------------------		
# Add line opacity
#--------------------------------
	def addLineOpacity(self, line):
		print('Adding line {0} at frequency {1} GHz - geff={2}'.format(line, self.lineInfo[line]['freq']/1e9, self.lineInfo[line]['geff']))
		constant = self.planck * self.lineInfo[line]['freq'] / (4.0*np.pi) * self.lineInfo[line]['gl_blu'] / self.partitionFunction(self.Tk)				
		opacity = constant * np.exp(-self.lineInfo[line]['elow'] / self.Tk) * (1.0 - np.exp(-self.lineInfo[line]['hnu_k'] / self.Tk)) * self.nCN
		return opacity		
	
#--------------------------------		
# Add background dust opacity for grains of 0.1 micron and density 2 g/cm^3
#--------------------------------
	def addBackgroundOpacity(self, line):
		return pydust.kappaDust(self.lineInfo[line]['freq'], 1, 0.1) * self.nH2 * self.amu * self.dustToGas * 2.0
	
#--------------------------------		
# Planck function
#--------------------------------
	def planckFunction(self, T, frequency):
		out = np.zeros_like(T)
		if (np.isscalar(T)):
			out = self.hc2 * frequency**3 / (np.exp(self.hk * frequency / T) - 1.0)
		else:
			out[T > 0] = self.hc2 * frequency**3 / (np.exp(self.hk * frequency / T[T>0]) - 1.0)
		return out

#--------------------------------		
# Inverse Planck function
#--------------------------------
	def invPlanckFunction(self, I, frequency):
		out = np.zeros_like(I)
		if (np.isscalar(I)):
			out = self.hk * frequency / np.log(1.0 + self.hc2 * frequency**3 / I)
		else:
			out[I > 0] = self.hk * frequency / np.log(1.0 + self.hc2 * frequency**3 / I[I>0])
		return out

		
#--------------------------------
# Synthesize a given impact parameter
#--------------------------------
	def synthesizeLines(self, impactParameter):
		self.extractAtmosphere(impactParameter)
		
		stokes = []
		for region in range(self.nRegions):
			print("Synthesizing region {0}".format(region))
			stokesLocal = np.zeros((5,int(self.ranges[2,region])))
			
# Generate frequency axis for the region
			stokesLocal[0,:] = self.ranges[0,region] + np.arange(self.ranges[2,region]) / (self.ranges[2,region]-1.0) * (self.ranges[1,region] - self.ranges[0,region])
			
			opacity = np.zeros((self.nPoints,7,int(self.ranges[2,region])))
			emissivity = np.zeros((self.nPoints,4,int(self.ranges[2,region])))
			
# Go line by line computing the opacity
			for line in self.linesInRegion[region]:
				opacityL = self.addLineOpacity(line)
				opacityB = self.addBackgroundOpacity(line)
							
				boundary = np.asarray([0.0,0.0,0.0,0.0])
								
# Doppler width of the line
				doppler = np.sqrt(self.microturbulence**2 + 2.0*self.boltzmann*self.Tk / self.CNMass)
				deltaNu = self.lineInfo[line]['freq'] * doppler / self.lightSpeed
				
# Zeeman and Faraday profiles for each height and each frequency
				zeemanVoigt, zeemanFaraday = self.zeemanProfile(line, self.BField, self.vBulk, deltaNu, self.aDamp, stokesLocal[0,:])
												
# Coefficients of the propagation matrix
				coefficients = self.zeemanOpacity(zeemanVoigt, zeemanFaraday, self.sinGamma, self.cosGamma, self.sin2Chi, self.cos2Chi)
				
				opacity += (opacityL / (deltaNu * self.sqrtpi))[:,np.newaxis,np.newaxis] * coefficients

				planckGas = self.planckFunction(self.Tk, self.lineInfo[line]['freq'])
				planckDust = self.planckFunction(self.TDust, self.lineInfo[line]['freq'])
								
# Emissivity of the gas
				emissivity[:,0:4,:] = opacity[:,0:4,:] * planckGas[:,np.newaxis,np.newaxis]
				
# Emissivity of the dust
				emissivity[:,0,:] += opacityB[:,np.newaxis] * planckDust[:,np.newaxis]
				
# Add continuum opacity
				opacity[:,0,:] += opacityB[:,np.newaxis]
												
				for i in range(int(self.ranges[2,region])):
					if (self.boundaryCondition == 'star'):
						boundary[0] = self.planckFunction(self.Tk[0], stokesLocal[0,i])	
					if (self.boundaryCondition == 'cmb'):
						boundary[0] = self.planckFunction(2.73, stokesLocal[0,i])									
					
					stokesLocal[1:,i] = pydelopar.deloparFormal(self.x, opacity[:,:,i], emissivity[:,:,i], boundary)

			stokes.append(stokesLocal)
			
		return stokes
	
#--------------------------------
# Synthesize many impact parameters
#--------------------------------
	def synthesizeImpact(self):		
		n = len(self.impactParameters)
		
		self.stokesTotal = []
		for i in range(n):
			stokes = self.synthesizeLines(self.impactParameters[i])
			self.stokesTotal.append(stokes)
			
#--------------------------------
# Synthesize many impact parameters and integrate
#--------------------------------
	def synthesizeAsStar(self):
		
		impactParameters = np.concatenate((np.atleast_1d(0),self.atmosphere[:,0]))
		n = len(impactParameters)
		
# Total area of the star
		areaTotal = np.pi * impactParameters[-1]**2
		
# Synthesize all impact parameters
		self.stokesTotal = []
		for i in range(n):
			stokes = self.synthesizeLines(impactParameters[i])
			self.stokesTotal.append(stokes)
										
# Add them weighted by the area associated to each impact parameter
		self.stokesStar = []
		for region in range(self.nRegions):
			stokes = np.zeros(self.stokesTotal[0][0].shape)
			for i in range(1,n):
				area = np.pi*(impactParameters[i]**2 - impactParameters[i-1]**2)							
				fun = (self.stokesTotal[i][region][1:,:] + self.stokesTotal[i-1][region][1:,:]) * 0.5
				stokes[0,:] = self.stokesTotal[i][region][0,:]
				stokes[1:,:] += area * fun / areaTotal
			self.stokesStar.append(stokes)
				
#--------------------------------
# Plot the results line by line
#--------------------------------
	def plotLines(self, impactParameters):
		n = len(impactParameters)

		stop()
						
		stokes = ['I','Q','U','V']
		label = ['Intensity [K]','Q/I [%]','U/I [%]','V/I [%]']
		for indStokes in range(4):
			
			pl.close('all')
			
			fig, ax = pl.subplots(nrows=2, ncols=2, figsize=(8,10), sharex="col")
					
			ax = ax.flatten()
			
			loop = 0
						
			for i in range(2):
				for j in range(2):
					if (loop < len(self.lineInfo)):
						for k in range(n):
							region = self.lineInfo[loop]['region']
							v = (self.stokesTotal[k][region][0,:] - self.lineInfo[loop]['freq']) * self.lightSpeed / self.lineInfo[loop]['freq'] * 1e-5
							solidAngle = 2.0 * np.pi * (1.0 - np.cos(19.0 / 206265.0))
							factor = (self.lightSpeed / self.lineInfo[loop]['freq'])**2 / (2.0 * self.boltzmann)
							
							if (indStokes == 0):
								ax[loop].plot(v, factor * self.stokesTotal[k][region][indStokes+1,:], label='p='+str(impactParameters[k]))
							else:
								ax[loop].plot(v, 100.0*self.stokesTotal[k][region][indStokes+1,:] / np.max(self.stokesTotal[k][region][1,:]))
								
							ax[loop].set_title("{0:7.3f} GHz".format(self.lineInfo[loop]['freq']/1e9), fontsize=10)
						ax[loop].set_xlim((-10,10))
					
					if (i == 4):
						ax[loop].set_xlabel("v [km/s]")
						
					for l in (ax[loop].get_xticklabels()+ax[loop].get_yticklabels()):
						l.set_fontsize(8)
					loop += 1
								
			ax[0].set_ylabel(label[indStokes])
			ax[0].legend(loc='center left', fontsize=6)
			#pl.tight_layout()
			pl.savefig('{0}.{1}.pdf'.format(self.config['general']['output file'], stokes[indStokes]))

#--------------------------------
# Plot the results line by line
#--------------------------------
	def plotRegions(self):
		n = len(self.impactParameters)
						
		stokes = ['I','Q','U','V']
		label = ['Intensity [K]','Q/Imax [%]','U/Imax [%]','V/Imax [%]']

		pl.close('all')
			
		fig, ax = pl.subplots(nrows=4, ncols=self.nRegions, figsize=(17,10), sharex="col")

		if (self.nRegions == 1):
			ax = np.transpose(np.atleast_2d(ax))
					
		for i in range(self.nRegions):
			for indStokes in range(4):
				for k in range(n):					
					freq = np.mean(self.ranges[0:2,i])

					if (indStokes == 0):
						ax[indStokes,i].plot(self.stokesTotal[k][i][0,:] / 1e9, self.invPlanckFunction(self.stokesTotal[k][i][indStokes+1,:], freq), label='p='+str(self.impactParameters[k]))
					else:
						ax[indStokes,i].plot(self.stokesTotal[k][i][0,:] / 1e9, 100.0*self.stokesTotal[k][i][indStokes+1,:] / np.max(self.stokesTotal[k][i][1,:]))						
						# ax[indStokes,i].set_title("{0:7.3f} GHz".format(self.lineInfo[loop]['freq']/1e9), fontsize=10)
						# ax[loop].set_xlim((-10,10))
					
					if (indStokes == 3):
						ax[indStokes,i].set_xlabel("Frequency [GHz]")
						
					for l in (ax[indStokes,i].get_xticklabels()+ax[indStokes,i].get_yticklabels()):
						l.set_fontsize(8)
					# loop += 1
								
				ax[indStokes,0].set_ylabel(label[indStokes])
			ax[0,0].legend(loc='center left', fontsize=6)
			pl.tight_layout()
			pl.savefig('{0}.pdf'.format(self.config['general']['output file']))
			
#--------------------------------
# Plot the results line by line
#--------------------------------
	def plotLinesAsStar(self):
		

		stokes = ['I','Q','U','V']
		label = ['Intensity [10$^{16}$ cgs]','Q/I [%]','U/I [%]','V/I [%]']
		for indStokes in range(4):
			
			pl.close('all')
			
			fig, ax = pl.subplots(nrows=5, ncols=4, figsize=(8,10), sharex="col")
					
			ax = ax.flatten()
						
			loop = 0			
						
			for i in range(4):
				for j in range(5):
					if (loop < len(self.lineInfo)):
						region = self.lineInfo[loop]['region']
						v = (self.stokesStar[region][0,:] - self.lineInfo[loop]['freq']) * self.lightSpeed / self.lineInfo[loop]['freq'] * 1e-5
						solidAngle = 2.0 * np.pi * (1.0 - np.cos(19.0 / 206265.0))
						factor = (self.lightSpeed / self.lineInfo[loop]['freq'])**2 / (2.0 * self.boltzmann * solidAngle)
						stop()
						
						if (indStokes == 0):
							ax[loop].plot(v, self.stokesStar[region][indStokes+1,:] * 1e16)
						else:
							ax[loop].plot(v, 100.0*self.stokesStar[region][indStokes+1,:] / self.stokesStar[region][1,:])
								
						ax[loop].set_title("{0:7.3f} GHz".format(self.lineInfo[loop]['freq']/1e9), fontsize=10)
					
					ax[loop].set_xlim((-10,10))
					
					if (i == 2):
						ax[loop].set_xlabel("v [km/s]")
						
					for l in (ax[loop].get_xticklabels()+ax[loop].get_yticklabels()):
						l.set_fontsize(8)
					loop += 1
								
			ax[0].set_ylabel(label[indStokes])
			pl.tight_layout()
			#pl.savefig('figs/{0}.{1}.star.pdf'.format(self.config['general']['output file'], stokes[indStokes]))

	
if (len(sys.argv) < 2):
	print("You need to give the configuration file")
	sys.exit()
	
print("Using configuration file {0}".format(sys.argv[1]))

oldErrorSettings = np.geterr()
np.seterr(invalid = 'ignore')

f = open(sys.argv[1], 'r')
inputLines = f.readlines()
f.close()
inputLower = ['']
for l in inputLines:
	inputLower.append(lowerToSep(l)) # Convert keys to lowercase
	
config = ConfigObj(inputLower)


lte = lteSynthesizer(config)
lte.synthesizeImpact()
lte.plotRegions()

np.seterr(**oldErrorSettings)