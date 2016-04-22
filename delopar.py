import numpy as np
import pdb

def fillPropagationMatrix(etaI, etaQ, etaU, etaV, rhoQ, rhoU, rhoV):
	"""
	Fill the propagation matrix generating a (4,4,nz) matrix
	"""
	absMat = np.zeros((4,4,etaI.size))
	absMat[0,0,:] = etaI
	absMat[1,1,:] = etaI
	absMat[2,2,:] = etaI
	absMat[3,3,:] = etaI
	
	absMat[0,1,:] = etaQ
	absMat[1,0,:] = etaQ
	
	absMat[0,2,:] = etaU
	absMat[2,0,:] = etaU
	
	absMat[0,3,:] = etaV
	absMat[3,0,:] = etaV
	
	absMat[1,2,:] = rhoV
	absMat[2,1,:] = -rhoV
	
	absMat[3,1,:] = rhoU
	absMat[1,3,:] = -rhoU
	
	absMat[2,3,:] = rhoQ
	absMat[3,2,:] = -rhoQ
	
	return absMat

def linearCoef(dtM):
	"""
	Return the coefficients for the linear interpolation of the source function for one frequency
	"""
	out = np.zeros(2)
	if (dtM != 0.0):
		if (dtM > 1e-4):
			exu = np.exp(-dtM)
			u0 = 1.0 - exu
			u1 = dtM - 1.0 + exu
			out[0] = u1 / dtM
			out[1] = u0 - out[0]			
		else:
			out[0] = 0.5 * dtM - dtM**2 / 6.0
			out[1] = 0.5 * dtM - dtM**2 / 3.0
			
	return out

def parabCoef(dtM, dtP):
	"""
	Return the coefficients for the parabolic interpolation of the source function for one frequency
	"""
	out = np.zeros(3)
	if (dtM != 0.0 and dtP != 0.0):
		if (dtM > 1e-4):
			exu = np.exp(-dtM)
			u0 = 1.0 - exu
			u1 = dtM - 1.0 + exu
			u2 = dtM**2 - 2.0*dtM + 2.0 - 2.0*exu
		else:			
			u0 = dtM - 0.5 * dtM**2
			u1 = 0.5 * dtM**2 - dtM**3 / 6.0
			u2 = dtM**3 / 3.0 - dtM**4 / 12.0
				
		out[0] = (u2 - u1*(dtP+2.0*dtM)) / (dtM*(dtM+dtP)) + u0
		out[1] = (u1*(dtM+dtP) - u2) / (dtM * dtP)
		out[2] = (u2 - u1*dtM) / (dtP*(dtM+dtP))
			
	return out
	
def deloparFormal(z, opacity, emissivity, boundary):
	"""
	delopar formal solution of the polarized radiative transfer equation.
	The inputs are the absorption, emission and magneto-optical effects terms
	and the boundary condition, plus the distance along the ray
	"""
	
	nz = z.size
	
	etaI = opacity[0,:]
	etaQ = opacity[1,:]
	etaU = opacity[2,:]
	etaV = opacity[3,:]
	rhoQ = opacity[4,:]
	rhoU = opacity[5,:]
	rhoV = opacity[6,:]
	

# Fill the propagation matrix and the source vector	
	absMat = fillPropagationMatrix(etaI, etaQ, etaU, etaV, rhoQ, rhoU, rhoV)
	
	sourceVec = np.zeros((4,nz))
	sourceVec[0,:] = emissivity[0,:] / etaI
	sourceVec[1,:] = emissivity[1,:] / etaI
	sourceVec[2,:] = emissivity[2,:] / etaI
	sourceVec[3,:] = emissivity[3,:] / etaI
	
# Transform K into K* and then into K'
	for i in range(4):
		for j in range(4):
			absMat[i,j,:] = absMat[i,j,:] / etaI
		absMat[i,i,:] = absMat[i,i,:] - 1.0
		
	stokesVec = boundary
	
	identity = np.identity(4)
	
	for k in range(1,nz):
		if (k != nz-1):
			km = k - 1
			kp = k + 1
			
			chiM = etaI[km]
			chiO = etaI[k]
			chiP = etaI[kp]
			
			sM = sourceVec[:,km]
			sO = sourceVec[:,k]
			sP = sourceVec[:,kp]
			
			dM = np.abs(z[k]-z[km])
			dP = np.abs(z[kp]-z[k])
		else:
			km = k - 1
			
			chiM = etaI[km]
			chiO = etaI[k]
			
			sM = sourceVec[:,km]
			sO = sourceVec[:,k]
			
			dM = np.abs(z[k]-z[km])
		
		dtM = 0.5*(chiM + chiO) * dM
		dtP = 0.5*(chiP + chiO) * dP
		
		if (dtM > 1e-4):
			exu = np.exp(-dtM)
		else:
			exu = 1.0 - dtM + 0.5*dtM**2
		
		psiLin = linearCoef(dtM)
		
		mat1 = exu * identity - psiLin[0] * absMat[:,:,km]
		mat2 = identity + psiLin[1] * absMat[:,:,k]
		mat2 = np.linalg.inv(mat2)
		
		if (k != nz-1):
			psi = parabCoef(dtM, dtP)
			t1 = psi[0] * sM + psi[1] * sO + psi[2] * sP
		else:
			psi = linearCoef(dtM)
			t1 = psi[0] * sM + psi[1] * sO
		
		t2 = np.dot(mat1, stokesVec)
		stokesVec = np.dot(mat2, t2) + t1
		
	return stokesVec

z = np.arange(100)
opacity = np.zeros((7,100))
emissivity = np.zeros((4,100))
opacity[0,:] = np.ones(100)*0.01
opacity[1,:] = np.zeros(100)
opacity[2,:] = np.zeros(100)
opacity[3,:] = np.ones(100) * 0.001
opacity[4,:] = np.zeros(100)
opacity[5,:] = np.zeros(100)
opacity[6,:] = np.zeros(100)
emissivity[0,:] = np.ones(100) * 0.001
emissivity[1,:] = np.zeros(100)
emissivity[2,:] = np.zeros(100)
emissivity[3,:] = np.zeros(100)
boundary = np.asarray([1.0,0,0,0])
res = deloparFormal(z, opacity, emissivity, boundary)