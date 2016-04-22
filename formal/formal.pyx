from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray

cdef extern from "pyformal.h":
	void c_delopar_source(int* n, double* height, double* opacity, double* source, double* boundary, double* stokesOut)
	
	void c_delopar(int* n, double* height, double* opacity, double* emissivity, double* boundary, double* stokesOut)
	

def deloparFormalSource(ar[double] height, ar[double, ndim=2] opacity, ar[double] source, ar[double] boundary):
	"""
	stokes = deloparFormalSource(z, opacity, source, boundary)
	
	--- Input ---
	z : vector of length n. The RT equation is solved from the fist point to the last one
	opacity : array of size [n,7] containing etaI, etaQ, etaU, etaV, rhoQ, rhoU, rhoV
	source : vector of length n with the source function at each point
	boundary : vector of length 4 with the boundary condition for the four Stokes parameters
	
	--- Output ---
	stokes : vector of length 4 with the emergent Stokes parameters
	"""
	
	cdef:
		int n = height.size		
		ar[double] stokesOut = empty(4, order='F')
		ar[double, ndim=2, mode="c"] opa
   
# Make sure that the 2D array is C_CONTIGUOUS
	opa = ascontiguousarray(opacity)
	c_delopar_source(&n, &height[0], &opa[0,0], &source[0], &boundary[0], <double*> stokesOut.data)
	
	return stokesOut
	
	
def deloparFormal(ar[double] height, ar[double, ndim=2] opacity, ar[double, ndim=2] emissivity, ar[double] boundary):
	"""
	stokes = deloparFormal(z, opacity, emissivity, boundary)
	
	--- Input ---
	z : vector of length n. The RT equation is solved from the fist point to the last one
	opacity : array of size [n,7] containing etaI, etaQ, etaU, etaV, rhoQ, rhoU, rhoV
	emissivity: array of size [n,4] containing epsI, epsQ, epsU, epsV at each point
	boundary : vector of length 4 with the boundary condition for the four Stokes parameters
	
	--- Output ---
	stokes : vector of length 4 with the emergent Stokes parameters
	"""
	
	cdef:
		int n = height.size		
		ar[double] stokesOut = empty(4, order='F')
		ar[double, ndim=2, mode="c"] opa
		ar[double, ndim=2, mode="c"] emi
   
# Make sure that the 2D array is C_CONTIGUOUS
	opa = ascontiguousarray(opacity)
	emi = ascontiguousarray(emissivity)
	c_delopar(&n, &height[0], &opa[0,0], &emi[0,0], &boundary[0], <double*> stokesOut.data)
	
	return stokesOut