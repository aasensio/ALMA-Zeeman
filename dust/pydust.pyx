cdef extern from "pydust.h":
	void kappa_dust(double *freq, int *iddust, double *adust, double *kdust)
			
def kappaDust(double freq, int iddust, double adust):
	cdef:
		double kdust
	kappa_dust(&freq, &iddust, &adust, &kdust)
    
	return kdust