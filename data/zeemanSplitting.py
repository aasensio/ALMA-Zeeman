import numpy as np
def zeemanSplitting(S, I, N, J, F):
	"""
	Return the Zeeman splitting factor in Hz/microG
	"""
	alphaJ = (F*(F+1) + J*(J+1) - I*(I+1)) / (2*F*(F+1))
	alphaI = (F*(F+1) + I*(I+1) - J*(J+1)) / (2*F*(F+1))
	gJ = (J*(J+1)+S*(S+1)-N*(N+1)) / (J*(J+1))
	gI = 0.4036#*0.0
		
	return -2*(gJ * alphaJ + gI * alphaI) * 1.3996
	
Nu = [1,1,1,1,1,1,1]
Nl = [0,0,0,0,0,0,0]
Ju = [0.5,0.5,0.5,1.5,1.5,1.5,1.5]
Jl = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
Fu = [0.5,1.5,1.5,1.5,2.5,0.5,1.5]
Fl = [1.5,0.5,1.5,0.5,1.5,0.5,1.5]

for i in range(7):
	maxim = 0.0
	print 'hey'
	sl = zeemanSplitting(0.5, 1.0, Nl[i], Jl[i], Fl[i])
	su = zeemanSplitting(0.5, 1.0, Nu[i], Ju[i], Fu[i])
	Mu = np.arange(-Fu[i],Fu[i]+1)
	Ml = np.arange(-Fl[i],Fl[i]+1)
	for u in Mu:
		for l in Ml:
			if (np.abs(u-l) == 1):
				print l, u, sl*l-su*u
print "{0} Hz/microG".format((su-sl))

#sl = zeemanSplitting(0.5, 1.0, 0.0, 0.5, 0.5)
#su = zeemanSplitting(0.5, 1.0, 1.0, 0.5, 1.5)
#print "{0} Hz/microG".format(2*(su-sl))

#sl = zeemanSplitting(0.5, 1.0, 0.0, 0.5, 1.5)
#su = zeemanSplitting(0.5, 1.0, 1.0, 0.5, 1.5)
#print "{0} Hz/microG".format(2*(su-sl))

#sl = zeemanSplitting(0.5, 1.0, 0.0, 0.5, 1.5)
#su = zeemanSplitting(0.5, 1.0, 1.0, 1.5, 2.5)
#print "{0} Hz/microG".format(2*(su-sl))
