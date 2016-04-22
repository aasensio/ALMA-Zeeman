import numpy as np
import pydelopar

z = np.arange(100, dtype=np.float64)
opacity = np.zeros((100,7), dtype=np.float64)
source = np.zeros(100, dtype=np.float64)
opacity[:,0] = np.ones(100)*0.01
opacity[:,1] = np.zeros(100)
opacity[:,2] = np.zeros(100)
opacity[:,3] = np.ones(100) * 0.001
opacity[:,4] = np.zeros(100)
opacity[:,5] = np.zeros(100)
opacity[:,6] = np.zeros(100)
boundary = np.asarray([1.0,0,0,0], dtype=np.float64)
res = pydelopar.deloparFormal(z, opacity, source, boundary)
print res