#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import asciitable

rstar = 4.0e13
data = asciitable.read("RES_MOL1_d.DAT",header_start=None,data_start=0,delimiter=" ",comment="!",guess=False)
r_lte  = data['col3']
cn_lte = np.power(10.0e0,data['col46']) * 2.0e0
data = asciitable.read("irc10216_smo_out_01.dat",header_start=None,data_start=0,delimiter=" ",comment="!",guess=False)
r_kin  = np.power(10.0e0,data['col2']) / rstar
cn_kin = np.power(10.0e0,data['col87']) * 2.0e0

rmin = 1.0
rmax = 10000.0
amin = 1.0e-10
amax = 1.0e-3
plt.figure(1)
plt.subplot(111)
plt.axis([rmin,rmax,amin,amax])
plt.title('CN abundance in IRC+10216')
plt.xlabel("radius (R*)")
plt.ylabel("abundance relative to H2")
plt.loglog(r_lte,cn_lte,'m-',label="CN LTE")
plt.loglog(r_kin,cn_kin,'b-',label="CN kinetics")
plt.legend(loc="upper left")
plt.savefig("cn.png")
plt.show()
