import numpy as np

fIn = open('SO_V1.FRE', 'r')
fOut = open('sov1.linelist', 'w')
fOut.write(' LINE   Nu  Ju  Nl  Jl     FREQ(Hz)     Error     Eupp     ul      Sij       gu\n')
fOut.write('-----  --- --- --- --- -------------- ---------- ------- --------- ---------- ---\n')

for i in range(706):
    line = fIn.readline()
    dat = line.split()

    Nu = int(dat[3])
    Ju = int(dat[4])
    Nl = int(dat[5])
    Jl = int(dat[6])
    freq = float(dat[1])
    error = float(dat[2])
    Eup = float(dat[9])
    Aul = float(dat[10].split('=')[1])
    Sij = float(dat[7].split('=')[1][:-1])
    gu = 2.0 * Ju+1
    
    fOut.write("{0}   {1}   {2}   {3}   {4}   {5}   {6}   {7}    {8}   {9}    {10}\n".format(i,Nu,Ju,Nl,Jl,freq,error,Eup,Aul,Sij,gu))
    print(i,Nu,Ju,Nl,Jl,freq,error,Eup,Aul,Sij,gu)

fOut.close()
fIn.close()