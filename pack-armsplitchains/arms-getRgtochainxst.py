import sys 

N = 40
for ch in range(N):
	alldir = file('execall.in','r').readlines()

	for diret in alldir:
		print file(diret.split()[0]+'/Rg2.dat','r').readlines()[ch].split()[1]
	print '#\n'
