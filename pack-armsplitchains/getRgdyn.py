N = 40
for ch in range(N):
	alldir = file('paths.in','r').readlines()
        for diret in alldir:
		if len(diret.split())== 0:
			break
		name = diret.split()[0]+'/Rgm.dat'
                
		print file(name,'r').readlines()[ch].split()[1]
        print '#\n'

