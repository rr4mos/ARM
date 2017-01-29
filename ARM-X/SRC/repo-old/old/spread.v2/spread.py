#!/usr/bin/python
# separa arquivos pdbs em arquivos independentes
# de tamanho dado no dir. chains.

#
# generalizacao para lista de sistemas. v.2. // RR. set.2010


import sys
from os import system

def char6(string):
    hetatm = ''
    if len(string) >= 6:
        for i in range(6):
            hetatm+=string[i]
    else:
        hetatm = string
    return hetatm

def acerta(string,N):
    while len(string) < N:
        string = '0' + string
    return string

#
# arq. info
info = []
for st in file(sys.argv[1], 'r').readlines():
    info.append(st.split()[0])
#
#
N = len(info[2])



if len(sys.argv) < 3:
    print 'falta um argumento, parando...'
    exit()


for num in range(int(info[2])):

    inputfile = info[0]+'/'+ info[1]+ acerta(str(num + 1),N) + '.pdb'    
    print 'input:', inputfile
    ncad      = int(info[3])    
    dirnew    = info[0]+'/'+ acerta(str(num + 1),N)

    mkdir =  info[0]+'/'+ acerta(str(num + 1),N)
    system('mkdir -p ' + mkdir)
    print 'mkdir:', mkdir, 'e ' + sys.argv[2]

    dirout  = dirnew + '/' + sys.argv[2] + '/'
    system('mkdir -p ' + dirout)
    
    lines = file(inputfile,'r').readlines()

    saida    = []
    atm      = 0
    chainnum = 0 

    ct = 0
    while ct < (len(lines)):

        if char6(lines[ct].split()[0]) == 'HETATM':
            saida.append( lines[ct] )
            atm += 1
            if (atm  % ncad ) == 0:

                chainnum +=1
                name = dirout + acerta(str(chainnum),3) + '.pdb'
                fout = file(name, 'w')

                for any in saida: 
                    fout.write(any)

                fout.close()
                saida = []
        
        ct += 1
    print 'geradas:', chainnum, 'cadeias'

