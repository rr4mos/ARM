#!/usr/bin/python
# separa arquivos xyzs em arquivos independentes
# de tamanho dado (cadeias) no dir. chains. e chainsout (sem Hs)

#
# RR. set.2010 // generalizacao para lista de sistemas. v.2. // 
# RR. out.2011 // trabalhando com xyz por simplicidade: numeracao implicita.
#

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
    N = 3
    while len(string) < N:
        string = '0' + string
    return string

#
# arq. info
info = []
for i in file('info', 'r').readlines():
    info.append(i)

N = len(info[2])

for num in range(int(info[2].split()[0])):

    inputfile = info[0].split()[0]+'/'+ info[1].split()[0]+ acerta(str(num + 1),N) + '.xyz'    
    print 'input:', inputfile
    ncad      = [int(info[3].split()[0]),int(info[3].split()[1])]
    ltype     = [int(info[4].split()[0]),int(info[4].split()[1])]
    dirnew    = info[0].split()[0]+'/'+ acerta(str(num + 1),N)

    mkdir =  info[0].split()[0]+'/'+ acerta(str(num + 1),N)
    system('mkdir -p ' + mkdir)
    print 'mkdir:', mkdir, 'chains, e chainsout'

    dirout  = dirnew + '/chains' + '/'
    system('mkdir -p ' + dirout)

    dirout2  = dirnew + '/chainsout' + '/'
    system('mkdir -p ' + dirout2)
    
    lines = file(inputfile,'r').readlines()

    saida    = []
    catoms   = []
    cabxyz   = 2
    atm      = 0
    itype = 1
    ntype = len(ltype)
    total = 0
    atmtype = ['p3ht', 'c60'] 
    chainnum = 0

    print 'ntype:', ntype
    print 'ltype/ncad:',ltype, ncad
 
    for itype in range(ntype):
       print 'itype:',itype, ltype[itype],ncad[itype]
       print ltype[itype]*ncad[itype]
       while atm < (total+ltype[itype]*ncad[itype]):
               saida.append( lines[atm+cabxyz] ) 

               if saida[len(saida)-1].split()[0] != 'H':
                   catoms.append(atm)
               atm += 1
  
               if ((atm-total)  % ncad[itype] ) == 0:
                   print 'okatm:',atmtype[itype], atm
                   chainnum +=1
                   name  = dirout + acerta(str(chainnum),3) + '.xyz'
                   name2 = dirout2 + acerta(str(chainnum),3) + '.xyz'

# escrevendo o chains
                   fout = file(name, 'w')
                   fout.write(str(ncad[itype])+'\n\n')
                   for any in saida:
                       fout.write(any)
                   fout.close()

# escrevendo o chainsout
                   fout = file(name2, 'w')
                   fout.write(str(len(catoms))+'\n\n')
                   for any in catoms:
                        fout.write(lines[any+cabxyz])
                   fout.close()
                   saida  = []
                   catoms = []
       total += atm       
       print 'geradas:', chainnum, 'cadeias'
