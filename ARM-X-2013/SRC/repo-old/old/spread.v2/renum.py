#!/usr/bin/python

#recompondo um arquivo com manipulacao do formato
#e.g. corrigindo labels e retirando Hs
# >> reorganizadoNN.pdb
# v.2. - RR/ set.2010

from os import system
import sys

def acerta(string,N):
    while len(string) < N:
        string = '0' + string
    return string

def acertaV(string,N):
    while len(string) < N:
        string = ' ' + string
    return string

def corinthians(renum,linha):
  
    split = linha.split()
    output.write( 'HETATM' + ' ')
    output.write( acertaV(str(renum),4)+ ' ')


    for j in range(12,20):
        output.write( linha[j])

    output.write('          ')

    for j in range(30,len(linha)):
        output.write( linha[j])

    return 

print 'gerando reorganizado.pdb'
print 'CAUTION: maybe you need edit this prog.'
print '--------------------------------------'


#
# arq. info
info = []
for st in file('info', 'r').readlines():
    info.append(st.split()[0])
#
#
N = len(info[2])

for num in range(int(info[2])):

    print num + 1

    dir = info[0]+'/'+ acerta(str(num+1),N)+'/chains/'

    system('ls -1 '+ dir +'*.pdb  > dados')

    dados = file ('dados', 'r').readlines()

    output = file(dir+'/../../reorg'+ acerta(str(num+1),N) +'.pdb','w')
    renum = 1       # indice para renumeracao geral.

    for dado in dados:
    
        lines = file(dado.split()[0], 'r').readlines()

        lastinp = len(lines[0]) - 4
        par = []
        ct  = 0
        while ct < len(lines):

            if lines[ct][lastinp] == 'C':
                par.append(lines[ct])
        
            ct +=1
        
        for k in range(len(par)):    

            corinthians(renum,par[k])
            renum += 1


    output.write("END \n")
