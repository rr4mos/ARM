#!/usr/bin/python

def acerta(string,N):
    while len(string) < N:
        string = '0' + string
    return string

import sys
from os import system

system("awk '{print $3, $4, $5}' intra.vec.1nb.dat > .aux")

N = 3

if len(sys.argv) > 1: 
    system("awk '{print $3, $4, $5, $6, $7, $8, $9, $10, $11}' all.vec.vinil"+sys.argv[1]+".dat > .aux2")
    system('paste all.vec.rg'+sys.argv[1]+'.dat .aux .aux2 > x; mv x centro.dat')    
else:
    system("awk '{print $3, $4, $5, $6, $7, $8, $9, $10, $11}' all.vec.vinil.dat > .aux2")
    system('paste all.vec.rg.dat .aux .aux2 > x; mv x centro.dat')
    print 'se quiser com pbc use: brk -pbc'

system('rm -f .aux*')

arq = file('centro.dat', 'r').readlines()

system('mkdir -p centros')

chain = acerta(arq[0].split()[0],N)
out = file('./centros/' + chain + '.dat','w')

for line in arq:
    ccur = acerta(line.split()[0],N)

    if ccur != chain:        
        out = file('./centros/' + ccur + '.dat','w')    
        chain = ccur

    out.write(line)
