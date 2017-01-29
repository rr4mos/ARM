#!/usr/bin/python

def acerta(string,N):
    while len(string) < N:
        string = '0' + string
    return string

import sys
from os import system

system("awk '{print $3, $4, $5}' intra.vec.1nb.dat > .aux")

system('paste all.vec.rg-pbc.dat .aux > x; mv x centro.dat')

arq = file('centro.dat', 'r').readlines()

system('mkdir -p centros')

chain = acerta(arq[0].split()[0],2)
out = file('./centros/' + chain + '.dat','w')

for line in arq:
    ccur = acerta(line.split()[0],2)

    if ccur != chain:        
        out = file('./centros/' + ccur + '.dat','w')    
        chain = ccur

    out.write(line)
