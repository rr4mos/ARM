#!/usr/bin/python

from visual import *
import sys

# colocando a celula unitaria

mat = file('../MAT.in','r').readlines()

a1 = vector(float(mat[0].split()[0]),\
            float(mat[0].split()[1]),\
            float(mat[0].split()[2]))

a2 = vector(float(mat[1].split()[0]),\
            float(mat[1].split()[1]),\
            float(mat[1].split()[2]))

a3 = vector(float(mat[2].split()[0]),\
            float(mat[2].split()[1]),\
            float(mat[2].split()[2]))

aresta = 0.2

box1 = cylinder(pos=vector(0,0,0), axis=a1, radius=aresta, color=(255,0,0))
box2 = cylinder(pos=vector(0,0,0), axis=a2, radius=aresta, color=(0,255,0))
box3 = cylinder(pos=vector(0,0,0), axis=a3, radius=aresta, color=(0,0,255))
box4 = cylinder(pos=a1, axis=a2, radius=aresta)
box5 = cylinder(pos=a2, axis=a1, radius=aresta)
box6 = cylinder(pos=a1, axis=a3, radius=aresta)
box7 = cylinder(pos=a3, axis=a1, radius=aresta)
box8 = cylinder(pos=a2, axis=a3, radius=aresta)
box9 = cylinder(pos=a3, axis=a2, radius=aresta)
box10 = cylinder(pos=(a1+a3), axis=a2, radius=aresta)
box11 = cylinder(pos=(a1+a2), axis=a3, radius=aresta)
box12 = cylinder(pos=(a2+a3), axis=a1, radius=aresta)

if len(sys.argv) < 2:
    print 'falta argumento... saindo :('
    print 'chka [atoms/map/all]'
    exit()

list = file('chainsplot.in', 'r').readlines()

N = len(list[0].split()[0]) - 4

vai = sys.argv[1]

for nchain in list:
        
#colocando os atomos do arquivo pdb.

    if vai=='atoms' or vai=='all':

        Cnchain = ''
        for i in range(N):
            Cnchain = Cnchain + nchain[i]

        arq = file('../chainsout/'+ Cnchain +'.pdb','r') 

        for line in arq:
            atom =   vector( float(line.split()[4]),\
                         float(line.split()[5]),\
                         float(line.split()[6]) )

            ball = sphere(pos=atom, radius=0.6, color=(1,1,1))

        arq.close()

#colocando elementos de mapeamento.
 
    if vai=='map' or vai=='all':
        arq = file(nchain.split()[0],'r') 
        for line in arq:
            center = vector( float(line.split()[2]),\
                         float(line.split()[3]),\
                         float(line.split()[4]) )

            normal = 4*vector( float(line.split()[5]),\
                           float(line.split()[6]),\
                           float(line.split()[7]) )

            vnn    = vector( float(line.split()[8]),\
                         float(line.split()[9]),\
                         float(line.split()[10]) )

            ball = sphere(pos=center, radius=1.0, color=(1,2,2))
            pointerN  = arrow(pos=center, axis=normal)
            pointerNN = arrow(pos=center, axis=vnn,color=(255,0,0))

        arq.close()

