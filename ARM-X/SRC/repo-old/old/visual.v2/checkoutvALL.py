#!/usr/bin/python

from visual import *
import sys
import random


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

scene.center=(a1+a2+a3)*0.5 # camera para o centro da caixa.

if len(sys.argv) < 2:
    print 'falta argumento... saindo :('
    print 'chka [atoms/map/all]'
    exit()

list = file('chainsplot.in', 'r').readlines()

N = len(list[0].split()[0]) - 4

vai = sys.argv[1]

colorshift = [-1,0,1]
shf = 0
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

        fct1 = 1.0 + colorshift[0]*shf
        fct2 = 0.0 + colorshift[1]*shf
        fct3 = 0.0 + colorshift[2]*shf

        fcf1 = 1.0 + colorshift[0]*shf
        fcf2 = 1.0 + colorshift[1]*shf
        fcf3 = 0.0 + colorshift[2]*shf

        Ndeg = len(file(nchain.split()[0],'r').readlines())
        dfct1 = (fcf1-fct1)/Ndeg
        dfct2 = (fcf2-fct2)/Ndeg
        dfct3 = (fcf3-fct3)/Ndeg

        ball=[]
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

            fctnext1 = fct1+dfct1 #random.random()
            fctnext2 = fct2+dfct2 #random.random()
            fctnext3 = fct3+dfct3 #random.random()
    
            color1 = vector(fct1,    fct2,    fct3    ) 
            color2 = vector(fctnext1,fctnext2,fctnext3)

            print color1
            print color2
            print '---' 
            

            ball.append( sphere(pos=center, radius=1.0, color=color1))
            pointerN  = arrow(pos=center, axis=normal, color=(255,255,255))
            pointerNN = arrow(pos=center, axis=vnn,color=color2)

            print color1

            fct1 = fctnext1
            fct2 = fctnext2
            fct3 = fctnext3

        arq.close()
        shf+=1.0
