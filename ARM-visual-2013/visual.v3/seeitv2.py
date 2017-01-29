#!/usr/bin/python

#gera imagem dos vetores de um arquivo em grupos de cores


from visual import *
import sys
from sys import exit

colors = [(255,000,000),\
          (000,255,000),\
          (000,000,255),\
          (255,255,000),\
          #(255,000,255),\
          (255,255,255)]

arq = 'intra.vec.ang-rg.dat'
vecs = file(arq,'r').readlines()

pointer = []

center = vector( 0.0,0.0,0.0)
delta  = vector( 1.0,0.0,0.0)
i = 0
for each in vecs:

    axisv = vector( float(each.split()[0]),\
                    float(each.split()[1]),\
                    float(each.split()[2]))
                           
    pointer.append(arrow(pos=center, axis=axisv, color=colors[i%len(colors)]))
    
    if ((i-1)%len(colors))== 0:
        eixo = axisv

    if ((i+1)%len(colors)==0):
        center = center + 1.5*eixo

    i+=1
