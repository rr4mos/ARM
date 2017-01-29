#!/usr/bin/python

#gera imagem dos vetores de um arquivo em grupos de cores


from visual import *
import sys
from sys import exit

colors = [(255,000,000),\
          (000,255,000),\
          (000,000,255),\
          (255,255,000),\
          (255,000,255)]


arq = 'fort.666'
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

    if ((i+1)%len(colors)==0):
        center = center+delta

    i+=1
