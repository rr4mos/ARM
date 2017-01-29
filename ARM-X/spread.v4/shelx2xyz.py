#
# shelx2xyz.py
# converte arquivo de posicao em coordenadas reduzidas em termos dos vet. de rede para coord. cartesianas.
# RR. oct. 2011
#

from numpy import *
from math import sqrt
from math import pi
from math import cos
from math import sin
from sys import exit
import sys

arq = sys.argv[1]

data = file(arq,'r').readlines()
cabecalho = 3 

a = float(data[1].split()[2])
b = float(data[1].split()[3])
c = float(data[1].split()[4])
al = float(data[1].split()[5])*pi/180.0
bt = float(data[1].split()[6])*pi/180.0
gm = float(data[1].split()[7])*pi/180.0

#print  a,b,c,al,bt,gm

latta = array([ a, 0, 0 ])
lattb = array([ b*cos(gm), b*sin(gm), 0 ])

V = a*b*c*sqrt(1 - cos(al)**2 - cos(bt)**2 - cos(gm)**2 + 2.0*b*c*cos(al)*cos(bt)*cos(gm))
lattc = array([ c*cos(bt), c*(cos(al)-cos(bt)*cos(gm))/sin(gm), V/(a*b*sin(gm)) ])

M = array([ latta,lattb,lattc ])

matinout = file(arq+'-MAT.in', 'w')

matinout.write( '%5.10f '%M[0][0] + '%5.10f '%M[0][1] + '%5.10f'%M[0][1] + '\n' +\
                '%5.10f '%M[1][0] + '%5.10f '%M[1][1] + '%5.10f'%M[1][2] + '\n' +\
                '%5.10f '%M[2][0] + '%5.10f '%M[2][1] + '%5.10f'%M[2][2] + '\n')
matinout.close()

N = len(data) - cabecalho -1

atomsout = file(arq+'.xyz','w')

atomsout.write(str(N)+'\n\n')

for i in range(N):
	j = i + cabecalho
	type = data[j].split()[0][0]

	ata = float(data[j].split()[2])
	atb = float(data[j].split()[3])
	atc = float(data[j].split()[4])

	atfct = array([ ata, atb, atc ])
	atxyz = dot(atfct,M)

#	print 'atfct:', atfct
#	print 'atxyz:',atxyz

	atomsout.write( type + ' %.10f ' %atxyz[0]+ ' %.10f ' %atxyz[1] + ' %.10f\n'%atxyz[2])

atomsout.close()

