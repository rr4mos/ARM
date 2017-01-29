#!/usr/bin/python

#falta: 1) cortar soh no sitio (e nao cadeia como atual)
#       2) colorir diferente cada sitio
#       3) colorir vetor intersitio diferente dos sitios e dos vetores normais


# visualiza o sitio "sitio" incluindo sitios conectados e setas de conexao
# RR. dez/2011
# uso: seeit [-pbc] "sitio"
#      no diretorio onde rodou o ARM (stst)

import sys
from os import system
from visual import *

# variaveis externas implicitas (tamanho da cadeia)

#Ncad = 40 # no. de cadeias
#Nanl = 26 # no. de aneis

Ncad = 96 # no. de cadeias
Nanl = 5 # no. de aneis

# variaveis externas explicitas: sitio

def acerta(string,N):
	# acerta nome da string, completando com 0's ateh N.
    while len(string) < N:
        string = '0' + string
    return string

def tiraredund(vec):
	# elimina termos redundantes do vetor vec
	vecok = [vec[0]]
	for every in vec:
		flag = 1
		for i in range(len(vecok)):
			if every == vecok[i]:
				flag = 0
		if flag == 1:
			vecok.append(every)
	return vecok

def aloca(Ncad,Nanl):
	# aloca matriz Ncad x Nanl 
	centros = []
	for i in range(Ncad):
		centros.append([])
		for j in range(Nanl):
			centros[i].append([])
	return centros
	
def caixa(): # colocando a celula unitaria

	mat = file('./MAT.in','r').readlines()

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



def visualizacadeia(nchain,cor,imprimeatomos):

		d1 = -vector(255,0,0)
		d2 = -vector(0,255,0)
		d3 =  vector(0,0,0)

		cor1 = cor + d1
		cor2 = cor + d2
		cor3 = cor + d3
		if imprimeatomos == 1:
			arq = file('./chainsout/'+ nchain +'.xyz','r') 

		# alterado para tratar xyz
			for line in arq:
				if len(line.split()) == 4:
					atom =   vector( float(line.split()[1]),\
							         float(line.split()[2]),\
								     float(line.split()[3]) )
    
					ball = sphere(pos=atom, radius=0.5, color=(0.5,0.5,0.5))
    
			arq.close()

		#colocando elementos de mapeamento

		arq = file('./centros/'+nchain +'.dat','r') 
		Ndeg = len(file('./centros/'+nchain +'.dat','r').readlines())

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


			ball.append( sphere(pos=center, radius=0.5,  color=cor1))
			pointerN  =   arrow(pos=center, axis=normal, color=cor2, shaftwidth= 0.1)
			pointerNN =   arrow(pos=center, axis=vnn,    color=cor3, shaftwidth= 0.5)

#
## programa principal
#

if len(sys.argv) == 3: 
	arquivo1  = "all.vec.rg"+sys.argv[1]+".dat" # melhor usar com pbc
	sitio     =  int(sys.argv[2])
	system('brk -pbc')
else:
	if len(sys.argv) == 2:
		arquivo1  = "all.vec.rg"+".dat" # melhor usar com pbc
#		print 'se quiser com pbc use: brk -pbc'
		sitio     =  int(sys.argv[1])
		system('brk')
	else:
		print 'problemas com o numero de entradas do programa'
		exit()
	
arquivo2  = "intra.conj"    +".dat" 
arquivo3  = "inter.site.vec"+".dat"

arqaneis  = file(arquivo1, 'r').readlines()
asitemap  = file(arquivo2, 'r').readlines()

centrosx = aloca(Ncad,Nanl)
centrosy = aloca(Ncad,Nanl)
centrosz = aloca(Ncad,Nanl)

sitemap = aloca(Ncad,Nanl)
invsitemap = aloca(Ncad*Nanl,0)

corcon = vector(255,0,0)

for i in range(Ncad*Nanl):

	# pega coordenadas dos centros
	chainID = int(arqaneis[i].split()[0])-1
	ringID  = int(arqaneis[i].split()[1])-1

	centrosx[chainID][ringID] = float(arqaneis[i].split()[2])	
	centrosy[chainID][ringID] = float(arqaneis[i].split()[3])
	centrosz[chainID][ringID] = float(arqaneis[i].split()[4])
#	print chainID,ringID, centrosx[chainID][ringID],centrosy[chainID][ringID],centrosz[chainID][ringID]
	# pega mapeamento centros <-> sitios
	
	sitemap[chainID][ringID] = int(asitemap[i].split()[2])
	invsitemap[int(asitemap[i].split()[2])].append([chainID,ringID])
#	print chainID,ringID, sitemap[chainID][ringID]

# reduzindo consumo de memoria
arqaneis  = [] 
asitemap  = []

# le conexoes

arqconex  = file(arquivo3, 'r').readlines()

centrosconex = aloca(Ncad,Nanl)

sitant = 0
caixa()
chligada = []
for line in arqconex:
		
	chi = int(line.split()[0])-1
	rgi = int(line.split()[1])-1
	sti = int(line.split()[2])
	
	sitat = sti
	
	if sti == sitio:
		chj = int(line.split()[3])-1
		rgj = int(line.split()[4])-1
	
		print '(', chi,',',rgi,') '+ '---> ', '(',chj,',', rgj,')'
	
		chligada.append(str(chj+1))
		cadeia_sitio = str(chi+1)

		centrosconex[chi][rgi].append( [chj,rgj] )

		vx = float(line.split()[6])
		vy = float(line.split()[7])
		vz = float(line.split()[8])

		conexion = vector((vx,vy,vz))
		centeri  = vector((centrosx[chi][rgi],centrosy[chi][rgi],centrosz[chi][rgi]))
		centerj  = vector((centrosx[chj][rgj],centrosy[chj][rgj],centrosz[chj][rgj]))

		ball     = sphere(pos=centeri, radius=0.6   , color=corcon)
		ball     = sphere(pos=centerj, radius=0.6   , color=corcon)
		pointerN =  arrow(pos=centeri, axis=conexion, color=corcon)

chligada = tiraredund(chligada)

print sitio, 'cadeia:', cadeia_sitio, ' --- ligado a:', chligada

cor = vector(1,1,1)
imprimeatomos = 1

visualizacadeia(acerta(str(cadeia_sitio),3),cor,imprimeatomos)
#print chligada, tiraredund(chligada)

ligadas = []
if len(chligada) != 0:
	for eachligada in chligada:
		ligadas.append(eachligada)
		visualizacadeia(acerta(eachligada,3),cor,imprimeatomos)

#ETAPA DE VISUALIZACAO DE SITIOS

#1. quais sao os aneis de um dado sitio? vide invsitemap[sitio] (construido anteriormente)

cadeiaprint = []
sitei = []
for each in invsitemap[sitio]:
	sitei.append( each )

cadeiaprint.append(sitei)
#print 'sitios-i:', cadeiaprint
#print invsitemap[sitio]

#2. adicionando sitios ligados

print centrosconex[6][0]

lista = []
for each in invsitemap[sitio]:
	print each

	for every in centrosconex[each[0]][each[1]]:
		print every
		lista.append(sitemap[every[0]][every[1]])

lista = tiraredund(lista)

for every in lista:
	pega = invsitemap[every]
	cadeiaprint.append( pega )

# 3. mapa de todos os atomos <--> centros 
# cabeca: 6 aneis de C
# unidade: 6 aneis de C (benzeno) + 2 aneis de C (vinil)

atomshead    = 6 
atomsperunit = 8
	
indicesporanel = [[0,atomshead-1]]
ind = (atomshead - 1) + 1
for i in range(Nanl):
	indicesporanel.append([ind,ind+atomsperunit-1])
	ind = ind + atomsperunit 

# imprime todos os aneis do vetor cadeiaprint
first = 1
for siteprint in cadeiaprint:
	
	corapp = cor - vector(random.random(),\
						  random.random(),\
						  random.random())
	if corapp == vector(0,0,0):
		corapp = vector(1,1,1)
	if first == 1:
		corapp = vector(1,1,1)
		first = 0

	for center in siteprint:

		nchainarq = acerta(str(center[0]+1),3)
		cadarq    = './chainsoutpbc/'+ nchainarq +'.xyz'
		print cadarq
		arq = file(cadarq,'r').readlines()

		atom = []
		for line in arq:
			if len(line.split()) == 4:
				atom.append( vector( float(line.split()[1]),\
							         float(line.split()[2]),\
								     float(line.split()[3]) ))
		i = center[1]
		for j in range(indicesporanel[i][0],indicesporanel[i][1]+1):
			ball = sphere(pos=atom[j], radius=0.6, color=corapp)
	
	
