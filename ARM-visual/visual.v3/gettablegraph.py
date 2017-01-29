#!/usr/bin/python

#pega dados do arquivo graphml, gerando tabela geral.

def acerta(strg):
	N = 10
	while len(strg) < N:
		strg += ' '

	return strg

import sys
from os import system

arqfile = sys.argv[1]

parameters = ['label','degree','charge','conjug']

count = 1
for parameter in parameters:
	system('grep \'key=\"'+ parameter + "\"\' "+ arqfile + " > .data" + str(count))
	count += 1 

geraldados = []
for i in range(count-1):
	arqfile = '.data' + str(i+1)
	arq = file(arqfile,'r')
	
	dados = []	
	for line in arq:
		for j in range(len(line)):			
			if line[j] == '>':
				char = ''
				k=j+1
				while line[k] != '<':
					char += line[k]
					k+=1
				dados.append( char )
				break
		ndados=len(dados)
	geraldados.append(dados)

# imprime

output=file('tabelageral.dat','w')

for j in range(ndados):
	for i in range(count-1):	
		if i < count-2:
			output.write(acerta(geraldados[i][j])+'  ')		
		else:
			output.write(acerta(geraldados[i][j])+'\n')
		


