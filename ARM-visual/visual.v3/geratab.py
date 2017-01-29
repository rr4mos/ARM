#!/usr/bin/python
from os import system

# lendo comp. de conjug.

conjug = []
arqc = file('../intra.conj.dat','r')

csant = 'x'
for each in arqc:
	cs = each.split()[2]
	cl = each.split()[3]
	if cs != csant:
		conjug.append(cl)
		#print cl, len(conjug)
	csant = cs	

arqc.close()

diretin = './detalhe/'
inputs=file('probflux.in','r').readlines()

diretout='./ok-detalhe/'
system('mkdir '+diretout)

for each in inputs:
	arq=each.split()[0]
	arqin =file(diretin +arq,'r').readlines()

	arqout=file(diretout+arq,'w')
	i = 0
	for data in arqin:

		data = data.split()
		datap = ''
		for k in range(1,len(data)):
			datap+= ' '+ data[k]

		arqout.write(data[0]+' '+conjug[i]+ ' '+datap+'\n')
		i+=1
	print 'gerado ==> '+ diretout+arq
