# acerta labels e corrige os numeros usados nos arquivos.
#
# uso: coloca os arquivos .res a serem convertidos no arquivo.
# lista-arquivos.in, e executa python acerta.py
# RR. nov. 2011.

from os import system

def acertanumb(strg):
	N = 3
	while len(strg) < N:
		strg = '0' + strg
	return strg

def extrainumero(strg,preamb):
	N = len(preamb)
	numb = ''
	count = 0
	for i in range(N,len(strg) - 4):
		numb = numb + strg[i]
		count += 1
	return numb 


arqs = file('lista-arquivos.in','r').readlines()

preamb = 'nvt-step-'

command = 'x=1; for i in *.xyz; do counter=$(printf  $x);'          # $(printf %03d $x) p/ acerto de ind.
command += 'ln -s "$i" ./'+ preamb+'-' +'"$counter".xyz; x=$(($x+1)); done'

for each in arqs:
	name   = each.split()[0]

	numb   = extrainumero(name,preamb)
	numbok = acertanumb(numb)
	
	print numb,'-->' , numbok
#	system('python shelx2xyz.py '+ name)
	system('mv '+ name+ ' '  + preamb + numbok + '.xyz')
#	system('mv '+ name+ '-MAT.in '+ preamb + numbok + '-MAT.in')
	
#system(command)



