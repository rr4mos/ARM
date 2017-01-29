#!/usr/bin/python

from math import sqrt
from sys import exit

#
#gera graphml incluindo as cargas e tamanho do sitio como propriedades do grafo.
#

# lendo cargas
charges = []
arqq = file('charges.dat', 'r')
for each in arqq:
	if len(each) >= 2:
		charges.append(each.split()[1])
arqq.close()

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

print '<?xml version="1.0" encoding="UTF-8"?>'
print '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">'
print '<key attr.name="label" attr.type="string" for="node" id="label"/>'
print '<key attr.name="Edge Label" attr.type="string" for="edge" id="edgelabel"/>'
print '<key attr.name="weight" attr.type="double" for="edge" id="weight"/>'
print '<key attr.name="Edge Id" attr.type="string" for="edge" id="edgeid"/>'
#print '<key attr.name="color" attr.type="integer" for="node" id="color"/>'
print '<key attr.name="x" attr.type="float" for="node" id="x"/>'
print '<key attr.name="y" attr.type="float" for="node" id="y"/>'
print '<key attr.name="size" attr.type="float" for="node" id="size"/>'
print '<key attr.name="local charge" attr.type="double" for="node" id="charge"/>'
print '<key attr.name="conjug. length" attr.type="double" for="node" id="conjug"/>'
print '<graph edgedefault="undirected">'

mobilines = file('mobilimax.input','r').readlines()

lastnode = int(mobilines[len(mobilines)-1].split()[0])

Nsqrt = int(sqrt(lastnode)) 
Nsqrt = Nsqrt + lastnode/Nsqrt
idy = 1

for id in range(lastnode):
        idx= id +1 - ((id)/Nsqrt )*Nsqrt
        
        print '<node id="'+str(id+1)+'">'
        print '<data key="label">'+str(id+1)+'</data>'
        print '<data key="size">10.0</data>'
        print '<data key="x">'+str(idx)+'</data>'
        print '<data key="y">'+str(idy)+'</data>'
        print '<data key="charge">'+str(charges[id])+'</data>'
        print '<data key="conjug">'+str(conjug[id])+'</data>'
        print '</node>'
        if (idx%Nsqrt == 0):
            idy += 1
        id+=1

id =0
stan = 999999
enan = 999999
for line in mobilines:
    ok = 0
    st = int(line.split()[0])
    en = int(line.split()[1])
    if st < en:
#        print en, enan, st, stan
        if ((en == enan)and(st != stan)):
            ok = 1
        if (en != enan):
            ok = 1
        
    if ok == 1:
        print '<edge source="'+str(st)+'" target="'+ str(en)+'">'
        print '<data key="weight">1.0</data>'
        print '</edge>'
    enan = en
    stan = st

print '</graph>'
print '</graphml>'
