#!/usr/bin/python

# script completamente arbitrario por que apenas
# ajusta nomes dos arqs. especificidade dos arquivos que geramos nos tcl que 
# nao sao da forma 00001,00002,10000 mas 1,2, 10000
#
from os import system
#
def acerta(string,N):
    N = 3
    while len(string) < N:
        string = '0' + string
    return string
#
info = []
for st in file('info', 'r').readlines():
    info.append(st.split()[0])
#
N = len(info[2])
#
system('rm -rf '+ info[0])
system('mkdir '+ info[0])
#system('cp -r analisa/*.xyz ' + info[0])
#
for num in range(int(info[2])):
#   
    nameWR = info[1]+ str(num + 1)         + '.xyz'
    nameOK = info[1]+ acerta(str(num + 1),N) + '.xyz'
#
#    system('cp '+ 'analisa/'+ nameWR + ' '\
#                + info[0]+'/'+nameOK)

    system('cp '+ 'analisa/'+ nameOK + ' '\
                + info[0]+'/'+nameOK)


