#!/usr/bin/python

# 1. ajusta nomes dos arqs. (no. de chars do N=2, em acerta)

from os import system

def acerta(string,N):
    while len(string) < N:
        string = '0' + string
    return string

#
# arq. info
info = []
for st in file('info', 'r').readlines():
    info.append(st.split()[0])

N = len(info[2])
#
#

system('rm -rf '+ info[0])
system('cp -r analisa ' + info[0])

system( 'mkdir '+ info[0] + '/novos/')

for num in range(int(info[2])):
   
    nameWR = info[1]+ str(num + 1)         + '.pdb'
    nameOK = info[1]+ acerta(str(num + 1),N) + '.pdb'
    system('wc -l ' + info[0] + '/' + nameWR + ' > .wcl') 
    wcl = int(file('.wcl', 'r').readlines()[0].split()[0])
    system('tail -n ' + str(wcl-5) + ' '\
                      + info[0]+'/'+ nameWR+ '> '\
                      + info[0]+'/novos/'+nameOK)
                   

system('rm -f ' + info[0] + '/*.pdb ' )
system('mv '    + info[0] + '/novos/* ' + info[0])
system('rm -rf '+ info[0] + '/novos')

