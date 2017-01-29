#!/usr/bin/python
# ajusta taxas de hopping pela media para rodar o mobili

from os import system

system('sort -n -k1 -k2 mobili.input > mobili.input-sort')
system('mv mobili.input-sort mobili.input')

mobili   = file('mobili.input','r').readlines()
mobiliok = file('mobilimax.input','w')

nlin = len(mobili)


outb  = mobili[0].split()[0]
inb   = mobili[0].split()[1]
hpMax = float(mobili[0].split()[6])
linemax = mobili[0]

i = 1
while i < nlin:

    outn = mobili[i].split()[0]
    inn = mobili[i].split()[1]

    if ((inn == inb) and (outn == outb)):
        hpn = float(mobili[i].split()[6])
        if hpn > hpMax:
            hpMax = hpn
            linemax = mobili[i]

    if ((inn != inb) or (outn !=outb)):
        mobiliok.write(linemax)
        linemax = mobili[i]
        hpMax = float(mobili[i].split()[6])
    i+=1
    outb = outn
    inb  = inn
mobiliok.write(linemax)
