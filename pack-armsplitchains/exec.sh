#!/bin/bash

pwd

#organiza
cd analisa
allpdb2xyz
#mkdir PDB
#mv *.pdb PDB
#mkdir MATs
#mv MAT-* MATs

#linka pra acertar indexacao
ct=0
for i in $( ls *.xyz); do
let ct=ct+1
ln $i stretch-$ct.xyz
done

#separa
cd ..; pwd;
cp ../info .
spread

# prepara
cp  ./analisa/MAT-*-000 P26V25STR/MAT.in
cp ../chains.in P26V25STR
cp ../info.stat P26V25STR
cp ../execall.in P26V25STR


#roda
pwd

cd P26V25STR
rm -f *.xyz

arms-rodageral.py 1



