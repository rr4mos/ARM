#!/bin/bash

srcsp=`pwd`

BINsp="$srcsp/spread.v2/"

echo $BINsp

echo 'scripts em: '$BINsp

echo '#ajustando inputs'
$BINsp/adj.py

echo '# separacao com Hs e sem edicao.'
$BINsp/spread.py info chains

echo '# edicao dos outputs < reorganizado.pdb'
echo '# tirando Cs e corrigindo numeracao'
$BINsp/renum.py

echo '# separacao com Hs e sem edicao.'
$BINsp/spread.py info2 chainsout

echo '#limpando...'
rm -f dados
rm -f .wcl
