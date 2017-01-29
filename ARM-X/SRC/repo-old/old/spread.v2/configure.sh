src=`pwd`

echo '-----------'
echo 'comandos: '
echo '------------------------------------------------------------------'

alias brk=$src/break.py
echo 'brk  - separa arquivo centros.dat'

alias chk=$src/checkoutv.py
echo 'chk  - checa com visual python descrição de cadeias isoladamente'

alias chka=$src/checkoutvALL.py
echo 'chk  - checa com visual python descrição de todas cadeias'


alias cpl="gfortran -o $src/stat.v8 $src/AngAneis_v8.0_pbc_beta.f90"
echo 'cpl  - compila o fonte AngAneis v8'

alias stat="$src/stat.v8"
echo 'stat - Versao 8 da analise estatistica'

echo '------------------------------------------------------------------'
echo 
echo 'versao 8/ Set.2010 / Analise estatistica.'
echo 

