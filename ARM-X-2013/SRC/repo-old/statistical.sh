#!/bin/bash

#curr=`pwd`
#cd ../../
#source configure.sh
#cd $curr

for DIR in $( ls );
do

cp ../MAT.in     $DIR
cp ../info.stat  $DIR
cp ../chains.in  $DIR

cd $DIR
../../../SRC/stat.v10

echo $DIR'  --> ok'
cd ../
done
