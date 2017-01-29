#!/bin/bash

rm -rf chain-$1
rm -f chain-$1.tar
rm -f chain-$1-all.xyz

mkdir chain-$1
ct=1
for i in $(ls */P26V25STR/*/chainsout/$1.xyz); do
	mct=$(printf %03d "$ct")
	cp $i ./chain-$1/$mct.xyz
	cat $i >> chain-$1-all.xyz
	let ct=$ct+1
done

tar cvf chain-$1.tar chain-$1


