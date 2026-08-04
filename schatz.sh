#!/usr/bin/bash

ok=0
for file in results/*
do
	if ! grep -q  $file goodies.txt
	then
		ok=1
		break
	fi
done
if [ $ok = 1 ]; then
	echo $file
fi

nb=$( nproc )

if [ ! -f $file ]; then
	exit
fi

if grep $file goodies.txt ; then
	echo goodies ?
	exit
fi

echo $file  
echo $file  >> goodies.txt


for(( j = 0; j < nb ; j++ )) ; do 
	./schatz.exe -f $file -j$j:$nb & 
done

wait

echo $file  : done >> goodies.txt

