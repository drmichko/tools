#!/usr/bin/bash
#
iter=20
while [ $iter -gt 0 ];
do
for (( num=178; num > 0; num-- ))
do
	file=results/res-$num.txt
	if [ ! -f  $file ] ; then
		break
	fi
done

	if [  $num =  0 ] ; then
		exit
	fi

echo iter=$iter  $file

file=data/$num-3-3-7.zip
if [ !  -f  $file  ]; then
	echo making $file takes a while...
	./ziptable.exe  $num   >  results/res-$num.txt 
fi

./nlqfast.exe  -t  $num   -g38        
sleep 10
let iter--
done
