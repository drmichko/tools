for x in $*
do
	./rho.exe  -f ../boole/data/B-4-5-7.dat -t $x -w   &
done
