for x in $*
do
	./small.exe  -f ../boole/data/B-4-5-7.dat -t $x -w  -p64  &
done
