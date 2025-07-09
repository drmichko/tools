CFLAGS = -Wall -g -I../boole/src -L../boole/src

all : invariant.exe dyadic.exe nnl.exe ab.exe anfload.exe print.exe anfsimple.exe rd.exe
	
anfsimple.exe : anfsimple.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

anfload.exe : anfload.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp -lm

rd.exe : rd.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

nnl.exe : nnl.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp -lm

ab.exe : ab.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp -lm

print.exe : print.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

invariant.exe :  invariant.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

dyadic.exe :  dyadic.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

clean :
	rm -f *.exe
