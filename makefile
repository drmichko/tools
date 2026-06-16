
OPTION = -I../boole/src -L../boole/src


ifeq ($(DEBUG),1)
	CFLAGS = -Wall -g  $(OPTION) 
else

	CFLAGS = -O3  -march=native   -Wno-unused-result $(OPTION)

endif


all : schatz.exe action.exe stabredo.exe stab.exe ft.exe  nl.exe test.exe invariant.exe nnl.exe ab.exe anfload.exe print.exe anfsimple.exe rd.exe regroup.exe 
	

oldall : invariant.exe dyadic.exe nnl.exe ab.exe anfload.exe print.exe anfsimple.exe rd.exe
	
debug:
	$(MAKE)  -B DEBUG=1

schatz.exe : schatz.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

action.exe : action.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

stabredo.exe : stabredo.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

stab.exe : stab.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp
anfsimple.exe : anfsimple.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

anfload.exe : anfload.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp -lm

rd.exe : rd.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

nl.exe : nl.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp
ft.exe : ft.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

test.exe : test.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

regroup.exe : regroup.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp -lm

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

install :
	port.sh  -w tools -h rayol
clean :
	rm -f *.exe
