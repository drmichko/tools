CFLAGS = -Wall -g -I../boole/src -L../boole/src

all : anfload.exe print.exe anfsimple.exe rd.exe
	
anfsimple.exe : anfsimple.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

anfload.exe : anfload.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp -lm

rd.exe : rd.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

print.exe : print.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp



clean :
	rm -f *.exe
