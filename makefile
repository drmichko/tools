CFLAGS = -Wall -g -I../boole/src -L../boole/src

all : anfsimple.exe
	
#alphasel.exe mapstat.exe select.exe last.exe check.exe testing.exe  

#header.exe anfsimple.exe building.exe statboole.exe orbital.exe pairing.exe 


anfsimple.exe : anfsimple.c
	gcc $(CFLAGS) $^  -o $@  -lboole  -lgmp

alphasel.exe : alphasel.c
	gcc $(CFLAGS) $^  -o $@  -lboole
building.exe : building.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

initmap.exe : initmap.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

clean.exe : clean.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

testing.exe : testing.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

merge.exe : merge.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

last.exe : last.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

reduce.exe : reduce.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp


filter.exe : filter.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp



check.exe : check.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

mapstat.exe : mapstat.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

select.exe : select.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp
final.exe : final.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp -lm

pairing.exe : pairing.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

header.exe : header.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp

orbital.exe : orbital.c
	gcc $(CFLAGS) $^  -o $@  -lboole -lgmp


statboole.exe : statboole.c
	gcc $(CFLAGS) $^  -o $@  -lboole

clean :
	rm -f *.exe
