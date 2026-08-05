#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"
#include "code.h"
#include "distrib.h"


int main(int argc, char *argv[])
{
    FILE *src = NULL;
    int opt;

    while ((opt = getopt(argc, argv, "a:m:f:r:hiv:s:")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'f':
	    src = fopen(optarg, "r");
	    if (!src) {
		perror(optarg);
		return 1;
	    }
	    break;

	}
    }

   int all = 0;
   boole  f;
    while ((f = loadBoole(src))) {
	    initboole( ffdimen-1);
	    panf( stdout, f );
	    initboole( ffdimen+1);
	    all++;
	    free(f);
    }
    fclose(src);

    printf("\nall=%d\n", all );



    return 0;
}
