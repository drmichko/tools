
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "boolean.h"


int main(int argc, char *argv[])
{
    int opt, count = 0, verbe = 0;
    while ((opt = getopt(argc, argv, "a:m:f:r:hiv:s:")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'v':
	    verbe++;
	    break;
	default:		/* '?' */
	    fprintf(stderr, "Usage: %s [-a anf ] [-r log iter]\n",
		    argv[0]);
	    exit(EXIT_FAILURE);
	}
    }

    boole f;
    for( int i = optind; i < argc; i++ )  {
    	FILE *src = fopen( argv[i], "r" );
        if ( ! src ) {
		perror( argv[i] );
		exit( 1 );
	}
	while ((f = loadBoole(src))) {
		pTT( stdout, f );
		free( f );
		count++;
	}
	fclose(src);
    }
    printf("\n#count=%d\n", count);
    return 0;
}
