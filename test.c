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

void doit( boole f, boole g )
{
code Q = rmcode( 2, 2, 5 );
int cpt = 1, limite = 1 << Q.nbl;
int * a = calloc( limite, sizeof(int ) );
boole  t = getboolecpy( f );
int wt = ( ffsize - linearity( f ) ) / 2;
a[ 0 ] = wt;
while (  cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= Q.fct[ i ][x];

                wt = ( ffsize - linearity( t ) ) / 2;;
                a[ cpt ] = wt;
                cpt++;
        }

pdistrib("f", a, limite);
printf("\n");
}



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


    boole f;
    boole t[2];
    int i = 0;
	while ((f = loadBoole(src))) {
	    panf( stdout, f  );
	    putchar('\n');
	    t[ i++ ] = f;; 
	}
    fclose(src);

    doit( t[0], t[0] );
    doit( t[0], t[1] );
    doit( t[1], t[1] );

    return 0;
}
