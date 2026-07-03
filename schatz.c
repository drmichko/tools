
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"
#include "distrib.h"

aglGroup grp ;
int64_t size ;

int *table( boole f )
{ code cc = rmcode( 2, 2, 7 );
  int * A = calloc( 1 << cc.nbl , sizeof(int ) );
	
	int  cpt=1, limite = 1 << cc.nbl;
        boole  t = getboolecpy( f );
        A[0] = ( ffsize - linearity( f ) ) / 2;
        int q = 0;
        while ( cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= cc.fct[i][x];
                q ^=  ( 1 << i );
                A[ q  ] = ( ffsize - linearity( t ) ) /2;
                cpt++;
        }

        free( t );
        freecode( cc );
      	pdistrib( "\n W=", A, limite  ); 
  return A;
}


int test( int* F, int * G)
{
	int limite = rmdimen( 2, 2, 7 );
	limite = 1 << limite;
	for( int g = 0; g < limite; g++ ) {
		int q = 0;
		while ( (F[ q ]  + G[ q ^ g]  > 88 ) &&  ( q < limite ) ) q++; 
		if ( q  == limite ) {
		       	puts(" YES !!! ");
			exit(0);
		}
	}
}

int main(int argc, char *argv[])
{  

    char *fn = NULL;


    aglGroup g = NULL;
    int num = 0, cls = -1;
    int opt, optw = 0, optr=-1, optinit = 0;
    int deg = 0, dimen = 7;
    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:")) != -1) {
	switch (opt) {
	case 'w':
	    optw++;
	    break;
	case 'd':
	    deg = atoi(optarg);
	   break;
	case 'c':
	    cls = atoi(optarg);
	    break;
	case 'i': 
		optinit = atoi( optarg);
	break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    optr=atoi( optarg );
	    break;
	default:
	    exit(0);
	}

    }


    initboole(  dimen );
    initagldim( dimen );

    boole f;
    num = 0;
    uint64_t grpSize, orbSize;
    int val = 0;
    
    FILE * src = fopen(  fn , "r" );
    if ( ! src ) {
		perror( fn  );
		exit(1);
	}
    
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
		 panf( stdout, f );
		 int * F = table( f );
	    	 boole g;
	    	 int nb = 0;
    	    	 while ((g = loadBoole(src ) ) ) {
			for( int x = 0; x < ffsize; x++ )
				g[x] ^= f[x];
		 	int * G = table( g );
			test( F, G );
			free( G );
	    		free( g );
			nb++;
	    	}
	        printf("\n#nb: %d\n", nb);
	    free(f);
            aglfreeGroup( grp );
	    num++;
        }

     fclose(src);

     printf("\n#maps : %d\n", num );
     return 0;
}
