#include <stdlib.h>
#include <stdio.h>
#include "boolean.h"
#include <assert.h>


#include <unistd.h>
#include "distrib.h"


#include "orbitools.h"


    #define JMAX 5
int verb = 0;

basis_t bder;
basis_t bres;


int   countfr = 0;
void *rootfr = NULL;
int invtfr( boole f )
{
    int tp[ ffsize ];
    int x;
    for( x = 0; x < ffsize; x++)
	    tp[x] = f[x] ? -1 : +1;

     int last = countfr;
     int val  = findtable(tp, ffsize, &rootfr, &countfr, 1);
     if (verb == 2   && val == last )
               printf("counttfr :%d\n", countfr);

    return val;
}

int derprepare(  int dim )
{
    initboole( dim );
    initagldim( dim );
    bder = monomialBasis(1, dim, ffdimen);
    loadBasename( "lift", &bder );
    return 1;
    int r = orbitBasic( mkaglGroup( dim ) , &bder);
    printf("\norbit(1,%d) : %d\n", r , dim);
    saveBasename( "lift", bder );
    return 2;
}

int resprepare(  int dim )
{
    initboole( dim );
    initagldim( dim );
    bres = monomialBasis(2, dim, ffdimen);
    loadBasename( "res", &bres );
    return 1;
    int r = orbitBasic( mkaglGroup( dim ) , &bres);
    printf("\norbit(1,%d) : %d\n", r , dim);
    saveBasename( "res", bres);
    return 2;
}
int countd = 0;
void *rootd = NULL;

int countr = 0;
void *rootr = NULL;

int countj= 0;
void *rootj = NULL;

int J( boole f ) 
{
	int R[ 2 ];
        R[ 0 ] = invSimpleDerivation( f , ffsize, &bder, &rootd, &countd);
        R[ 1] = invRestriction( f , ffsize, &bres, &rootr,  &countr );
	int val  = findspltable(R, 2, &rootj, &countj);
	return val;
}

int main(int argc, char *argv[])
{
    FILE *src;
    boole fct;
    int num = 0;
    char fn[64];
    int dim = 6;
    
    int opt, optdeg = 0, optfr=0, optall=0, optlift=0, optres=0;
    while ((opt = getopt(argc, argv, "ad:lvtr")) != -1) {
	switch (opt) {
	case 'a': optall = 1; break;
	case 'd': optdeg = atoi( optarg ); break;
	case 't': optfr  = 1;break;
	case 'r': optres  = 1;break;
	case 'l': optlift = 1;break;
	case 'v':
	    verb++;
	    break;
	default:		/* '?' */
	    fprintf(stderr, "Usage: %s [-a anf ] [-r log iter]\n",
		    argv[0]);
	    exit(EXIT_FAILURE);
	}
    }

   
    derprepare(  dim - 1 );
    resprepare(  dim - 1 );

    initboole( 6 );
    initagldim( 6 );
    sprintf(fn, "../data/class-2-%d.txt", ffdimen);
    src = fopen(fn, "r");
    if (!src) {
	perror("");
	exit(1);
    }
    uint64_t size;
    aglGroup g;
    int D[ 8 ] = {0};
    boole f = getboole();
    while ( (fct = loadaglboolesize( src, &g, &size ) )) 
    	if ( optdeg==0  || degree(fct) == optdeg ) {
	    	int tfr [ffsize];
		int x;
	    	for( x = 0; x < ffsize; x++ )
			tfr[x] = fct[x];	    
	   	 Fourier( tfr, ffsize );
	    	for( x = 0; x < ffsize; x++ )
                	if ( tfr[x] < 0 ) tfr[x] = (ffsize-1) * (-tfr[x]) ;
		printf("\ndeg %d :", degree(fct) );
		int r;
	    	for( r = 0; r <= ffdimen; r++) {
		    for( x = 0; x < ffsize; x++ )
			    f[x] = (( 1 << r ) & tfr[x] ) > 0;
		    printf("%3d", degree(f) );
		    D[r] = degree(f);
	    	    }
		r = 0;
		while ( D[r] <= D[r+1] ) r++;
		printf( "r=%d", r );
		//printf("\ninv %d :", J(fct) );
		printf("\ninv :" );
	    	for( r = 0; r <= ffdimen; r++) {
		    for( x = 0; x < ffsize; x++ )
			    f[x] = (( 1 << r ) & tfr[x] ) > 0;
		    printf(" %d", J(f) );
	    	    }
		}
	    num++;
	    free( fct );
	    aglfreeGroup( g );
    fclose(src);
    printf("\ncount  : %d \n", countj);
    return 0;
}
