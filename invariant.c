#include <stdlib.h>
#include <stdio.h>
#include "boolean.h"
#include <assert.h>


#include <unistd.h>
#include "distrib.h"


#include "orbitools.h"


    #define JMAX 5
int verb = 0;
void *rootj = NULL;
int  countj=0;

basis_t bder;
basis_t bres;

void showdtable(int *t, int q)
{
    int i;
    for (i = 0; i < q; i++)
        if (t[i])
            printf("%3d [%d]", t[i], i);
}


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
    //loadBasename( "lift", &bder );
    //return 1;
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

int main(int argc, char *argv[])
{
    FILE *src;
    boole f;
    int num = 0;
    char fn[64];
    int dim = 6;

    int opt, optdeg = 0, optfr=0, optall=0, optlift=0, optres=0;
    while ((opt = getopt(argc, argv, "adlvtr")) != -1) {
	switch (opt) {
	case 'a': optall = 1; break;
	case 'd': optdeg = 1; break;
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

   
    if ( optlift) derprepare(  dim -1 );
    if ( optres ) resprepare(  dim -1 );

    initboole( 6 );
    initagldim( 6 );
    sprintf(fn, "../boole/data/class-2-%d.txt", ffdimen);
    src = fopen(fn, "r");
    if (!src) {
	perror( fn );
	exit(1);
    }
    uint64_t size;
    aglGroup g;
    int R[ JMAX ];
    int nbj = 0;
    while ( (f = loadaglboolesize( src, &g, &size ) )) {
	    nbj = 0;
	    if ( optall || optdeg ) {
		    R[ nbj++] = degree( f );
	    }
	    if ( optall || optfr ) {
		    R[ nbj++] = invtfr( f );
	    }
	    if ( optall || optlift ) {
		    R[ nbj++] = invSimpleDerivation( f , ffsize, &bder, &rootd, &countd);
	    }
	    if ( optall || optres ) {
		    R[ nbj++] = invRestriction( f , ffsize, &bres, &rootr,  &countr );
	    }
	    int last = countj;
	    int val  = findspltable(R, nbj, &rootj, &countj);
	    num++;
	    if (verb  && val == last )
	       printf("countj: %d (%d)\n", countj, num);
	
	    free( f );
	    aglfreeGroup( g );
    }
    fclose(src);
    printf("countfr: %d\n", countfr);
    printf("count  : %d\n", countj);
    printf("maps   : %d\n", num);
    return 0;
}
