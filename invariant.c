#include <stdlib.h>
#include <stdio.h>
#include "boolean.h"
#include <assert.h>


#include <unistd.h>
#include "distrib.h"


#include "orbitools.h"


typedef struct {
	int     num;
	int     cpt;
	boole   fct;
	int64_t size;
} repres;

repres * table;

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

int countd = 0;
void *rootd = NULL;

int countr = 0;
void *rootr = NULL;

int main(int argc, char *argv[])
{
    FILE *src = NULL;
    boole f;
    int num = 0;
    char fn[64];
    int dim = 6;

    int opt, optdeg = 0, optfr=0, optall=0, optlift=0, optres=0;
    while ((opt = getopt(argc, argv, "adlvtrf:")) != -1) {
	switch (opt) {
	case 'a': optall = 1; break;
	case 'd': optdeg = 1; break;
	case 't': optfr  = 1;break;
	case 'r': optres  = 1;break;
	case 'l': optlift = 1;break;
        case 'f': src = fopen( optarg, "r" ); 
    		if ( ! src ) {
		perror( optarg );
		exit(1);
    		}
		 break;
	case 'v':
	    verb++;
	    break;
	default:		/* '?' */
	    fprintf(stderr, "Usage: %s [-a anf ] [-r log iter]\n",
		    argv[0]);
	    exit(EXIT_FAILURE);
	}
    }

   
    if ( optlift) derprepare(  &bder, dim -1 );
    if ( optres ) resprepare(  &bres, dim -1 );

    initboole( 6 );
    initagldim( 6 );
    sprintf(fn, "../boole/data/class-2-%d.txt", ffdimen);
    if ( ! src )  src = fopen(fn, "r");
    if ( ! src ) {
	perror( fn );
	exit(1);
    }
    table = calloc( 160000, sizeof(repres) );
    uint64_t size = 0;
    int R[ JMAX ];
    int nbj = 0;
    //while ( (f = loadBooleStab( src , &size ) )) {
    while ( (f = loadBoole( src) )) {
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
	    int val  = findspltable(R, nbj, &rootj, &countj);
	    if ( table[ val ].cpt == 0  ) 
		   panf( stdout, f);
	    else   {
		    panf( stdout, f);
		    printf(" [doublon]");
	    }
	    table[ val ].num = num;
	    table[ val ].size = size;
	    table[ val ].fct  = f;
	    table[ val ].cpt++;
	    num++;
    }
    fclose(src);
    printf("\ncountfr: %d\n", countfr);
    printf("count  : %d\n", countj);
    printf("maps   : %d\n", num);

   for (int i = optind; i < argc; i++) { 
	   src = fopen( argv[i], "r" );
    		if ( ! src ) {
		perror( fn );
		exit(1);
    		}
    	   while ( (f = loadBoole( src ) )) {
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
            int val  = findspltable(R, nbj, &rootj, &countj);
	    table[ val ].cpt++;
            free( f );
	   }
	   fclose( src );
   }

   if ( optind != argc )
   for( int i = 0; i < 160000 ; i++ )
	   if ( table[i].cpt > 1 ) {
	    	   boole f = table[i].fct;
		   printf("\n\nnum=%d mult=%d", table[i].num, table[i].cpt -1 );
		   printf(" alpha=%.2f deg:%d size=%ld", alpha(f), degree(f), table[i].size ) ;
		   panf( stdout, f  );
	   }

    printf("\ninv counter=%d", countj );
    return 0;
}
