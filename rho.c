
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"

aglGroup grp ;
int64_t size ;
code rmh;


typedef struct info {
      boole f;
      int*  a;
      int deg;
      uint64_t stab;
      aglGroup grp;
} info;

info table[ 64 ];

void NL2( int item, uint64_t size , boole f )
{
        boole  t = getboolecpy( f );
        int *  a = calloc( ffsize+1, sizeof( int ) );
        int wt = ( ffsize - linearity( f ) ) / 2;
        a[ wt ]++;
        int nl2 = wt;
        int  cpt=1, limite = 1 << rmh.nbl;
        while ( cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= rmh.fct[i][x];

                wt = ( ffsize - linearity( t ) ) / 2;;
                a[ wt ]++;
		if ( wt < nl2) nl2 = wt;
                cpt++;
        }

        printf("item=%2d deg=%d  stab=%ld NL2=%2d :", item, degree( f  ), size, nl2 );
        for( int i = 0; i <=ffsize; i++ )
                        if ( a[i] ) printf(" %d [ %d ]", a[ i ], i );
                putchar('\n');
        table[item].a   = a;
	table[item].f   = f;
	table[item].stab   = size;
	table[item].deg = degree( f );;
        free( t );
}

void doit( int item, uint64_t size ,  boole f )
{
	
}


int checkstab( boole f, aglGroup g, int r )
{
while ( g ) {
    boole h = getboole( );
    for ( shortvec x = 0; x < ffsize; x++)
        h[x] = f[  aglImage(x, g->per ) ] ^ f[x];
    if( degree(h) >  r ) {
                printf("\ndegrees : f=%d  b=%d\n", degree(f),  degree(h) );
		panf( stdout, f );
		panf( stdout, h );
		puts("\n");
            return 0;
    }
    free( h );
    g = g -> next;
}

return 1;
}


int main(int argc, char *argv[])
{  

    char *fn = "../boole/data/B-4-5-7.dat";


    aglGroup g = NULL;
    int opt; 
    int dimen = 7;
    int target = 0;

    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:")) != -1) {
	switch (opt) {
	case 'w':
	    break;
	case 'd':
	   break;
	case 'c':
	    break;
	case 'i': 
	break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    break;
	case 't':
	    target=atoi( optarg );
	    break;
	default:
	    exit(0);
	}

    }


    FILE *src = fopen( fn, "r" );
    initboole(  dimen );
    initagldim( dimen );

    boole f;
    int item = 0;
    uint64_t size = 0;
    basis_t base = monomialBasis( 3, 3, 7 );

    while ((f = loadaglboolesize( src , &g, &size ))) {
		item++; 
    	        if ( ! target )
			printf( "item=%d degree=%d order=%ld\n", item, degree(f), size );
                if ( item == target ) break;
		aglfreeGroup( g );
		free(f);
    }
    if ( item != target ){
	    exit(0);
    }
    fclose(src);
    char name[128];
    sprintf( name, "rho/rho-%d.out", target );

   
    if ( target && ( access( name  , F_OK) == 0 ) ) {
	    puts("done ?");
	    exit(0);
    }

    FILE* dst;
   
   if ( ! target ) dst = stdout;
   else dst = fopen( name , "w");

    clock_t debut = clock();
 
    fprintf( dst, "item=%d degree=%d order=%ld\n", item, degree(f), size );
    panf( dst , f );
   // assert( degree( f ) > 4  );
    initBrowse( &base  );

    aglVectorGroup G = aglBoundaryGroupAction( f, g, &base );
    size_t  count = 0;
    for( vector v = 0; v < base.size; v++ )
	if ( ! base.table[v] ) {
		browse( v, G );
                base.table[ v ] = 1;
		count++;
		}
		
    fprintf( dst, "\n#orbits=%ld", count );
    for( vector v = 0; v < base.size; v++ )
	    base.table[v] =  ( base.table[v] == 1 );
    clock_t fin = clock();
    double temps = (double)(fin - debut) / CLOCKS_PER_SEC;
    fprintf(dst, "\nruntime : %.1f secondes\n", temps);
    char bname[64];
    sprintf( bname, "%d", item );
    saveBasename( bname , base );
    aglfreeGroup( g );
    aglVectorGroupFree(  G );
    free( base.table );
    
    free( f );
    fclose(dst);
    return 0;
}
