
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"

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
    initBrowse( &base );
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
    char name[64];
    sprintf( name, "rho/rho-%d.out", target );


    FILE* dst = stdout;
    if ( access( fn , F_OK) == 0) {
	    puts("done ?");
    } else dst = fopen( name , "w");
    clock_t debut = clock();
 
    fprintf( dst, "item=%d degree=%d order=%ld\n", item, degree(f), size );
    panf( dst , f );
    assert( degree( f ) > 4  );
    aglVectorGroup lda    = aglBoundaryGroupAction( f, g, &base );
    size_t  count = 0;
    for( vector v = 0; v < base.size; v++ ) {
		if ( ! base.table[ v ] ) {
			base.table[ v ] = 2;
			count++;
    			if ( v % (1<<30) == 0 )
				fprintf( dst, "\n#v=%ld count=%ld", v, count );
		}
         	aglVectorGroup aux = lda;
        	while (aux != NULL) {
            		uint64_t  y  = aglVectorImages( v  , aux->act, 35);
			if ( y > v ) base.table[y] = 1;
            		aux = aux->next;
                }
            }
    fprintf( dst, "\n#count=%ld", count );
    clock_t fin = clock();
    double temps = (double)(fin - debut) / CLOCKS_PER_SEC;
    fprintf(dst, "\nruntime : %.6f secondes\n", temps);
    char bname[64];
    sprintf( bname, "%d", item );
    saveBasename( bname , base );
    aglfreeGroup( g );
    aglVectorGroupFree(  lda );
    free( base.table );
    
    free( f );
    fclose(dst);
    return 0;
}
