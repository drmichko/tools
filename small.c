
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"



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

int aglrankvector (vector *s, int dim )
{   
int i, l;
vector  bit, x[ dim ]; 
for( i = 0;  i < dim; i++ )
	x[i] = s[i];
i=0;
while ( i < dim ) {
  if ( x[i] ) {
        bit = x[i];
        bit = bit ^ ( bit & ( bit - 1) );
        for( l = 0; l < dim; l++)
          if ( l != i )
                if ( bit & x[l] ) x[l] ^= x[i];
        i++;
   } else {
        dim = dim - 1;
        x[i] = x[ dim ];
   }
 }
return dim;
}

uint64_t nbfixed( vector * s, int  dim )
{   
uint64_t   x[ dim ], y[ dim ]; 
int i,  rank = dim;
for( i = 0;  i < dim; i++ ){
	x[i] = s[i]   ^   ( ( 1ULL <<i ) )  ;
        y[i] = ( s[ dim ] >> i ) & 1 ;
}
i = 0;
while ( i < rank  ) {
  if ( x[ i ] ) {
        int j = __builtin_ctzll( x[i] );
	uint64_t bit = ( ( 1ULL << j ) );
        for( int l = i + 1 ; l <  rank ; l++)
                if ( bit & x[ l ] ) {
			x[ l ] ^= x[ i ];
			y[ l ] ^= y[ i ];
			}
        i++;
   } else {
	if (  y[i] ) return 0;
        rank = rank - 1;
        x[ i ] = x[ rank];
	y[ i ] = y[ rank];
   }
 }
return   ( 1ULL  <<  (dim - rank) );
}



int main(int argc, char *argv[])
{  

    char *fn = "../boole/data/B-4-5-7.dat";


    aglGroup g = NULL;
    int opt, optw = 0; 
    int dimen = 7;
    int target = 0;

    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:")) != -1) {
	switch (opt) {
	case 'd':
	   break;
	case 'c':
	    break;
	case 'w':
	   	optw++; 
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
    char name[64];
    sprintf( name, "rho/rho-%d.out", target );

    FILE* dst = stdout;

    if ( ! optw || access( name, F_OK) == 0) {
	    puts("using stdout");
    } else dst = fopen( name, "w");
    clock_t debut = clock();
 
    fprintf( dst, "item=%d degree=%d order=%ld\n", item, degree(f), size );
    panf( dst , f );
    assert( degree( f ) > 4  );
    initBrowse( &base  );

    assert(  checkstab( f, g,  3) == 1 );

    agl*  table = devGroup( g, &size );
    aglVectorGroup G = aglBoundaryGroupAction( f, g, &base );
    size_t  count = 0;

    fprintf( stderr, "\norder=%ld\n", size );
    aglVector *L = calloc( size, sizeof( aglVector ) );
    for( size_t i = 0;  i < size; i++ ) {
    	L[i]   = boundaryAction(  f, table[ i ], &base );
	vector tmp = nbfixed( L[i] , 35  );
	count+= tmp;
    }
    assert( count %size == 0 );
    count /= size;
    fprintf( stderr, "\n#orbits=%ld\n", count);

    uint64_t   cpt = 0;
#pragma omp parallel for schedule(dynamic) reduction(+:cpt)
    for( vector v = 0; v < base.size; v++ ){
        size_t i = 0;
        for( i = 0 ; (i < size) &&  v < aglVectorImages( v, L[i], 35 );  i++ );;
        int accept = ( i == size );
        base.table[v] = accept;
        cpt += accept;
    }
		
    fprintf( dst, "\n#orbits=%ld cpt=%ld", count, cpt  );
    clock_t fin = clock();
    double temps = (double)(fin - debut) / CLOCKS_PER_SEC;
    fprintf(dst, "\nruntime : %.6f secondes\n", temps);

    return 0;
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
