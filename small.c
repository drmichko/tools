
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"
#include <omp.h>

#include "boolean.h"
#include "orbitools.h"

#include <stdatomic.h>   


int checks( boole f, agl L, int r )
{
    boole h = getboole( );
    for ( shortvec x = 0; x < ffsize; x++)
        h[x] = f[  aglImage( x,  L  ) ] ^ f[x];
    if( degree(h) >  r ) {
                printf("\ndegrees : f=%d  b=%d\n", degree(f),  degree(h) );
		panf( stdout, f );
		panf( stdout, h );
		puts("\n");
            return 0;
    }
    free( h );
return 1;
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

uint64_t nbfixed(const uint64_t *s, int dim)
{
    uint64_t x[64], y[64];
    for (int i = 0; i < dim; i++) {
        x[i] = s[i] ^ (1ULL << i);
        y[i] = (s[dim] >> i) & 1ULL;
    }
    int row = 0;
    for (int col = 0; col < dim && row < dim; col++) {
        int pivot = -1;
        for (int r = row; r < dim; r++)
            if ((x[r] >> col) & 1ULL) { pivot = r; break; }
        if (pivot < 0) continue;
        uint64_t t = x[row]; x[row] = x[pivot]; x[pivot] = t;
        t = y[row]; y[row] = y[pivot]; y[pivot] = t;
        for (int r = 0; r < dim; r++)
            if (r != row && ((x[r] >> col) & 1ULL)) {
                x[r] ^= x[row];
                y[r] ^= y[row];
            }
        row++;
    }
    int rank = 0;
    for (int r = 0; r < dim; r++) {
        if (x[r]) rank++;
        else if (y[r]) return 0;
    }
    return 1ULL << (dim - rank);
}
uint64_t mynbfixed( vector * s, int  dim )
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
        for( int l = i+1 ; l <  rank ; l++)
                if ( (bit & x[ l ])  ) {
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
return   ( ((uint64_t)1)  <<  (dim - rank) );
}

aglVector myboundaryAction(boole f, agl s, basis_t * b)
{
    // assuming  f o s = f + lower degree terms...
    int j;
    shortvec x;
    aglVector res = getaglVector( b->mat.nbl  + 1 );
    boole h = calloc( b->q, sizeof(uchar) );

    for (x = 0; x < b->q; x++)
        h[x] = f[  aglImage(x, s) ] ^ f[x];


    assert( degree(h) <=   b->t );
	projboole(b->s, b->t, h );
    res[ b->mat.nbl ] = booleVector( h , b );

    for (j = 0; j < b->mat.nbl; j++) {
        for (x = 0; x < b->q ; x++)
            h[x] = b->mat.fct[j][aglImage(x, s)];
	projboole( b->s, b->t, h );
        res[j] = booleVector(h, b );
        }

    free(h);
    return res;
}


int main(int argc, char *argv[])
{  

    char *fn = "../boole/data/B-4-5-7.dat";


    aglGroup g = NULL;
    int opt, optw = 0; 
    int dimen = 7;
    int target = 0;
    int nbthread = 1;
    while ((opt = getopt(argc, argv, "i:f:p:wb:d:m:r:t:")) != -1) {
	switch (opt) {
	case 'd':
	   break;
	case 'p':
	    nbthread = atoi( optarg );
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
    initBrowse( &base );
    while ((f = loadaglboolesize( src , &g, &size ))) {
		item++;
	        //panf( stderr, f );	
		fprintf(stderr,  "\nitem=%d degree=%d order=%ld\n", item, degree(f), size );
		//paglGroup( stderr, g );
                if ( item == target ) 
			break;
		aglfreeGroup( g );
		free(f);
		g=NULL;
    }
    if ( item != target ){
	    exit(0);
    }
    fclose(src);
    char name[64];
    sprintf( name, "rho/rho-%d.out", target );

    FILE* dst = stdout;
/*
    if ( ! optw || access( name, F_OK) == 0) {
	    puts("\nusing stdout");
    } else dst = fopen( name, "w");
  */
    clock_t debut = omp_get_wtime();
 
    panf( dst , f );
    fprintf( dst, "\nitem=%d degree=%d order=%ld\n", item, degree(f), size );
    paglGroup( stderr , g );
    assert( degree( f ) > 4  );

    assert(  checkstab( f, g,  3) == 1 );

    agl*  table = devGroup( g, &size );
    aglVectorGroup G = aglBoundaryGroupAction( f, g, &base );
    size_t  count = 0;

    fprintf( stderr, "\norder=%ld\n", size );
    aglVector *L = calloc( size, sizeof( aglVector ) );
    for( size_t i = 0;  i < size; i++ ) {
        assert(  checks( f, table[i],  3) == 1 );
    	L[i]   = myboundaryAction(  f, table[ i ], &base );
		//printf("X");pagl( stdout , table[i] );
	vector tmp = nbfixed( L[i] , 35  );
        if ( tmp ) 
	       fprintf( stderr, "#%ld tmp=%ld\n", i, tmp );
	count+= tmp;
    }
    fprintf( stderr, "\n#count=%ld\n", count);
    count /= size;
    fprintf( stderr, "\n#orbits=%ld\n", count);
    if ( count % size  )  
	    exit(0);
    return 1;
    omp_set_num_threads(nbthread);


uint64_t cpt = 0;
#pragma omp parallel for schedule(guided) reduction(+:cpt)
for( vector v = 0; v < base.size; v++ ){
    if ( base.table[v] ) continue;
	int is_rep = 1;
    	for( size_t i = 0; i < size; i++ ){
        	vector y = aglVectorImages( v, L[i], 35 );
		if ( y < v ) {                     // image strictement plus petite
            		is_rep = 0;
            		break;
        	}
        	if ( y > v ) {                     // y a un antécédent plus petit
            		base.table[y] = 1;
        	}
    	}
    	if ( is_rep ) {
	    cpt++; 
	    base.table[v] = 2;
    	}
   }

		
    fprintf( dst, "\n#orbits=%ld cpt=%ld", count, cpt  );
    clock_t fin = omp_get_wtime();
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
