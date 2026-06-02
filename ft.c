#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"
#include "math.h"
#include "code.h"


int limite;
code cc[8];

int rho[4][9];



int *F[ 9 ];

void mktable( void )
{
for( int v = 0 ; v <= ffdimen; v++ ) 
	F[ v ] = calloc( ffsize  , sizeof(int) );
}

int count ;

int Z[ 8] = { 0 };

#define FF(v, s, w) ( F[(v)][ ((s) << (v)) + (w) ] )

void inittable( boole f )
{
count = 0;
for( int s = 0; s < ffsize ; s++ ) 
		FF( 0, s, 0 ) = ( f[ s ]  ) ? -1 : 1;
return;
for( int s = 0; s < ffsize/2 ; s++ ) 
   for( int t = 0; t < 2; t++ ){
		FF( 1, s, t ) = ( f[ (s << 1) + t ] ^ t  ) ? -1 : 1;
}
}

int score = 0;

void zzz( boole f ) {
boole q = getboole( );
for( int v = 1; v <= ffdimen; v++ )
	for( int j = 0; j < v-1; j++ )
		if ( Z[ v ] & ( 1 << j ) )
			q[ ( 1 << (v-1) ) + ( 1 <<  j ) ]= 1;
TTtoANF( q );
for( int x = 0; x < ffsize; x++ )
	q[x] ^= f[x];
int tmp = linearity( q );
if ( tmp >  score ) {
	score = tmp;
	count = 0;
	printf("\nscore=%d", score );
}
if ( tmp == score) count++;
free(q);
}

int mygamma( boole f , int v, int s ) {
boole q = getboole( );
int size =  1 << v;
for( int i = 1;  i <= v; i++ )
	for( int j = 0; j < v-1; j++ )
		if ( Z[ i ] & ( 1 << j ) )
			q[ ( 1 << (i-1) ) + ( 1 <<  j ) ]= 1;
xform( q , size );
int tfr[256];
for( int x = 0; x < size; x++ )
	q[x] ^= f[ (s << v )   + x ];
for( int x = 0; x < size; x++ )
	tfr[x] = 1 - 2*q[x];

Fourier( tfr, size  );
int res = 0;
for( int x = 0; x < size; x++ )
   if ( abs( tfr[x] )  > res ) res = abs( tfr[x] ); 

free(q);
return res;
}

void doit( int v, boole f, int R   )
// au niveau v, on  liste  (v-1)-bit Q candidat...
// F est connu pour les niveaux < v
// propager au niveau v 
{
if ( v == ffdimen ) {
	zzz( f );
	return;
}
int * M = calloc( ffsize, sizeof(int ) );

//show( v - 1 );
int ret = 1 << (ffdimen - v );
for( int s = 0; s <  ( 1 << ( ffdimen - v - 1) ); s++ ) {
	for( int Q = 0;  Q< ( 1 << ( v - 1) ) ; Q++ ) {
		int max = 0;
		for( int u = 0;  u < ( 1 << ( v-1 ) ) ; u++ ) {
			int tmp = abs( FF( v-1, s,  u )  )  + abs( FF( v-1, s + ret , u^Q ) );    // t ? w ?
			if ( tmp > max ) max = tmp;
		}
		M[  ( s << ( v - 1 ) )  +  Q  ] = max;
	}
}
for ( int Q = 0; Q < ( 1 << (v-1) ) ; Q++ ) {
	int Gamma = 0;
 	 for( int s = 0; s < 1 << ( ffdimen - v ) ; s++ )
		Gamma += M[    ( s << ( v-1)  )  +  Q ];
	Z[ v ] = Q;
        int G = 0;
	for( int s = 0; s < 1 << ( ffdimen - v ) ; s++ )
		G += mygamma( f, v, s );
        printf("v=%d G=%d E=%d\n", v, G, Gamma ); 
        if ( G >= R && Gamma < R ) { printf("v=%d G=%d Gest=%d\n", v, G, Gamma ); }
	if (  G  >= R ) {
		int bit = 1 << ( v - 1 );
		for( int s = 0; s < 1 << ( ffdimen - v ) ; s++ ) {
		   for( int u = 0; u < 1 << ( v-1   ) ; u++ )
			FF(v, s, u )  =  FF( v-1, s, u ) +  FF( v-1, ret + s , Q ^ u );
		   for( int u = 0; u < 1 << ( v-1  ) ; u++ )
		   	FF( v, s, u+bit )  =  FF( v-1,        s, u  ) - FF( v-1, ret + s, Q ^ u  );
		}
		Z[ v ] = Q;
		doit( v+1, f, R  );
	}
}
free( M  );
}



int main(int argc, char *argv[])
{
    FILE *src = NULL;
    int opt,  R = 0;

    while ((opt = getopt(argc, argv, "a:k:m:f:R:hi:vsl:")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'R':
	    R = atoi(optarg);
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

    mktable( );
    boole f;
    printf("\nepsilon=%.2f",  (ffsize -R ) / 2.0 / 2 / ffsize ); 
    if (src) {
	while ((f = loadBoole(src))){ 
            if ( valuation(f) >= 0 ) 	{
		    inittable( f  );
		    score= 0;
		    doit( 1, f, R  );
	  	    printf("\nscore=%d %d", score , ffsize / 2 - score / 2);
		    printf( "\ncount=%d\n", count );
		   
	    }
	 free( f );
	//printf("\rtour=%d (%d)", tour, soluce  );
	fflush(stdout );
	}
	fclose(src);
    }

    printf("\n");
    printf("\n");
    return 0;
}
