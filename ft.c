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
typedef  struct _list_ {
	 boole  fct;
	 int wt;
	 struct _list_ * next;
} enrlist, *liste;
int limite;
code cc[8];

int rho[4][9];


void append( int wt, int m, boole f, liste *l )
{  liste aux = malloc( sizeof( *aux )  );
   aux->fct = calloc( 1 << m, 1 );
   for( int x = 0; x < (1<<m) ; x++ )
	   aux->fct[x] = f[ x ];
   aux->wt = wt;
   aux->next = *l;
   *l = aux;
}

void glue( int wt, int m, boole L, boole R,  liste *l )
{  liste aux = malloc( sizeof( *aux )  );
   int q = 1 << m;
   aux->fct = calloc( 2 * q, 1 );
   for( int x = 0; x < q ; x++ ){
	   aux->fct[x] = L[ x ];
	   aux->fct[x+q] = R[ x ];
   }
   aux->wt = wt;
   aux->next = *l;
   *l = aux;
}

int length( liste l )
{ int r = 0;
	while ( l ) {
		r++;
		l=l->next;
	}
	return r;
}

int freeliste( liste l )
{ int r = 0;
	while ( l ) {
		liste tmp = l;
		r++;
		l=l->next;
		free( tmp->fct );
		free( tmp );
	}
	return r;
}

int check( liste l , int R )
{ 
	while ( l ) {
		if ( l->wt < R ) return 0;
		l = l->next;
	}
	return 1;
}
int showliste( char *msg, liste l  )
{ int r = 0;
	printf("\n%s:", msg);
	while (  l ) {
		r++;
		printf(" %d", l->wt );
		l = l->next;
	}
	printf("\ncard:%d\n", r);
	return r;
}
int wtboole ( boole f, int m )
{
	int wt = 0;
	for( int x = 0; x < ( 1 << m ); x++ )
                wt += f[x];
	return wt;
}

int  addwtboole ( boole f, boole g, int m )
{
	int wt = 0;
	for( int x = 0; x < ( 1 << m ); x++ ) {
                f[x] ^=  g[x];
		wt += f[x]; 
	}
	return wt;
}


boole* rm[ 9 ];

boole *  mkcode ( int k, int v  )
{       
        code    c = rmcode( k , k, v );
        boole   t = getboole( );
        int  cpt = 1, limite = 1 << c.nbl;
	boole * res = calloc(  limite, sizeof( boole ) ) ;
        res[ 0 ] = getboolecpy( t );
        while (  cpt < limite ) {
                int i = __builtin_ctz( cpt  );
        	addwtboole( t, c.fct [i], v );
                res[ cpt ] = getboolecpy( t );
                cpt++;
        }
	free( t );
	printf("rm(%d,%d).", k, v );
	fflush(stdout );
	return res;
}


void initrmcode( int k, int m )
{
	for( int v = k - 1; v < m; v++ ) {
		rm[ v ]  = mkcode( k - 1, v );

	}
}

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
}

void show( int v )
{
printf("\nv=%d : ", v );
for( int x = 0; x < ffsize; x++ )
	printf("%3d", F[v][x] );
printf("\n");
}

void mshow( int v, int *M )
{
printf("\nm=%d : ", v );
for( int x = 0; x < ffsize; x++ )
	printf("%3d", M[x] );
printf("\n");
}

int  FFF( int v, int s, int w )
{
    int t =  ((s) << (v)) + (w);
    if ( t >= ffsize ) {
	    printf("\nv=%d s=%d w=%d", v, s, w );
	    exit(0);
    }
    return  FF(v, s, w);
}

int score = 0;

void zzz( boole f ) {
boole q = getboole( );
for( int v = 1; v <=ffdimen; v++ )
	for( int j = 0; j < v-1; j++ )
		if ( Z[ v ] & ( 1 << j ) )
			q[ (1<<(v-1)) + ( 1 <<  j ) ]= 1;
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
}
void doit( int v, boole f, int R, int k  )
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
	for( int w = 0; w < ( 1 << ( v - 1) ) ; w++ ) {
		int max = 0;
		for( int t = 0; t < ( 1 << ( v-1 ) ) ; t++ ) {
			int tmp = abs( FF( v-1, s, w )  )  + abs( FF( v-1, s + ret , w^t ) );
			if ( tmp > max ) max = tmp;
		}
		M[  ( s << ( v - 1 ) )  + w ] = max;
	}
}
//mshow( v, M  );
for ( int Q = 0; Q < ( 1 << (v-1) ) ; Q++ ) {
	int Gamma = 0;
	for( int s = 0; s < 1 << ( ffdimen - v ) ; s++ )
		Gamma += M[    ( s << ( v-1)  )  + Q ];
	if ( Gamma >= R ) {
		int bit = 1 << ( v - 1 );
		for( int s = 0; s < 1 << ( ffdimen - v ) ; s++ ) {
		   for( int u = 0; u < 1 << ( v-1   ) ; u++ )
			FF(v, s, u )  =  FFF( v-1, s, u ) +  FFF( v-1, ret + s , Q ^ u );
		   for( int u = 0; u < 1 << ( v-1  ) ; u++ )
		   	FF( v, s, u+bit )  =  FFF( v-1,        s, u  ) - FFF( v-1, ret + s, Q ^ u  );
		}
		Z[ v ] = Q;
		doit( v+1, f, R, k );
	}
}
free( M  );
}



int main(int argc, char *argv[])
{
    FILE *src = NULL;
    char *anf = NULL;
    int opt,  R = 0;
    int k = 2;
    limite = 5;

    while ((opt = getopt(argc, argv, "a:k:m:f:R:hi:vsl:")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'a':
	    anf = strdup(optarg);
	    break;
	case 'l':
	    limite = atoi(optarg);
	    break;
	case 'R':
	    R = atoi(optarg);
	    break;
	    case 'k' : k = atoi( optarg) ; break;
	case 'f':
	    src = fopen(optarg, "r");
	    if (!src) {
		perror(optarg);
		return 1;
	    }
	    break;

	}
    }

    initrmcode( k, ffdimen );
    mktable( );
    boole f;
    if (src) {
	while ((f = loadBoole(src))){ 
            if ( valuation(f) >= 0 ) 	{
		    inittable( f  );
		    boole g = getboole(  );
		    score= 0;
		    doit( 1, f, R, k );
	  	    printf("\nscore=%d %d", score , ffsize / 2 - score / 2);
		    printf( "\ncount=%d\n", count );
		    free( g );
		   
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
