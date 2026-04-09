#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"
#include "code.h"
#include "distrib.h"

basis_t base;

void  NLL( boole f , int goal)
{
        code   c = rmcode( 1 , 2, ffdimen  );
        boole  t = getboole( );
        int *  a = calloc( ffsize+1, sizeof( int ) );
        int wt ; // = ( ffsize - linearity( f ) ) / 2;
        wt = weightBoole( f );
        a[ wt ]++;
        int  cpt=1, limite = 1 << c.nbl;
                for( int x = 0; x < ffsize; x++ )
                        t[x] = f[x];
        while ( wt >= goal && cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= c.fct[i][x];

                wt = weightBoole( t );
                a[ wt ]++;
                cpt++;
        }

        if ( cpt == limite ) {
                printf("\nNL2 ( %d )  :", goal);
                for( int i = 0; i <=ffsize; i++ )
                        if ( a[i] ) printf(" %d [ %d ]", a[ i ], i );
        }  ; // else printf("\ngoal : %d /%d\n", wt, goal );
        free( a );
        free( t );
        freecode( c );
}


void check(boole L, vector v)
{
   
    boole R = vectortoboole  ( v , & base );       

    //panf( stdout, L );
    //panf( stdout, R );
    initboole( 6 );
    code Q = rmcode(0, 2, 6);
    boole  t = getboole( );
    for(int  x = 0; x < 32; x++ ){
	       t[x] = L [x];
	    t[x+32] = R [x];
    }
    panf(stdout, t );

    int cpt = 1, limite = 1 << Q.nbl;
    int *a = calloc( ffsize+1, sizeof(int ) );
    int wt = (ffsize - linearity( t )) / 2;
    NLL( t, 10 );
    free( R );
    initboole( 5 );
    free(a);
}

void doit(boole f, boole g)
{
    int num = base.table[booleVector(g, &base)];
    printf("\nclasse=%d %d\n", base.table[booleVector(g, &base)], num );

    for (int v = 0; v < base.size; v++)
	if (  1 ||  base.table[v] == num) {
	    check( f, v  );
	}
}

int main(int argc, char *argv[])
{
    FILE *src = NULL;
    int opt;
    initboole( 6 );
    boole h = strtoboole( "anf=ace+bce+bde+bcf+adf");
    panf( stdout, h );
    NLL( h, 18 );

    initboole( 5 );
    boole L = getboole(),  R = getboole();;
    for( int x = 0; x < 32; x++ ){
	    L[x] = h[x];
	    R[x] = h[x+32];
    }
    
    panf( stdout, L );
    panf( stdout, R );

    while ((opt = getopt(argc, argv, "a:m:f:r:hiv:s:")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
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


    base = monomialBasis(3, 4, 5);
    initagldim(5);
    orbitBasic(mkaglGroup(), &base);


    printf("\nclass=%d\n", base.table[booleVector(L, &base)] );
    printf("\nclass=%d\n", base.table[booleVector(R, &base)] );
    vector v = booleVector(R, &base ); 
    check( L, v ); 

    boole f, t[2];
    int i = 0;
    while ((f = loadBoole(src))) {
	panf(stdout, f);
	int num = base.table[booleVector(f, &base)];
	int size = 0;
	for (int v = 0; v < base.size; v++)
	    if (base.table[v] == num)
		size++;
	printf("\nclass=%d size=%d\n", num, size);
	t[i++] = f;
    }
    fclose(src);




    return 0;
}
