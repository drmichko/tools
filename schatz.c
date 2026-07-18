
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"
#include "distrib.h"

#define TARGET 88
aglGroup grp ;
int64_t size ;


 
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



int *table( boole f )
{ code cc = rmcode( 2, 2, 7 );
  int * A = calloc( 1 << cc.nbl , sizeof(int ) );
	
	int  cpt=1, limite = 1 << cc.nbl;
        boole  t = getboolecpy( f );
        A[0] = ( ffsize - linearity( f ) ) / 2;
        int q = 0;
        while ( cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= cc.fct[i][x];
                q ^=  ( 1 << i );
                A[ q  ] = ( ffsize - linearity( t ) ) /2;
                cpt++;
        }

        free( t );
        freecode( cc );
      	//pdistrib( "\n W=", A, limite  ); 
  return A;
}

int best = TARGET;

int test( int* F, int * G)
{
	int limite = rmdimen( 2, 2, 7 );
	limite = 1 << limite;
	for( int g = 0; g < limite; g++ ) {
		int q = 0;
		while ( (F[ q ]  + G[ q ^ g]  >  TARGET  ) &&  ( q < limite ) ) q++; 
		if ( q  == limite ) {
			puts("YES");
		        return 1;
		}
		int tmp = F[ q ]  + G[ q ^ g] ;
		if ( tmp > best ) {
			printf("best=%d\n", tmp );
			best = tmp;
		}
	}
        return 0;
}

int main(int argc, char *argv[])
{  

    char *fn = NULL;


    aglGroup g = NULL;
    int num = 0, cls = -1;
    int opt, optw = 0, optr=-1, optinit = 0;
    int deg = 0, dimen = 7;
    int job = 0, mode = 1;

    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:j:")) != -1) {
	switch (opt) {
	case 'w':
	    optw++;
	    break;
	case 'd':
	    deg = atoi(optarg);
	   break;
	case 'c':
	    cls = atoi(optarg);
	    break;
	case 'i': 
		optinit = atoi( optarg);
	break;
	case 'j': 
	        sscanf( optarg, "%d:%d", &job, & mode );
		break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    optr=atoi( optarg );
	    break;
	default:
	    exit(0);
	}

    }


    initboole(  dimen );
    initagldim( dimen );

    boole f;
    num = 0;
    uint64_t grpSize, orbSize;
    int val = 0;
    
    FILE * src = fopen(  fn , "r" );
    if ( ! src ) {
		perror( fn  );
		exit(1);
	}
    assert( mode > 0 );
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
		 panf( stdout, f );
		 int * F = table( f );
	    	 boole g;
	    	 int nb = 0;
		 boole c;
    	    	 while (( c = loadBoole(src ) ) ) {
			 if ( nb % mode == job ) {
				 boole g = getboole(  );
				 for( int x = 0; x < ffsize; x++ )
					 g[x] ^= f[x] ^ c[x] ;
				 int * G = table( g );
				 if ( test( F, G ) ) {
					 FILE *dst = fopen( "goodies.out", "a");
					 fprintf( dst, "%s\n", fn );
					 initboole( 8 );
					 boole h = getboole();

					 for( int x = 0; x < ffsize/2; x++ )
						 h[x] = f[x];

					 for( int x = 0; x < ffsize/2; x++ )
						 h[x+ffsize/2] = f[x] ^ c[x];
					 panf( dst, h );
					 fclose(dst);
					 // NLL( h, 90  );
					 initboole( 7 );
					 free( g );
				 }
				 free( G );
			 }
			 free( c );
			 nb++;
		 }
		 printf("\n#nb: %d\n", nb);
		 free(f);
            aglfreeGroup( grp );
	    num++;
        }

     fclose(src);

     printf("\n#maps : %d\n", num );
     return 0;
}
