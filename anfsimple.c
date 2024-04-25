
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "boolean.h"

int count = 0;
int score = 0;
int verbe = 0;

char * trace;


size_t hash( agl A )
{
size_t res = 0;
int i;
for(i=0; i <= agldimen; i++ )
	res = (res << agldimen) + A[i];
return res;
}



void randAction( boole f )
{ boole g = getboole();
	int x;
	agl A = getaglrand( random() % ffsize ); 
	for( x = 0; x < ffsize; x++ )
		g[x]= f[ aglImage( x, A ) ];

	for( x = 0; x < ffsize; x++ )
		f[x] = g[x];
	free(g);
	free(A);
}

void alphause( boole f )
{  int j, u;
   int cpt[26]={0};
	   ANFtoTT(f);
	   for( u = 0; u < ffsize; u++ )
		   if ( f[u] )
			   for( j = 0; j < ffdimen; j++ )
				   if ( (1<<j) & u ) cpt[j]++;
	   for( j = 0; j < ffdimen; j++ )
		   printf("\n%c : %d ", 'a'+j, cpt[j] );
	   ANFtoTT(f);
}


int doit( boole f , int val, size_t limite )
{
	int nbt = 1024;
	int u, wt;
	boole g = NULL;
	size_t iter = 0;
	while ( iter  < limite ) {
	   ANFtoTT(f);
	   wt = 0; 
	   for( u = 0; u < ffsize; u++ ){
		   if ( weight(u) < val ) f[u] = 0;
		   wt += f[u];
	   }
	   ANFtoTT(f);
	   if ( wt < nbt ) {
		   nbt = wt;
		   if ( verbe ) {
	   	   	   printf("\n");
			   panf( stdout, f );
		   	   printf("#iter=%lu nbt=%d\n", iter, nbt );
		   }
		   if ( verbe > 1 ) alphause( f );
		   if ( g ) free(g);
		    g = getboolecpy( f );
	   }
	   iter++;
	   randAction( f );
	}
	panf( stdout, g );
	printf("\n#nbt=%d", nbt);
	free( g );
	free( f );
	return nbt;
}


int main( int argc, char* argv[] )
{       int opt;
	int val = 0;
	long long limite = 1;
	char *anf = NULL;
	FILE *src = NULL;
        while ((opt = getopt(argc, argv, "a:m:f:r:hiv:s:")) != -1) {
               switch (opt) {
               case 'm':
                   initboole( atoi( optarg ) );
                   break;
               case 'r':
                   limite = limite << atoi( optarg );
                   break;
               case 'a':
                   anf = strdup(optarg);
                   break;
	       case 'f':
                   src = fopen( optarg, "r" );
		   if ( ! src ) {
			   perror( optarg );
			   return 1;
		   }
                   break;

               case 'v':
                   verbe++;
                   break;
               case 's':
                   val = atoi( optarg );
                   break;
               default: /* '?' */
                   fprintf(stderr, "Usage: %s [-a anf ] [-r log iter]\n",
                           argv[0]);
                   exit(EXIT_FAILURE);
               }
           }

	boole f;;
	srandom( time(NULL) + getpid() );
	score=1024;
	if ( anf ) {
		 int  nbv = anfdimen( anf  );
		 printf("\ndimension=%d\n",  nbv );
		 initboole( nbv );
		 initagldim( ffdimen );
        	 f = strtoboole( anf  );
		 if ( ! val  )
		   val = valuation( f );
		 score = doit( f, val , limite );
	}

	initagldim( ffdimen );
	if ( src ) {
		while ((f = loadBoole(src))) {
		      int r = doit( f, val, limite);
		      if ( r < score ) score = r;

		}
		fclose( src );
	}
	printf("\nscore=%d\n", score);
	return 0;
}
