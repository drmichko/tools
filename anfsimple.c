
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "boolean.h"

int count = 0;
int score = 0;

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
int main( int argc, char* argv[] )
{       int opt;
	int val = 0;
	long long iter = 0 , limite = 10;
	char *anf = "abc";;
        while ((opt = getopt(argc, argv, "a:r:hV:")) != -1) {
               switch (opt) {
               case 'r':
                   limite = 1 << atoi( optarg );
                   break;
               case 'a':
                   anf = strdup(optarg);
                   break;
               case 'V':
                   val = atoi( optarg );
                   break;
               default: /* '?' */
                   fprintf(stderr, "Usage: %s [-a anf ] [-r log iter]\n",
                           argv[0]);
                   exit(EXIT_FAILURE);
               }
           }

	boole f;
	int nbv = anfdimen( anf  );
	int nbt = 1024;
	int u, wt;


	srandom( time(NULL));
	printf("\ndimension=%d\n",  nbv );
	initboole( nbv );
	initagldim( ffdimen );
        f = strtoboole( anf  );
	panf(stdout, f );
	printf("\n");
	pwalsh( f);
	printf("\n");
	if ( ! val  )
		val = valuation( f );

	while ( iter  < limite ) {
	   ANFtoTT(f);
	   wt = 0; 
	   for( u = 0; u < ffsize; u++ ){
		   if ( weight(u) < val ) f[u] = 0;
		   wt += f[u];
	   }
	   ANFtoTT(f);
	   if ( wt <  nbt ) {
	   	   printf("\n");
	   	   panf( stdout, f );
		   alphause( f );
		   nbt = wt;
		   printf("\niter=%Ld nbt=%d\n", iter, nbt );
	   }
	   iter++;
	   randAction( f );
	}


	return 0;
}
