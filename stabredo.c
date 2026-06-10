#include "boolean.h"
#include "unistd.h"
#include "string.h"
#include "sstools.h"
#include <time.h>
#include "distrib.h"



typedef struct  plain  {
        shortvec img[ 256 ];
        agl  A;
        struct plain * next;
} enrplaingrp, plaingrp;


void * rootp = NULL;


int countp = 0;



void orbittreeboole( boole f, boole h,  aglGroup S )
{ 
  
   while ( S  != NULL ) {
	   boole g  = getaglboole( f,  S->per  ) ;
	   //projboole(2, 4, g );
	   if ( newtt( g  , &rootp ) ) {
		countp++;
		orbittreeboole( g,  h,  S );
	    } else free(g );
            S  = S->next;
        }
}


int check( boole h, aglGroup S)
{
while ( S ) {
	boole f = getboole();
	for( shortvec x = 0; x < ffsize; x++ )
		f[x]  = h [  aglImage( x, S->per  ) ];
	for( shortvec x = 0; x < ffsize; x++ )
		f[x] = f[x] ^ h[x];
	assert( degree(f) <= 3 ) ;
	free(f );
	S = S->next;

}
return 0;
}
void free_boole(void *nodep)
{
    free(*(char **)nodep);
}


int nbgenerator( aglGroup g )
{
	int r = 0;
	while ( g ) {
		r++; 
		g = g->next;
	}
	return r;
}

int main(int argc, char *argv[])
{

    int dimen = 8;
    int opt;
    int tour = 1024;
    int optRnk = -1;
    int verbe = 0;
    int optSize =  0;
   
    char *fn = NULL;    
    while ((opt = getopt(argc, argv, "m:n:vt:r:s:f:")) != -1) {
	switch (opt) {
	case 'f':
	    fn  = strdup(optarg);
	    break;
	case 'm':
	    dimen = atoi(optarg);
	    break;
	case 's':
	    optSize = atoi(optarg);
	    break;
	case 't':
	    tour = atoi(optarg);
	    break;
	case 'v':
	    verbe = 1;
	    break;
	case 'r':
	case 'n':
	    optRnk = atoi(optarg);
	    break;
	default:		/* '?' */
	    fprintf(stderr, "Usage: %s [-a anf ] [-m diemn ]\n", argv[0]);
	    exit(EXIT_FAILURE);
	}
    }



    initboole(dimen);
    initagldim(dimen);




    FILE *src = fopen( fn , "r");
    if ( ! src  ) {
	    perror( fn);
	    exit(1);
    }
    
    boole f;
    uint64_t grpSize, orbSize;
    aglGroup grp;

    int val = 0;
    int num   = 0;
    int count = 0;
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
    	panf( stdout, f);
	int nb  = nbgenerator( grp );
	if ( nb > 2 ) {
        	aglGroup   S = aglReduce( grp, tour  );
    		paglGroup( stdout, S );
    		assert(   ssGroupOrder( S ) == grpSize );
		int  n = nbgenerator( S ) ;
		fprintf( stdout, "\n#generator: %d --> %d\n", nb, n );
		if ( n < nb ) count++;
		
	} else {
		paglGroup( stdout, grp );
	        fprintf( stdout, "\n#generator: %d\n", nb );
	}
    	printf("\nstabSize=%ld", grpSize );
	free( f );
    }
    fclose(src);
    printf("\n#change=%d\n", count );
    return 0;

}
