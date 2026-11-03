#include "boolean.h"
#include "unistd.h"
#include "string.h"
#include "sstools.h"
#include <time.h>
#include "distrib.h"

int check( boole h, aglGroup S, int s )
{
while ( S ) {
	boole f = getboole();
	for( shortvec x = 0; x < ffsize; x++ )
		f[x]  = h [  aglImage( x, S->per  ) ];
	for( shortvec x = 0; x < ffsize; x++ )
		f[x] = f[x] ^ h[x];
	assert( degree(f) <  s ) ;
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

void help( void )
{
system( "sed -n '/^[a-zA-Z_].*\\boption\\s*(/,/^}/p' stabredo.c");
}

    int dimen = 8;
    int opt;
    int tour = 1024;
    int verbe = 0;
    int optgen  = 1024;  
    int optw = 0, item =0; 
    int optS = 0;
    char *fn = NULL;    
    char* destin = NULL;

int option(int argc, char *argv[])
{

    while ((opt = getopt(argc, argv, "m:s:vt:rs:f:g:w:hW")) != -1) {
	switch (opt) {
	case 'f':
	    fn  = strdup(optarg);
	    break;
	case 'm':
	    dimen = atoi(optarg);
	    break;
	case 't':
	    tour = atoi(optarg);
	    break;
	case 'v':
	    verbe = 1;
	    break;
	case 'r' :
		srandom( time(NULL ));
		break;
	case 'w':
	    destin = strdup( optarg );
	    break;
	case 'W':
	    optw++;
	    destin = "/tmp/stab.out";
	    break;
	case 'g':
		optgen = atoi( optarg );
	        break;
	case 's':
	    optS = atoi(optarg);
	    break;
	case 'h' : help();
		   return 0;
	default:		/* '?' */
	    fprintf(stderr, "check usage of %s\n", argv[0]);
	    exit(EXIT_FAILURE);
	}
    }
    return 1;
}


int main(int argc, char *argv[])
{

    if ( ! option( argc, argv ) ) 
		exit(0);
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
    int max = 0;
    FILE* dst = stdout;

    if (  destin ) {
		dst = fopen( destin , "w");
		if ( ! dst  ) {
		perror( destin );
		return 0;
          }
    }

    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
    	panf( dst, f);
	int nb  = nbgenerator( grp );
	if ( optS )
		check( f, grp, optS ) ;
	if ( nb > optgen  ) {
        	aglGroup   S = aglReduce( grp, tour  );
    		paglGroup( dst , S );
    		assert(   ssGroupOrder( S ) == grpSize );
		int  n = nbgenerator( S ) ;
		fprintf( dst, "\n#generator: %d --> %d\n", nb, n );
		if ( n < nb ) count++;
                nb = n;
		
	} else {
		paglGroup( dst , grp );
	}
	if ( nb > max ) max =  nb;
    	fprintf(dst , "\nstabSize=%ld", grpSize );
	assert( f != NULL );
	free( f );
	item++;
    }
    fclose(src);
    if  (destin ) 
		fclose(dst);
    fprintf(stderr,"%d class in %s", item, destin );
    fprintf(stderr, "\n#change=%d  max=%d\n", count, max );
    if ( destin   ) {
	printf("new file : %s\n", destin );
	if ( optw == 1 )
		puts("NO UPDATE : use -WW to updtate!" );
    }
    if ( optw > 1 ) {
	char cmd[128];
	sprintf(cmd, "mv %s %s",  destin, fn  );
	puts(cmd);
	system( cmd );
    }
    return 0;

}
