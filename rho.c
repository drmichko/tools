
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"

aglGroup grp ;
int64_t size ;

void doit( boole f )
{
orbitData df  = initData( f , 4, grp, size );
freeData( &df );
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

int main(int argc, char *argv[])
{  

    char *fn = NULL;


    aglGroup g = NULL;
    int num = 0, cls = -1;
    int opt, optw = 0, optr=0, optinit = 0;
    int deg = 0, dimen = 6;
    int target = 3;

    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:")) != -1) {
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
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    optr=atoi( optarg );
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
    num = 0;
    uint64_t grpSize, orbSize;

    basis_t base   = monomialBasis( 5 , ffdimen ,  ffdimen);
    int64_t aglSize = aglcard(ffdimen); 
    aglVectorGroup  ldg = aglVectorGroupAction( mkaglGroup(6),  & base );
    
    while ((f = loadBoole(src ))) {
            panf( stdout, f );
	    vector vec  = booleVector( f, & base );
            initBrowse( &base );
            size_t orbSize = browse( vec , ldg  );
            printf("\norbsize=%ld", orbSize );
            assert( 0 == aglSize %  orbSize );
            size_t stabSize  = aglSize / orbSize;
            aglGroup stab = boundStabilizer( vec ,  f, grp, & base, stabSize);
	    paglGroup( stdout, stab );
	    fprintf( stdout, "\nstabSize=%ld\n", stabSize ); 
            free( base.table);
	    aglVectorGroupFree( ldg );
            aglfreeGroup( stab );
	    free( f );
	    num++;
        }

     fclose(src);

     printf("\n#maps : %d\n", num );
     return 0;
}
