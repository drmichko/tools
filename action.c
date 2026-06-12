
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


int checkstab( boole f, aglGroup g, int r )
{
while ( g ) {
    boole h = getboole( );
    for ( shortvec x = 0; x < ffsize; x++)
        h[x] = f[  aglImage(x, g->per ) ] ^ f[x];
    if( degree(h) >  r ) {
                printf("\ndegrees : f=%d  b=%d\n", degree(f),  degree(h) );
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
    int deg = 0, dimen = 7;
    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:")) != -1) {
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
		perror( tmp );
		exit(1);
	}
    
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
            panf( stdout, f );
	    int t = degree( f );
            paglGroup( stdout, grp );
            fprintf(stdout, "\nstabSize=%ld\n", grpSize ); 
            basis_t base   = monomialBasis( optr, optr,  ffdimen);
	    
            aglVectorGroup  ldg = NULL;
            ldg = aglBoundaryGroupAction( f , grp , & base );
               
            initBrowse( &base );
            int rank = 0;
            for( vec = 0; vec < base.size; vec++ )
		if ( ! base.table[vec] ) {
            		size_t orbSize = browse( vec , ldg  );
            		printf("\norbsize=%ld", orbSize );
			rank++;
                        }
	
            printf("\nrank=%d", rank );
            free(f);
            aglfreeGroup( grp );
	    num++;
        }

     fclose(src);

     printf("\n#maps : %d\n", num );
     return 0;
}
