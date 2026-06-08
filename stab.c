
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

int main(int argc, char *argv[])
{  

    char *fn = NULL;


    aglGroup g = NULL;
    int num = 0, cls = -1;
    int opt, optw = 0, optr=0, optinit = 0;
    int deg = 0;
    while ((opt = getopt(argc, argv, "if:c:wb:d:m:r:")) != -1) {
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
		optinit = 1;
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


    initboole(7);
    initagldim(7);

    boole f;
    num = 0;
    uint64_t grpSize, orbSize;
    int val = 0;
    if ( optinit  ) {
    	FILE *src = fopen( fn , "r");
    	if (  ! src  ) {
		perror( fn );
		exit(1);
    	}

       sprintf(fn,  "stab/stab-%d.txt", optr );
       FILE * dst = fopen( fn, "w" ); 
	grp = mkaglGroup( ffdimen );
	int64_t size = aglcard(ffdimen);
    	while ((f = loadBooleValue(src, &val))) {
            panf( dst, f );
	    paglGroup( dst, grp );
	    fprintf(dst, "\nstabSize=%ld\n", size );
	    num++;
	    free(f);
	}
        fclose(dst);
        fclose(src);
       printf("\n#maps : %d\n", num );
    }
    if ( ! optr ) exit(0);

    char tmp[64];
    sprintf( tmp ,  "stab/stab-%d.txt", optr );
    FILE * src = fopen(  tmp , "r" );

    sprintf( tmp ,  "stab/stab-%d.txt", optr - 1 );
    FILE * dst = fopen( tmp , "w" );
        
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
            panf( stdout, f );
            paglGroup( stdout, grp );
            fprintf(stdout, "\nstabSize=%ld\n", size ); 
            basis_t base   = monomialBasis( optr, optr,  ffdimen);
            vector vec  = booleVector( f, & base );
            aglVectorGroup ldg  = aglBoundaryGroupAction( f , grp , &base );
            initBrowse( &base );
            size_t orbsize = browse( 0, ldg  );
            printf("\norbsize=%ld", orbsize );
            assert(   grpSize % orbsize == 0 );
            free( base.table);
            size_t stabSize = grpSize /orbSize ;
            printf("\nstabsize=%ld", stabSize );

            //res.W[r].route = transversale( vec ,  res.W[r].ldg, & res.W[r].base  );
	    //total *= res.W[r].base.orbSize;
	    aglVectorGroupFree( ldg );
            free(f);
            aglfreeGroup( grp );
        }

     fclose(dst);
     fclose(src);

     printf("\n#maps : %d\n", num );
     return 0;
}
