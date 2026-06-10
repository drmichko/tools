
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
    if ( optinit  ) {
    	FILE *src = fopen( fn , "r");
    	if (  ! src  ) {
		perror( fn );
		exit(1);
    	}

       sprintf(fn,  "stab/stab-%d.txt", optinit );
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
        return 0;
    }
    if ( ! optr ) return (0);

    char tmp[64];
    sprintf( tmp ,  "stab/stab-%d.txt", optr );
    FILE * src = fopen(  tmp , "r" );
    if ( ! src ) {
		perror( tmp );
		exit(1);
	}
    sprintf( tmp ,  "stab/stab-%d.txt", optr - 1 );
    FILE * dst = fopen( tmp , "w" );
    
    if ( ! dst) {
		perror( tmp );
		exit(1);
	}
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
            panf( stdout, f );
	    int t = degree( f );
            paglGroup( stdout, grp );
            fprintf(stdout, "\nstabSize=%ld\n", grpSize ); 
            basis_t base   = monomialBasis( optr, optr,  ffdimen);
	    vector vec;
             aglVectorGroup  ldg = NULL;
            if ( t <= optr ) {
		    vec  = booleVector( f, & base );
               		ldg = aglVectorGroupAction( grp , & base );
	    }  else {
               vec = 0;
               ldg = aglBoundaryGroupAction( f , grp , & base );
               }
            initBrowse( &base );
            size_t orbSize = browse( vec , ldg  );
            printf("\norbsize=%ld", orbSize );
            assert(   grpSize % orbSize == 0 );
            size_t stabSize = grpSize /orbSize ;
            printf("\nstabsize=%ld", stabSize );
            aglGroup stab = NULL;
            stab = plainStabilizer( vec  , grp, &base, stabSize);
            assert( checkstab( f, stab, optr - 1 ) == 1 );	    
            checkstab( f, stab, optr - 1 );	    
	    panf( dst, f );
	    paglGroup( dst, stab );
            fprintf(dst, "\nstabSize=%ld\n", stabSize ); 
            free( base.table);
	    aglVectorGroupFree( ldg );
            free(f);
            aglfreeGroup( grp );
            aglfreeGroup( stab );
	    num++;
        }

     fclose(dst);
     fclose(src);

     printf("\n#maps : %d\n", num );
     return 0;
}
