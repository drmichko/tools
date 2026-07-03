
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
    int deg = 0, dimen = 7;
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

    optr = 4;

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

    basis_t base   = monomialBasis( optr, optr,  ffdimen );
    size_t  taille[ 2000 ] = {0};
    int  no = 1;
    initBrowse( &base );
    aglVectorGroup  ldg = aglVectorGroupAction( grp , & base );
    while ((f = loadaglboolesize(src, &grp, &grpSize))) {
	    vector vec;
	    boole g = getboolecpy( f );
	    projboole( 4, 4, g );
	    vec  = booleVector( g, & base );
	    int cls =  base.table[ vec ];
	    size_t orbSize;
	    panf( stdout, f );
	    panf( stdout, g );
            printf("check:%d", checkstab( g, grp , optr ) );	    
            free( g );

            if ( cls == 0  ) {
		    orbSize = numbrowse( no, vec , ldg  );
	    	    fprintf(stdout, "\norbSize=%ld\n", orbSize); 
            	    taille[ no ] =  orbSize;
		    cls = no;
		    no++;
	    }
	    printf("\n#cls=%d / %ld ", cls, vec );
            orbSize = taille[ cls ];
            assert(   grpSize % orbSize == 0 );
            size_t stabSize = grpSize /orbSize ;
	    assert(   grpSize % orbSize == 0 );
	    fprintf(stdout, "\norbSize=%ld\n", orbSize); 
	    fprintf(stdout, "\nstabSize=%ld\n", stabSize ); 
            aglGroup stab = NULL;
            stab = plainStabilizer( vec , grp, & base, stabSize);

            assert( checkstab( f, stab, optr - 1 ) == 1 );	    
            checkstab( f, stab, optr - 1 );	    
	    panf( dst, f );
	    paglGroup( dst, stab );
	   
	    fprintf(dst, "\nstabSize=%ld\n", stabSize ); 
            aglfreeGroup( stab );
            free( f );
            aglfreeGroup( grp );
	    num++;
        }

     fclose(dst);
     fclose(src);

     printf("\n#maps : %d\n", num );
     printf("\n#orbs : %d\n", no );
     return 0;
}
