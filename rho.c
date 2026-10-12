
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
code rmh;


typedef struct info {
      boole f;
      int*  a;
      int deg;
      uint64_t stab;
      aglGroup grp;
} info;

info table[ 64 ];

void NL2( int item, uint64_t size , boole f )
{
        boole  t = getboolecpy( f );
        int *  a = calloc( ffsize+1, sizeof( int ) );
        int wt = ( ffsize - linearity( f ) ) / 2;
        a[ wt ]++;
        int nl2 = wt;
        int  cpt=1, limite = 1 << rmh.nbl;
        while ( cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= rmh.fct[i][x];

                wt = ( ffsize - linearity( t ) ) / 2;;
                a[ wt ]++;
		if ( wt < nl2) nl2 = wt;
                cpt++;
        }

        printf("item=%2d deg=%d  stab=%ld NL2=%2d :", item, degree( f  ), size, nl2 );
        for( int i = 0; i <=ffsize; i++ )
                        if ( a[i] ) printf(" %d [ %d ]", a[ i ], i );
                putchar('\n');
        table[item].a   = a;
	table[item].f   = f;
	table[item].stab   = size;
	table[item].deg = degree( f );;
        free( t );
}

void doit( int item, uint64_t size ,  boole f )
{
		NL2( item, size, f );
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

int accept( int i, int j )
{
if ( table[i].deg <=3 && table[j].deg <= 3) return 0;
int L = 0, R = 0; 
for( int r = 0; r < 40 ; r++ ){
	L=0;
        R=0;
        for( int u=0; u < r; u++ )
		L += table[i].a [ u ];
        for( int u=64; u >= 40 - r; u-- )
		R += table[j].a [ u ];
	if ( L > R ) return 0;
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
    rmh = rmcode(2,2,6);
   int item = 0;
   uint64_t size = 0;
    while ((f = loadaglboolesize( src , &g, &size ))) {
	if ( degree(f) <=4 ) {
		doit(item,  size, f );
		table[item].grp  = g;
		item++;
	}
    }
    printf("\nitem=%d\n", item );
    int count = 0;
    for( int i =0; i < item; i++ )
    for( int j =i ; j < item; j++ )
	if ( accept(i,j) )  {
		count++; 
    		printf("i=%d j=%d   (%ld,%ld) \n", i, j, table[i].stab, table[j].stab  );
	}
    printf("\npair=%d\n", count);
    return 0;

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
