/*
 * nq7_test.c : test NQ(f) >= T pour une fonction booléenne à 7 variables.
 * Table de vérité : f(x) = (tt[x >> 6] >> (x & 63)) & 1.
 *
 * Compilation : gcc -O3 -march=native -o nq7_test nq7_test.c
 */
#include <stdint.h>
#include <stdio.h>
    

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"


#define M      7
#define N      (1 << M)
#define NMONO  (M + M * (M - 1) / 2)   /* 28 monômes de degré 1 ou 2 */

static uint64_t mono[NMONO][2];
static int mono_ready = 0;

static void build_monomials(void)
{
    int k = 0;
    for (int i = 0; i < M; i++)
        for (int j = i; j < M; j++, k++) {      /* j == i : monôme x_i */
            mono[k][0] = mono[k][1] = 0;
            for (int x = 0; x < N; x++)
                if (((x >> i) & (x >> j)) & 1)
                    mono[k][x >> 6] |= 1ULL << (x & 63);
        }
    mono_ready = 1;
}

/*
 * Retourne 1 si NQ(f) >= T, 0 sinon.
 * Échec dès qu'une g de degré <= 2 vérifie d(f,g) < T ou d(f,g) > N - T
 * (cas g + 1, par complément).
 */
int nq_at_least_m7(const uint64_t tt[2], int T)
{
    if (!mono_ready) build_monomials();

    const int dmin = T;        /* distance admissible : T <= d <= N - T */
    const int dmax = N - T;

    uint64_t c0 = tt[0], c1 = tt[1];           /* f + g courant */

    for (uint32_t i = 1; ; i++) {
        int d = __builtin_popcountll(c0) + __builtin_popcountll(c1);
        if (d < dmin || d > dmax) return 0;    /* sortie anticipée */

        if (i >> NMONO) return 1;              /* 2^28 termes parcourus */
        const uint64_t *m = mono[__builtin_ctz(i)];
        c0 ^= m[0];
        c1 ^= m[1];
    }
}



uchar *loadZipname(int target , basis_t * b)
{
    FILE *srcFile;
    size_t lus;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;
    size = size / 8;
    uchar *res;
    sprintf(fn, "data/%d-%d-%d-%d.zip", target, b->s, b->t, b->m);
    srcFile = fopen(fn, "r");

    res = calloc(size, sizeof(uchar));
    if (!res) {
        perror("allocation table");
        exit(1);
    }

    if (!srcFile) {
        perror(fn);
        exit(1);
    }
    lus = fread(res, 1, size, srcFile);
    if (lus != size) {
        printf("load problem:-(\n");
        printf("lus=%ld spaceSize=%ld\n", lus, size);
        exit(1);
    }
    fclose(srcFile);
    return res;
}

#define _POSIX_C_SOURCE 200809L
int value(FILE *src, uint64_t *v )
{
    char  *line = NULL;
    size_t cap  = 0;
    while (getline(&line, &cap, src) != -1)
        for (const char *p = line; *p; p++)
            if (sscanf(p, "#orbits=%ld", v) == 1) {
                free(line);
                return  1;
            }

    free(line);
    return 0;
} 

int main(int argc, char *argv[])
{  

    int opt; 
    int dimen = 7;
    int target = 0;
    int goal = 38;
    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:g:")) != -1) {
	switch (opt) {
	case 'w':
	    break;
	case 'd':
	   break;
	case 'c':
	    break;
	case 'i': 
	break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
                break;
	case 'g':
	    goal  = atoi( optarg ); 
	    break;
	case 't':
	    target=atoi( optarg );
	    break;
	default:
	    exit(0);
	}

    }

    if ( ! target ) return 0;

    char name[128];
    sprintf( name, "rho/rho-%d.out", target );


    FILE *src = fopen( name, "r" );
    initboole(  dimen );
    initagldim( dimen );

    int item = 0;
    uint64_t size = 0;
    basis_t base = monomialBasis( 3, 3, 7 );

    boole h  = loadBoole( src );
    uint64_t  orbits;
    if ( ! value( src, &orbits  ) ) exit(1);
    fclose(src);

    clock_t debut = clock();
    printf("\nUsing no %d\n", target );fflush(stdout);
    size_t count = 0, soluce = 0;
    uchar * zip = loadZipname( target, &base );
    for( vector v = 0; v < base.size   &&   count < orbits ; v++ ) 
	    if ( 1  & ( zip[ v / 8 ] >> ( v % 8 ) ) )  {
		    count++;
		    boole  f = vectortoboole( v, &base );
		    uint64_t tt[2]={0,0};
		    for ( int x = 0; x < ffsize; x++ ){
			f[x] ^= h[x];
			if ( f[x]  ) tt[ x / 64] ^= ((uint64_t) 1 )  << ( x & 63 );
		    }
                    if ( nq_at_least_m7( tt , goal ) ) {
			    soluce++;
			    panf( stdout, f );
		    	    fflush( stdout);
		    }
		    free( f );
	    } 
    printf("\n#orbit=%ld count=%ld  %ld\n", orbits, count , soluce);
    clock_t fin = clock();
    double temps = (double)(fin - debut) / CLOCKS_PER_SEC;
    printf("\n#runtime : %.6f secondes\n", temps);
    return 0;
}
