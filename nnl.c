#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"
#include "math.h"

boole test(int k)
{
    boole f = getboole();
    for (int u = 0; u < ffsize; u++)
	if (weight(u) == k)
	    f[u] = 1;
    ANFtoTT(f);
    return f;
}

int  RD(int dr, boole f, int iter)
{
    int t = dr;

    shortvec trans, base[ffdimen / 2];

    int i, r = 0;
    shortvec p;
    int limite = 1 << t;

    int q = ffsize;
    int old = ffdimen;
    initboole(t);
    boole g = getboole();
    int best = 0;
    shortvec x, y;
    int cpt = 0;
    while (iter--) {
	r = 0;
	cpt++;
	while (r < t) {
	    y = random() % q;
	    for (i = 0; i < r; i++) {
		p = base[i] & (base[i] - 1);
		p ^= base[i];
		if (p & y)
		    y ^= base[i];
	    }
	    if (y != 0)
		base[r++] = y;
	}
	trans = random() % q;
	for (x = 0; x < limite; x++) {
	    y = 0;
	    for (i = 0; i < r; i++)
		if ((x & (1 << i)))
		    y ^= base[i];
	    g[x] = f[y ^ trans];
	}
	int d = linearity( g );
	if (d >  best) {
	    best = d;
	}

    }
    initboole(old);
    return  best;
}

boole symboole( int  d )
{ boole r = getboole();
  for( int u = 0; u < ffsize; u++ )
	  r[u] = weight(u) == d ? 1 : 0;
  ANFtoTT( r );
  return ( r );
}

int main(int argc, char *argv[])
{
    FILE *src = NULL;
    char *anf = NULL;
    int opt, r = 0;
    int iter=1000;
    int special= 0;
    while ((opt = getopt(argc, argv, "a:m:f:r:hi:vs")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'a':
	    anf = strdup(optarg);
	    break;
	case 'r':
	    r = atoi(optarg);
	    break;
	case 'i':
	    iter = atoi(optarg);
	    break;
	    case 's' : special++; break;
	case 'f':
	    src = fopen(optarg, "r");
	    if (!src) {
		perror(optarg);
		return 1;
	    }
	    break;

	}
    }


    boole f;

    if ( special ) {
	    f = symboole( 4 );
	    panf(stdout, f);
    	    printf("\niter=%d deg=%d lin=%d :", iter, degree(f) , linearity(f) );
	    float nnl = 0.5 * ( 1 - (float) linearity( f )  / ffsize  ) ;
	    printf(" %.4f",  nnl );
	    for( int i = ffdimen -1; i >=r; i-- ){
		    int lin = RD(r, f, iter);
		    float nnl = 0.5 * ( 1 - (float) lin / (1<<i) );
		    printf(" %.4f",  nnl );
	    }
    }
    if (r == 0)
	r = (ffdimen + 1) / 2;
    printf(" r= %d", r);


    if (anf) {
	f = strtoboole(anf);
	RD(r, f, iter);
    }

    if (src) {
	while ((f = loadBoole(src))) {
	    panf(stdout, f);
    	    printf("\niter=%d deg=%d lin=%d :", iter, degree(f) , linearity(f) );
	    float nnl = ( 1 - (float) linearity( f )  / ffsize  ) / 2;
	    printf(" %.4f",  nnl );
	    int old = linearity(f);
	    for( int i = ffdimen -1; i >=r; i-- ){
		    int lin = RD(r, f, iter);
		    float nnl = ( 1 - (float) lin / (1<<i) ) / 2;
		    printf(" %.4f",  nnl );
	    }
	}
	fclose(src);
    }


    return 0;
}
