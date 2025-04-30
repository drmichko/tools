#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"

boole test(int k)
{
    boole f = getboole();
    for (int u = 0; u < ffsize; u++)
	if (weight(u) == k)
	    f[u] = 1;
    ANFtoTT(f);
    return f;
}

void RD(int dr, boole f, int iter)
{
    int t = dr;

    shortvec trans, base[ffdimen / 2];

    int i, r = 0;
    shortvec p;
    int limite = 1 << t;

    int q = ffsize;
    int old = ffdimen;
    panf(stdout, f);
    initboole(t);
    boole g = getboole();
    int best = t;
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
	int d = degree(g);
	if (d < best) {
	    best = d;
	}

    }
    initboole(old);
    printf("\niter=%d score=%d\n", cpt, best);
}

int main(int argc, char *argv[])
{
    FILE *src = NULL;
    char *anf = NULL;
    int opt, r = 0;

    while ((opt = getopt(argc, argv, "a:m:f:r:hiv:s:")) != -1) {
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

    if (r == 0)
	r = (ffdimen + 1) / 2;
    printf(" r= %d", r);


    if (anf) {
	f = strtoboole(anf);
	RD(r, f, 100000);
    }

    if (src) {
	while ((f = loadBoole(src))) {
	    RD(r, f, 100000);
	}
	fclose(src);
    }


    return 0;
}
