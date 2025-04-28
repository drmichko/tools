#include <stdlib.h>
#include <stdio.h>
#include "boolean.h"
#include <assert.h>

void distrib(int *F)
{
    int x;
    int *mult = calloc(2 * ffsize, sizeof(int));
    for (x = 0; x < ffsize; x++)
	mult[F[x] + ffsize]++;
    for (x = 0; x < 2 * ffsize; x++)
	if (mult[x])
	    printf(" %d [%d]", mult[x], x - ffsize);
    free(mult);

}

int limite = 0;

int test(boole f, int w)
{
    int x;
    int best = ffsize;
    int *F = calloc(ffsize, sizeof(int));
    for (x = 0; x < ffsize; x++)
	F[x] = f[x];
    Fourier(F, ffsize);
    int cpt = 0;
    for (x = 0; x < ffsize; x++)
	if (abs(F[x]) >= 4)
	    cpt++;
    if (cpt <= limite) {
	panf(stdout, f);
	distrib(F);
    }

    free(F);
    return 1;
}

int main(int argc, char *argv[])
{
    FILE *src;
    boole f;
    int num = 0;

    initboole(atoi(argv[1]));

    src = fopen(argv[2], "r");
    if (!src) {
	perror("");
	exit(1);
    }
    int w = atoi(argv[3]);
    int theta = 8 - ffdimen;
    limite = 1 << (ffdimen - 2 * theta);
    while ((f = loadBoole(src))) {
	int x, sum = 0;
	for (x = 0; x < ffsize; x++)
	    sum += f[x];
	if (sum == w)
	    test(f, w);
	if (sum == ffsize - w) {
	    for (x = 0; x < ffsize; x++)
		f[x] = 1 - f[x];
	    test(f, w);
	}
	num++;
	free(f);
    }
    fclose(src);
    printf("\n#%d Boolean functions", num);
    return 0;
}
