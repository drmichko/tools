
#define _XOPEN_SOURCE 500

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include "boolean.h"
#include "boolean.h"
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "mapping.h"
#include "code.h"
#include "stdint.h"
#include "math.h"

/* Boolean function are loaded sequentially
*  from a file, the degree and norme are printed.
*  dimension should be small!
*/
#define SIZE 256
int xvalue[32];
int rvalue[32];

int Xvalue[32];
int Rvalue[32];
int tfr[SIZE];
int cross[SIZE];
int linear;
int nonzero;
int stab = 0, optval = 0;
int optalpha = 0, oplin = 0, optx = 0, optr = 0, optyp = 0, optdeg =
    0, optbal = 0, opres = 0, opwt = 0, optnum = 0, opt2c = 0, opt2w=0;
int triphase = 0;
int optX = 0, optR = 0;
int optmod = 0;
int degmin, degmax, optsha = 0, shamin, shamax;
int opdiv = 0, divmin = 0, divmax = 0;
int linmin = 0, linmax = 0;
int optz = 0;
int zmax = 1;

typedef int64_t nombre;

nombre agsize, fixsize, total = 0;
float alphamin, alphamax;

int valcross[SIZE + 1] = { 0 };
int valtfr[SIZE + 1] = { 0 };

int optl = 2;
int valmin = 0, valmax = 0;

int balanced, nonzero;

int modulo(int x, int r)
{
    if (x >= 0)
	return x % r;
    x = (-x) % r;
    return (r - x) % r;
}

void test(boole f)
{
    int x, y, t, a, b;
    int *F = calloc(2 * ffsize, sizeof(int));
    for (x = 0; x < ffsize; x++)
	for (y = 0; y < ffsize; y++) {
	    a = x ^ y;
	    b = f[x] ^ f[y];
	    F[(a << 1) + b]++;
	}
    for (x = 0; x < 2 * ffsize; x++)
	for (y = x + 1; y < 2 * ffsize; y++)
	    if (F[x] > F[y]) {
		t = F[x];
		F[x] = F[y];
		F[y] = t;
	    }
    x = 0;
    printf("\ndis:");
    while (x < 2 * ffsize) {
	t = x;
	while (t < 2 * ffsize && F[x] == F[t])
	    t++;
	printf(" %2d [%d]", t - x, F[x]);
	x = t;
    }

}

void doit(boole f)
{
    galois x, y;
    for (x = 0; x < ffsize; x++)
	tfr[x] = 1 - 2 * f[x];

    for (x = 0; x < ffsize; x++)
	cross[x] = 0;
    for (x = 0; x < ffsize; x++)
	for (y = 0; y < ffsize; y++)
	    cross[x ^ y] += tfr[x] * tfr[y];

    Fourier(tfr, ffsize);
    nonzero = 0;
    for (x = 0; x < ffsize; x++)
	if (tfr[x])
	    nonzero++;
    balanced = nonzero < ffsize;
    linear = 0;
    for (x = 0; x < ffsize; x++) {
	if (abs(tfr[x]) > linear) {
	    linear = abs(tfr[x]);
	}
    }

}

boole myanfloadboole(FILE * src, int mode)
//-get the current boolean function
{
    boole res;
    char line[1024];
    while (!feof(src)) {
	line[0] = '\0';
	fgets(line, 1024, src);
	switch (*line) {
	case '#':
	    break;
	case 'a':
	case 'h':
	    res = strtoboole(line);
	    if (mode)
		return res;
	    break;
	case 'f':
	    sscanf(line, "fix=%ld", &fixsize);
	    return res;
	case 's':
	    sscanf(line, "stabSize=%ld", &fixsize);
	    return res;
	case 'c':
	    sscanf(line, "card=%ld", &fixsize);
	    return res;
	default:;
	}
    }
    return NULL;
}


int poweroftwow(boole f)
{
    galois x;
    for (x = 0; x < ffsize; x++) {
	if ((tfr[x] != 0)) {
	    int t = log(abs(tfr[x])) / log(2);
	    if (1 << t != abs(tfr[x]))
		return 0;
	}
    }
    return 1;
}

int poweroftwoc(boole f)
{
    galois x;
    for (x = 0; x < ffsize; x++) {
	if (( cross[ x ] != 0)) {
	    int t = log(abs( cross[x])) / log(2);
	    if (1 << t != abs(cross[x]))
		return 0;
	}
    }
    return 1;
}
int testphase(boole f)
{
    galois x = 0, y;
    while (tfr[x] == 0)
	x++;
    for (y = x; y < ffsize; y++)
	if (tfr[y] && (abs(tfr[y]) != abs(tfr[x])))
	    return 0;
    return 1;
}

int exclude(int t[], int k, int v[])
{
    galois x;
    for (int i = 0; i < k; i++) {
	for (x = 0; x < ffsize; x++)
	    if (abs(t[x]) == v[i])
		break;
	if (x < ffsize) {
	    return 0;
	}
    }
    return 1;
}

int allin(int t[], int k, int v[])
{
    galois x;
    for (x = 0; x < ffsize; x++) {
	int i;
	for (i = 0; i < k; i++)
	    if (abs(t[x]) == v[i])
		break;
	if (i == k) {
//              printf("<%d>", t[x] );
	    return 0;
	}
    }
    return 1;
}

double moment( int r , int n )
{
    int64_t sum = 0;
    for (int a = 0; a < ffsize; a++) {
	int64_t tmp = 1;
	for (int i = 0; i < r; i++)
	    tmp *= tfr[a];
	sum += tmp;
    }
    double res = sum;
    while ( n-- ) res /= ffsize;
    return res;
}

int sha(boole f)
{
    int res = 0, tmp;
    galois x, y, z, t;
    for (x = 0; x < ffsize; x++)
	for (y = x + 1; y < ffsize; y++)
	    for (z = y + 1; z < ffsize; z++) {
		t = x ^ y ^ z;
		tmp = f[x] ^ f[y] ^ f[z] ^ f[t];
		res += (1 - 2 * tmp);
	    }
    return res / 24;
}


int accept(boole f, int optnum, int num)
{
    int ok = 1;
    int wt, x;
    int tmp;
    if (stab)
	if (fixsize % stab)
	    return 0;

    if (optnum)
	return optnum == num;

    if (oplin) {
	int tmp = linear;
	ok = (linmin <= tmp) && (tmp <= linmax);

    }

    if (triphase && ok) {
	ok = testphase(f);
    }

    if (opwt && ok) {
	wt = 0;
	for (x = 0; x < ffsize; x++)
	    wt += f[x];
	if (wt != opwt)
	    ok = 0;
    }

    if (ok && optbal) {
	ok = balanced;
    }

    if (ok && optz) {
	ok = nonzero <= zmax;
    }
    if (ok && optalpha) {
	double alfa = moment(4, 3);
	ok = (alphamin <= alfa) && (alfa <= alphamax);
    }

    if (ok && opt2c) {
	ok = poweroftwoc(f);
    }
    if (ok && opt2w) {
	ok = poweroftwow(f);

    }
    if (ok && optx) {
	ok = exclude(tfr, optx, xvalue);

    }
    if (ok && optr) {
	ok = allin(tfr, optr, rvalue);

    }
    if (ok && optX) {
	ok = exclude(cross, optX, Xvalue);

    }
    if (ok && optR) {
	ok = allin(cross, optR, Rvalue);

    }
    if (ok && optdeg) {
	tmp = degree(f);
	ok = (degmin <= tmp) && (tmp <= degmax);
    }
    if (ok && optsha) {
	tmp = sha(f);
	ok = (shamin <= tmp) && (tmp <= shamax);
    }
    return ok;
}

int intcmp(const void *p1, const void *p2)
{
    return *((int *) p1) - *((int *) p2);
}

void absdistrib(char *msg, int *t, int n)
{
    int  f[ n ];
    for( int i= 0; i < n; i++ )
	    f[i] = abs( t[i] );

    printf("%s", msg);
    qsort(f, n, sizeof(int), intcmp);

    int i = 0, j;
    while (i < n) {
	j = i;
	while (j < n && f[i] == f[j])
	    j++;
	printf(" %d [ %d ]", j - i, f[i]);
	i = j;
    }
}
void distribmod(char *msg, int *t, int n, int r)
{
    int  f[ n ];
    for( int i= 0; i < n; i++ )
	    f[i] = ( t[i] < 0 ) ? (r - (- t[i] ) % r ) % r : t[i] % r;

    printf("%s/%d:", msg, r);
    qsort(f, n, sizeof(int), intcmp);

    int i = 0, j;
    while (i < n) {
	j = i;
	while (j < n && f[i] == f[j])
	    j++;
	printf(" %d [ %d ]", j - i, f[i]);
	i = j;
    }
}

void distribution(char *msg, int *f, int n)
{

    printf("%s", msg);
    qsort(f, n, sizeof(int), intcmp);

    int i = 0, j;
    while (i < n) {
	j = i;
	while (j < n && f[i] == f[j])
	    j++;
	printf(" %d [ %d ]", j - i, f[i]);
	i = j;
    }
}
void bigFourier(int64_t * f, unsigned int n)
// Transformation de Fourier sur place.
{
    int x, y;
    int64_t z;
    if (n > 1) {
	bigFourier(f, n >> 1);
	bigFourier(&(f[n >> 1]), n >> 1);
	for (x = 0, y = (n >> 1); x < (n >> 1); x++, y++) {
	    z = f[x];
	    f[x] = z + f[y];
	    f[y] = z - f[y];
	}
    };
}


void Vcor(boole f, int t)
{
    int *fct = calloc(ffsize, sizeof(int));
    int64_t *big = calloc(ffsize, sizeof(int64_t));
    galois x;
    for (x = 0; x < ffsize; x++)
	big[x] = f[x] ? -1 : +1;

    bigFourier(big, ffsize);
    for (x = 0; x < ffsize; x++) {
	int64_t tmp = 1;
	for (int j = 1; j <= t; j++)
	    tmp *= big[x];
	big[x] = tmp;
    }
    bigFourier(big, ffsize);

    for (x = 0; x < ffsize; x++)
	fct[x] = big[x] / ffsize;

    for (x = 0; x < ffsize; x++) {
	int tmp = fct[x];
	int v = 0;
	if (tmp < 0)
	    tmp = -tmp;
	if (tmp) {
	    while (0 == (tmp & 1)) {
		v++;
		tmp /= 2;
	    }
	} else
	    v = 1024;
	fct[x] = v;
    }
    fct[0] = 1024;
    distribution("val:", fct, ffsize);
    free(fct);
}



void correlation(boole f, int t)
{
    int64_t tmp;
    int *fct = calloc(ffsize, sizeof(int));
    int64_t *F = calloc(ffsize, sizeof(int64_t));
    galois x;
    for (x = 0; x < ffsize; x++)
	F[x] = f[x] ? -1 : +1;

    bigFourier(F, ffsize);
    for (x = 0; x < ffsize; x++) {
	tmp = 1;
	for (int j = 1; j <= t; j++)
	    tmp *= F[x];
	F[x] = tmp;
    }
    bigFourier(F, ffsize);
    for (x = 0; x < ffsize; x++)
	fct[x] = F[x] / ffsize;

    distribution("cor:", fct, ffsize);
    free(fct);
    free(F);
}



double alfa(boole f)
{
    long long int sum = 0, tmp;
    int *fct = calloc(ffsize, sizeof(int));
    galois x;
    for (x = 0; x < ffsize; x++)
	fct[x] = f[x] ? -1 : +1;

    Fourier(fct, ffsize);
    for (x = 0; x < ffsize; x++) {
	tmp = fct[x] * fct[x];
	sum += tmp * tmp;
    }
    free(fct);
    return (double) sum / ffsize / ffsize / ffsize;
}


void usage(char *str)
{
    puts("GENERAL:");
    puts("\td    :intervalle degree ");
    puts("\tm    :dimension");
    puts("\tf    :file of Boolean function");
    puts("SELECTION:");
    puts("\tb    :balanced");
    puts("\tz max:at most max non zero Walsh");
    puts("OUTPUT:");
    puts("\t%w    : Walsh distribution");
    puts("\t%wm8  : Walsh distribution modulo 8");
    puts("\t%c    : correlation distribution");
    puts("\t%c+   : correlation distribution absolute");
    puts("\t%d    : degree");
    puts("\t%l    : linearity");
    puts("\t%z    : number of non zero Walsh");
    puts("\t%n    : new line");
    puts("\t%s    : stabilizer size");
    puts("SAMPLES:");
    puts("./anfload.exe  -d3:4 -a1.75:2 '%d%a%x%n'");
    puts("./anfload.exe  -m 6  -b -z8  -p '%d %z %w%n'");
    puts("./anfload.exe  -m 8  -b -z15 -f /home/drmichko/web-docs/data/bst/ag-1-3-8.txt -p'%d %S%n%w%n%c%n'");
    puts("./anfload.exe  -m6 -d4 -b -p'%d %wm8 %n'");
}

void derivative ( boole f )
{ int x,u; 
  boole g = getboole();
  printf("\nD:");
  for( u = 1; u < ffsize; u++ ){
	  for ( x = 0; x < ffsize; x++ )
		  g[x] = f[x] ^ f[x^u];
	  float a = alpha( g );
	  printf(" %.2f", a );
  }
  free(g);
}

void pfboole(FILE * dst, char *format, boole f)
{
	int tempo;

    while (*format) {
	if (*format == '%') {
	    format++;
	    switch (*format) {
	    case 'x':
		panf(dst, f);
		break;
	    case 'd':
		fprintf(dst, "deg=%d", degree(f));
		break;
	    case 'a':
		fprintf(dst, "alpha=%.4f", moment( 4, 3 ));
		break;
	    case 'C':
		format++;
		correlation(f, *format - '0');
		break;
	    case 'V':
		format++;
		Vcor(f, *format - '0');
		break;
	    case 's':
		fprintf(dst, "size=%ld", fixsize);
		break;
	    case 'l':
		printf(" lin=%d", linear);
		break;
	    case 'p':
		format++;
		int cpt = 0;
		for (int i = 0; i < ffsize; i++)
		    if (cross[i] < 0)
			cpt++;
		printf(" neg=%d", cpt);
		break;
	    case 'w':
		switch ( format[1]  ) {
			case '+' :
				absdistrib( " walsh", tfr, ffsize);
				format++;
				break;
			case 'm' :
				format++;
				sscanf( format, "m%d", &tempo );
				distribmod( " walsh", tfr, ffsize, tempo);
				format++;
				while ( isdigit(*format) ) format++; 
				format--;
				break;
			default :
		   		distribution(" walsh", tfr, ffsize);
		};
		break;
	    case 'c':
		switch ( format[1]  ) {
			case '+' :
				absdistrib( " cross", &cross[1], ffsize - 1);
				format++;
				break;
			case 'm' :
				format++;
                                sscanf(format, "m%d", & tempo );
				distribmod( " cross", &cross[1], ffsize -1, tempo);
				format++;
                                while ( isdigit(*format) ) format++;
                                format--;
				break;
			default :
				distribution( " cross", &cross[1], ffsize - 1);
		};
		break;
	    case 'n':
		fprintf(dst, "\n");
		break;
	    case 'z':
		printf(" nbnz=%d", nonzero);
		break;
	    case 'D':
		derivative( f );
		break;
	    case 'M':
		format++;
		int r =  *format - '0' ;
		format++;
		int n = *format - '0' ;
		fprintf(dst, "M%d=%.4f", r, moment( r , n ));
		break;
	    default:
		fprintf(dst, "?!");
		break;
	    }
	} else
	    fprintf(dst, "%c", *format);
	format++;
    }
}

void biftok( int t[], int* pos, char *s )
{
for (char *p = strtok(s,":"); p != NULL; p = strtok(NULL, ":")){
      t[ *pos ] = atoi( p );
      printf("%s\n", p );
      *pos = *pos + 1;
}
for( int i = 0; i < *pos; i++ )
	printf("value:%d\n", t[ i ] );

}
int main(int argc, char *argv[])
{
    FILE *src;
    boole f;
    int num = 0, count = 0;
    int opt;
    extern char *optarg;
    extern int optind, opterr, optopt;

    int dim = 6;
    char *fn = NULL, *format = "%x%n%a%n";
    printf("\n#command line : ");
    for (opt = 0; opt < argc; opt++)
	printf(" %s", argv[opt]);
    int optM = 0;
    while ((opt =
	    getopt(argc, argv,
		   "a:x:r:bt:d:i:m:f:hw:p:P:l:n:s:v:z:MS:2:3R:X:%:")) !=
	   -1) {
	switch (opt) {
	case 'a':
	    optalpha = 1;
	    if (1 == sscanf(optarg, "%f:%f", &alphamin, &alphamax))
		alphamax = alphamin;
	    break;
	case 'x':
	    xvalue[optx++] = atoi(optarg);
	    break;
	case 'r':
	    biftok( rvalue, &optr, optarg );
	    break;
	case 'X':
	    Xvalue[optX++] = atoi(optarg);
	    break;
	case 'R':
	    if ( isdigit( *optarg) ) {
	    	Rvalue[optR++] =  atoi(optarg);
	    	Rvalue[optR++] = -atoi(optarg);
	    }
	    else Rvalue[optR++] = atoi(optarg);
	    break;

	case '2':
	    if ( *optarg == 'w' ) opt2w  = 1;
	    if ( *optarg == 'c' ) opt2c  = 1;
	    break;
	case '3':
	    triphase = 1;
	    break;
	case 'b':
	    optbal = 1;
	    break;
	case 'M':
	    optM = 1;
	    break;
	case 'P':
	    optl = atoi(optarg);
	    break;
	case 'm':
	    dim = atoi(optarg);
	    break;
	case 's':
	    stab = atoi(optarg);
	    break;
	case 'n':
	    optnum = atoi(optarg);
	    break;
	case 'f':
	    fn = strdup(optarg);
	    break;
	case 'p':
	    format = strdup(optarg);
	    break;
	case 't':
	    optyp = atoi(optarg);
	    break;
	case 'l':
	    oplin = atoi(optarg);
	    if (1 == sscanf(optarg, "%d:%d", &linmin, &linmax))
		linmax = linmin;
	    break;
	case 'S':
	    optsha = 1;
	    if (1 == sscanf(optarg, "%d:%d", &shamin, &shamax))
		shamax = shamin;
	    break;
	case 'd':
	    optdeg = 1;
	    if (1 == sscanf(optarg, "%d:%d", &degmin, &degmax))
		degmax = degmin;
	    break;
	case 'v':
	    optval = 1;
	    if (1 == sscanf(optarg, "%d:%d", &valmin, &valmax))
		valmax = valmin;
	    break;
	case 'w':
	    opwt = atoi(optarg);
	    break;
	case '%':
	    optmod = atoi(optarg);
	    break;
	case 'z':
	    optz = 1;
	    zmax = atoi(optarg);
	    break;
	case 'h':
	    usage(argv[0]);
	    exit(0);
	    break;
	default:		/* '?' */
	    fprintf(stderr, "%s : check usage!\n", argv[0]);
	    exit(EXIT_FAILURE);
	}
    }
    if (!fn) {
	fn = malloc(64);
	sprintf(fn, "../boole/data/class-2-%d.txt", dim);
    }
    src = fopen(fn, "r");
    if (!src) {
	perror(fn);
	exit(1);
    }
    initboole(dim);


    agsize = aglCardinality(ffdimen);
    printf("\n#AG size = %ld\n", agsize);

    while ((f = myanfloadboole(src, optM))) {
	doit(f);
	if (accept(f, optnum, num)) {
	    pfboole(stdout, format, f);
	    for (int i = 0; i < ffsize; i++)
		valtfr[abs(tfr[i])] = 1;
	    for (int i = 0; i < ffsize; i++)
		valcross[abs(cross[i])] = 1;
	    count++;
	    if (fixsize)
		total += agsize / fixsize;
	}
	num++;
    }
    fclose(src);
    printf("\n# %ld Boolean functions  in %d classes among %d\n", total,
	   count, num);

    printf("\n  cross:");
    for (int i = 0; i <= ffsize; i++)
	if (valcross[i])
	    printf(" %d", i);

    printf("\nfourier:");
    for (int i = 0; i <= ffsize; i++)
	if (valtfr[i])
	    printf(" %d", i);
    printf("\n");
    return 0;
}
