
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
    char fn[256];

    const char *home = getenv("HOME");
    sprintf(fn, "%s/gitub/nonlinearity/tuyau/NL-2-7-36.dat" , home );

    FILE *src = fopen( fn , "r");
    if (!src) {
	perror( fn );
	exit(1);
    }
    aglGroup g = NULL;
    int num = 0, cls = -1;
    int opt, optw = 0, optr=0;
    float beta = 0;
    int deg = 0;
    while ((opt = getopt(argc, argv, "ic:wb:d:m:r")) != -1) {
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

	case 'r':
	    optr++;
	    break;
	case 'b':
	    beta = atof(optarg);
	    break;
	default:
	    exit(0);
	}

    }


    initboole(7);
    initagldim(7);

	grp = mkaglGroup( ffdimen );
	size = aglcard(ffdimen);

    boole f;
    num = 0;
    uint64_t size;
    int val = 0;
    while ((f = loadBooleValue(src, &val))) {
            panf( stdout, f );
            doit( f );
	    num++;
	    free(f);
	}
    fclose(src);
    printf("\n#maps : %d\n", num );
    return 0;
}
