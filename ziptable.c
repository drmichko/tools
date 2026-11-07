
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "boolean.h"
#include "orbitools.h"
#include "vectorial.h"


void saveZipname(char* dirname, char *name, basis_t * b  )
{
    FILE *dstFile;
    char fn[128];
    size_t size = ((size_t) 1) << b->mat.nbl;
    sprintf(fn, "%s/%s-%d-%d-%d.zip", dirname, name, b->s, b->t, b->m);
    dstFile = fopen(fn, "w");
    if (!dstFile) {
        perror(fn);
        exit(1);
    }
    uint64_t count = 0;
    uchar *z = calloc(size / 8, 1);
    for (vector v = 0; v < size; v += 8) {
        uchar val = 0;
        for (unsigned int i = 0; i < 8; i++)
            if ( b->table[v + i]  ) {
                val |= ((uchar) 1) << i;
		count++;
	    }
        z[v / 8] = val;
    }
    size = size / 8;
    if (fwrite(z, 1, size, dstFile) != size)
        printf("problem\n");
    fclose(dstFile);
    free(z);
    printf("\n#count:%ld", count );
}



int option(int argc, char *argv[]);

int optt, opts, dimen;
char *name  = NULL;
char *oname = "foo";

int main(int argc, char *argv[])
{

    assert( dimen > 0);
	initboole( dimen );
    assert( opts  > 0);
    assert( optt  >= opts);
    assert( name != NULL );
    basis_t base = monomialBasis( opts, optt, ffdimen);
    char cname[64];
    loadBasename( name , &base);
    saveZipname( "data", cname, &base );

    return 0;
}

void help( void )
{
system( "sed -n '/^[a-zA-Z_].*\\boption\\s*(/,/^}/p' stabredo.c");
}

int option(int argc, char *argv[])
{   int opt;
    while ((opt = getopt(argc, argv, "m:s:vt:rs:f:g:w:hW")) != -1) {
	switch (opt) {
	case 'f':
	    name = strdup(optarg);
	    break;
	case 'w':
	    oname = strdup(optarg);
	    break;
	case 'm':
	    dimen = atoi(optarg);
	    break;
	case 't':
	    optt = atoi(optarg);
	    break;
	case 's':
	    opts = atoi(optarg);
	    break;
	case 'h' : help();
		   return 0;
	default:		/* '?' */
	    fprintf(stderr, "check usage of %s\n", argv[0]);
	    exit(EXIT_FAILURE);
	}
    }
    return 1;
}

