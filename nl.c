#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"
#include "math.h"
#include "code.h"
typedef  struct _list_ {
	 boole  fct;
	 struct _list_ * next;
} enrlist, *liste;


void append( boole f, liste *l )
{  liste aux = malloc( sizeof( *aux ) );
   aux->fct = getboolecpy( f );
   aux->next = *l;
   *l = aux;
}

int length( liste l )
{ int r = 0;
	while ( l ) {
		r++;
		l=l->next;
	}
	return r;
}

int freelength( liste l )
{ int r = 0;
	while ( l ) {
		liste tmp = l;
		r++;
		l=l->next;
		free( tmp->fct );
		free( tmp );
	}
	return r;
}
liste  listing( boole f , int k, int m, int goal )
{       liste  l = NULL;
        code    c = rmcode( 0, k, ffdimen  );
        boole   t = getboole( );
        int wt = weightBoole( f );
        if ( wt <= goal ) 
		append( t, &l );
        int  cpt=1, limite = 1 << c.nbl;
        for( int x = 0; x < ffsize; x++ )
                        t[x] = f[x];
        while (  cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < ffsize; x++ )
                        t[x] ^= c.fct[i][x];

                wt = weightBoole( t );
                if ( wt < goal )
			append( t , & l);
                cpt++;
        }
        freecode( c );
	return l;
}

liste doit( boole f, int k, int m, int r )
{
if ( m == 6 ) {
	liste l = listing( f, k, m, r );
	printf("\nlen=%d\n", freelength( l ) );
	return l;
}
return NULL;
}

int main(int argc, char *argv[])
{
    FILE *src = NULL;
    char *anf = NULL;
    int opt, r = 0;
    int iter=1000;
    int k = 2;
    while ((opt = getopt(argc, argv, "a:k:m:f:r:hi:vs")) != -1) {
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
	    case 'k' : k = atoi( optarg) ; break;
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

    if (src) {
	while ((f = loadBoole(src))) {
	    panf(stdout, f);
    	    printf("\niter=%d deg=%d lin=%d :", iter, degree(f) , linearity(f) );

	    doit( f, k,  ffdimen, r );
	}
	fclose(src);
    }


    return 0;
}
