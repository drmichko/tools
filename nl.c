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

int wtboole ( boole f, int m )
{
	int wt = 0;
	for( int x = 0; x < ( 1 << m ); x++ )
                wt += f[x];
	return wt;
}

int  addwtboole ( boole f, boole g, int m )
{
	int wt = 0;
	for( int x = 0; x < ( 1 << m ); x++ ) {
                f[x] ^=  g[x];
		wt = f[x]; 
	}
	return wt;
}

liste  listing( boole f , int k, int m, int goal )
{       liste  l = NULL;
        code    c = rmcode( 0, k, m );
        boole   t  = getboole( );
	for( int x = 0; x < ( 1 << m ); x++ )
		t[x] =f[x];
        int wt =  wtboole( t, m );
        if ( wt <= goal ) 
		append( t, &l );
        int  cpt=1, limite = 1 << c.nbl;
        for( int x = 0; x < ( 1 << m ) ; x++ )
                        t[x] = f[x];
        while (  cpt < limite ) {
                int i = __builtin_ctz( cpt  );
        	int wt =  addwtboole( t, c.fct [i], m );
                if ( wt <=  goal )
			append( t , & l);
                cpt++;
        }
        freecode( c );
	return l;
}

liste doit( boole f, int k, int m, int r )
{
if ( m == 5 ) {
	liste l = listing( f, k, m, r );
	return l;
}
liste ll = listing( f, k, m-1, r/2 );
printf("\nL=%d\n", length( ll ) );
uchar s[ 256 ];
int q   = 1 << ( m-1) ;
int bad = 0;
while ( ll ) {
	for( int x = 0; x < q; x++  )
			s[x] = ll->fct [x] ^ f[ x + q ] ^ f[x];
	int wt = wtboole( ll->fct, m - 1 );
	liste  tmp = listing( s , k-1, m-1, r - wt  );
/*	
        while ( tmp ) {
		if (   wt + wtboole( tmp->fct, m-1 )  <  r ) 
			printf(" %d %d\n", wt,  wtboole( tmp->fct, m-1 ) ) ;
		tmp = tmp->next;
	}
	*/
	int nb = freelength( tmp );
	if ( nb  > 0  )  bad++;
	ll = ll->next;
}

printf("\nbad=%d", bad  );

liste lr = listing( & (f[ 1 << (m-1) ]) , k, m-1, r/2 );
printf("\nR=%d\n", freelength( lr ) );

return NULL;
}

int main(int argc, char *argv[])
{
    FILE *src = NULL;
    char *anf = NULL;
    int opt,  R = 0;
    int iter=1000;
    int k = 2;
    while ((opt = getopt(argc, argv, "a:k:m:f:R:hi:vs")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'a':
	    anf = strdup(optarg);
	    break;
	case 'R':
	    R = atoi(optarg);
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
	    doit( f, k,  ffdimen,  R );
	}
	fclose(src);
    }


    return 0;
}
