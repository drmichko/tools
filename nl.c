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
	 int wt;
	 struct _list_ * next;
} enrlist, *liste;
int limite;
code cc[9][9];



int K[5][9] = {
    {0, 0, 0, 4, 8, 16, 32, 64, 128},
    {0, 0, 0, 2, 6, 12, 28, 56, 120},
    {0, 0, 0, 1, 2, 6, 18, 40,  96},
    {0, 0, 0, 0, 1, 2, 8, 20,   60},
    {0, 0, 0, 0, 0, 1, 2, 8,    26}
};

void append( int wt, int m, boole f, liste *l )
{  liste aux = malloc( sizeof( *aux )  );
   aux->fct = calloc( 1 << m, 1 );
   for( int x = 0; x < (1<<m) ; x++ )
	   aux->fct[x] = f[ x ];
   aux->wt = wt;
   aux->next = *l;
   *l = aux;
}

void glue( int wt, int m, boole L, boole R,  liste *l )
{  liste aux = malloc( sizeof( *aux )  );
   int q = 1 << m;
   aux->fct = calloc( 2 * q, 1 );
   for( int x = 0; x < q ; x++ ){
	   aux->fct[x] = L[ x ];
	   aux->fct[x+q] = R[ x ];
   }
   aux->wt = wt;
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

int freeliste( liste l )
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

int check( liste l , int R )
{
       if ( l == NULL   ) return 0;	
	while ( l ) {
		if ( l->wt < R ) return 0;
		l = l->next;
	}
	return 1;
}
int showliste( char *msg, liste l  )
{ int r = 0;
	printf("\n%s:", msg);
	while (  l ) {
		r++;
		printf(" %d", l->wt );
		l = l->next;
	}
	printf("\ncard:%d\n", r);
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
		wt += f[x]; 
	}
	return wt;
}

int    count = 0;

liste  listing( boole f , int k, int m, int R )
{       liste   l = NULL;
        code    c = cc[ k ][ m ];
        boole   t = getboole( );
	for( int x = 0; x < ( 1 << m ); x++ )
		t[x] =f[x];
        int wt =  wtboole( t, m );
        if ( wt <=  R ) 
		append( wt, m, t, &l );
        int  cpt = 1, limite = 1 << c.nbl;
        while (  cpt < limite ) {
                int i = __builtin_ctz( cpt  );
        	int wt =  addwtboole( t, c.fct [i], m );
                if ( wt <=  R ) {
			append( wt, m, t , & l);
		}
                cpt++;
        }
	free( t );
	return l;
}

liste doit( boole f, int k, int m, int  rho )
{
if ( m == limite  ||  k == 1 ) {
	liste l = listing( f, k, m,  rho );
	return l;
}

liste  L = doit( f, k, m - 1, rho / 2 );
uchar sum[ ffsize ];
int q   = 1 << ( m-1) ;
liste res = NULL;


liste  l = L;
while ( l ) {
	for( int x = 0; x < q; x++  )
	   sum[x] = l->fct [x] ^ f[ x + q ] ^ f[x];
	int wt = wtboole( l->fct, m - 1 );
	liste  lr = doit( sum , k - 1, m - 1, rho - wt  );
		liste tmp = lr;
        	while ( tmp ) {
			glue( tmp->wt + wt, m - 1, l->fct, tmp->fct , & res  );
			tmp = tmp->next;
		}
		freeliste( lr );
	l = l->next;
}
freeliste( l );

liste R = doit( & f[q], k, m - 1, rho / 2 );
l = R;
while ( l ) {
        for( int x = 0; x < q; x++  )
           sum[x] = l->fct [ x ] ^ f[ x + q ] ^ f[x];
        int wt = wtboole( l->fct, m - 1 );
        liste  lr = doit( sum , k - 1, m - 1, rho - wt  );
        liste tmp = lr;
        while ( tmp ) {
                glue( tmp->wt + wt, m - 1, tmp->fct, l->fct  , & res  );
                tmp = tmp->next;
        }
        freeliste( lr );
        l = l->next;
}
freeliste( l );

return res;
}



int main(int argc, char *argv[])
{
    FILE *src = NULL;
    char *anf = NULL;
    int opt,  R = 0;
    int k = 2;
    limite = 5;
    int job = 0;
    while ((opt = getopt(argc, argv, "a:k:m:f:R:hij::vsl:")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	case 'a':
	    anf = strdup(optarg);
	    break;
	case 'j':
	    job = atoi(optarg);
	    break;
	case 'l':
	    limite = atoi(optarg);
	    break;
	case 'R':
	    R = atoi(optarg);
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

    for( int k = 0; k < 5; k++ ){
	    for( int m = 0; m < 9; m++ )
		    printf(" %3d", K[k][m] );
	    printf("\n");
    }
    for( int k = 0; k <=  4 ; k++ )
    for( int m = limite; m < 9 ; m++ )
		    cc[ k ][ m ]  = rmcode( 0, k,  m  );
    boole f;

    int tour = 0;
    int soluce = 0;
    if (src) {
	while ((f = loadBoole(src))){ 
            if ( valuation(f) >= 0  ) 	{
	    if ( tour == job ) {
		    liste l = doit( f, k,  ffdimen,  R );
	   		 if ( check( l, R )  ) {
		    soluce++;
		    printf("\n");
	    	    panf(stdout, f);
		    //showliste( "wt", l );
	    }
	    freeliste( l );
	    }
	    }
	 free( f );
	tour++;
	fflush(stdout );
	}
	fclose(src);
    }

    printf("\n");
    printf("\n");
    return 0;
}
