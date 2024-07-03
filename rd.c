#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include "boolean.h"

void RD( boole f, int iter )
{
  int t  = ffdimen/2;

  shortvec trans, base[ ffdimen/2 ];

  int i, r = 0;
  shortvec   p;
  int limite  = 1 << t;

  int q = ffsize;
  initboole( t );
  boole g = getboole();
  int best = t;
  shortvec x, y;
  int cpt = 0;
  while ( iter--) {
	  r = 0;
	  cpt++;
	  while ( r < t ) {
		 y = random() %  q ;
           	 for( i = 0; i < r; i++) {
                    p = base[i] & ( base[i] - 1 );
                    p ^= base[i] ;
                    if ( p & y ) y ^= base[i];
                 }
                 if ( y != 0 ) base[r++] = y;
	  }
	trans = random() % q;
  	for( x = 0; x < limite; x++ ){
         	y = 0;
          	for( i = 0; i < r; i++ )
                  if ( (x & ( 1<<i ) ) )   
			  y^= base[ i ];
                g[x] =  f[y ^ trans ];
                }
	int  d = degree( g ) ;
	if ( d < best ) {
		best = d;
		printf("\niter=%d score=%d\n", cpt, best );
	}

  }
 initboole( 2 * t );
}

int main(int argc, char *argv[])
{


    initboole(atoi(argv[1]));
    srandom(getpid() + time(NULL));

    boole f;
    FILE *src = fopen(argv[2], "r");

    if (!src)
	return 1;

    while ((f = loadBoole(src))) {
	panf(stdout, f);
	printf("\ndeg=%d linearity=%d\n", degree(f), linearity(f) );
	RD( f, 100000 );

    }
    return 0;
}
