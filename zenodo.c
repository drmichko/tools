
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "distrib.h"
#include "orbitools.h"

boole strdigitoboole( char *bfr )
{ char   *ptr;
  galois x = 0;
  int v = 0;
  boole res;
  res = getboole();
  ptr = bfr;
  while ( *ptr ) {
    if ( isdigit( *ptr) ) x ^= 1 << (*ptr - '0');
    if ( *ptr == '+' ) {
      if ( x == 0 ) res[ x ] = v; else res[ x ] = 1;
      x = 0;
    }
    ptr++;
  }
  if ( x == 0 ) res[ x ] = v; else res[ x ] = 1;
  xform( res, ffsize );
  return res;
}

int qlinearity( boole f, int qsize )
{ int tfr[ qsize ];
  int x;
  int res = 0;
  for( x = 0; x < qsize; x++ )
          tfr[x] = f[x] ? -1 : +1;
  Fourier( tfr, qsize );
  for( x = 0; x < qsize; x++ )
     if ( abs( tfr[x] ) >  res )
             res =  abs( tfr[x] );
  return res;
}


code   rmq;


int  myNQ( boole f, int  qsize )
{ 
        boole  t = getboole( );
        int demi = qsize / 2;
        int wt = 0;
        for( int x = 0; x < qsize; x++ )
		wt+= f[x];
        int  cpt=1, limite = 1 << rmq.nbl;
                for( int x = 0; x < qsize; x++ )
                        t[x] = f[x];
        while ( wt == demi && cpt < limite ) {
                int i = __builtin_ctz( cpt  );
                for( int x = 0; x < qsize; x++ )
                        t[x] ^= rmq.fct[i][x];

                wt = 0;
        	for( int x = 0; x < qsize; x++ )
			wt+= f[x];
                cpt++;
        }

        free( t );
        if ( wt < 0 ) wt = -wt;
        return wt;
}


int* devQuadDerivation( boole f , int q  )
{
shortvec u, p, msk, x, y;
int tmp;
int *t = calloc( q, sizeof(int ) );
vector v;
boole der;

der = calloc( q / 2  , sizeof(uchar) );

for( u = 1; u < q; u++ ){
        p = 1;
        while ((p & u) == 0)
                p <<= 1;
        msk =  p ^ ( q / 2) ;
        for( x = 0 ; x < q/2 ; x++ ){
                y = x;
                if ( y & p ) y^=msk;
                der[ x ] = f[ y ] ^ f[y ^ u ];
        }
        //tmp = qlinearity( der, q / 2 );
        tmp = myNQ( f, q/2 );
        t[u] = tmp;
}
free(der);
return t;
}
int invkersize( boole f, int r, int s, int t, int m  )
{
int  k;
int tmp;
vector v;
code cc = rmcode( 1 , r,  m );
code prd = multicode( cc, f  );
for( k=0; k < prd.nbl; k++ )
        rsreduction( prd.fct[k], r+s  , r+t, 1 << m ); 
code ker = kernel( cc, prd );
freecode( cc );
freecode( prd );
int res = ker.nbl;
freecode( ker );
return res;
}


int invQuadDerivation( boole f , int q,  void **r1, int *c1, int mode)
{
int tmp;
int* t;
t = devQuadDerivation( f, q  );
tmp = findtable( t, q, r1, c1 , mode);
return tmp;
}


int invKerDerivation( boole f , int q,  void **r1, int *c1, int mode)
{
int tmp;
int* t;
t = devQuadDerivation( f, q  );
tmp = findtable( t, q, r1, c1 , mode);
return tmp;
}


int countbis = 0;
void *rootbis = NULL;

int* devliftRestriction( boole f )
{
shortvec u, p, msk, x, y, w;
int t1, t2;
int *tp = calloc( ffsize, sizeof(int ) );
vector v;
boole F, G;
F = calloc( ffsize, sizeof(uchar) );
G = calloc( ffsize, sizeof(uchar) );
tp[0] = 0;
int tq[2];
int q = ffsize;
for( u = 1; u < q; u++ ){
        p = 1;
        while ((p & u) == 0)
                p <<= 1;
        msk =  p  ^  ( q/2 );
        for( x = 0 ; x < q/2; x++ ){
                y = x;
                if ( y & p ) y^=msk;
                w = weight( y & u ) & 1;
                if ( w ) y^= p;
                F[x] = f[y];
                G[x] = f[y^p];
        }
        t1 = invkersize( F, 3, 3, 3, 8 );
        t1 = invkersize( G, 3, 3, 3, 8 );

        if ( t1 > t2 ) {
                  int t = t1;
                  t1 = t2;
                  t2 = t;
         }
          //tp[ u ] =  t1<t2 ? t1 * 100 + t2 : t2 * 100 + t1;
          tq[0] = t1;
          tq[1] = t2;
          tp[ u ] =  findspltable( tq, 2, &rootbis, &countbis );
}
free(F);
free(G);
return tp;
}


int invliftRestriction( boole f, void **root, int * count , int mode)
{
int tmp;
int *t = calloc( ffsize, sizeof(int ) );
t = devliftRestriction( f  );
tmp  = findtable( t, ffsize , root, count, mode);
return tmp;
}



int main(int argc, char *argv[])
{  

    char *fn = NULL;


    aglGroup g = NULL;
    int num = 0, cls = -1;
    int opt, optw = 0, optr=0, optinit = 0;
    int deg = 0, dimen = 7;
    int target = 3;
    int verb   = 0;
    int mode = 0;
    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:vx:")) != -1) {
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
	case 'i': 
		optinit = atoi( optarg);
	break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'x': 
		mode   = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    optr=atoi( optarg );
	    break;
	case 't':
	    target=atoi( optarg );
	    break;
	case 'v':
	    verb++;
	    break;
	default:
	    exit(0);
	}

    }


    initboole(  dimen );
    initagldim( dimen );

 	rmcode( 1 , 2, ffdimen - 1 );
    
    num = 0;
    uint64_t stabSize, orbSize;
   void  *rootj = NULL;int   countj = 0;
   void  *rootQ = NULL;int   countQ = 0;
   void  *rootR = NULL;int   countR = 0;

 
    char line[ 2048 ];
    FILE * src = fopen( argv[optind], "r" );

    if ( ! src ) { perror( argv[optind] ); exit(1);}
  
    while ( fgets( line, 2048, src ) ) {
         int no;
         char tmp[1024];
         sscanf( line, "%d %s %ld", &no, &tmp, &stabSize);
	 //printf("no=%d %s %ld\n", no, tmp, stabSize );
	 boole f = strdigitoboole( tmp  );
         if ( verb  ) {
			panf( stdout, f );
	 }
         free( f );
	int R[8];
        int n = 0;
		R[n++] = stabSize;
		R[ n++ ] = invkersize( f, 3 , 3, 3, 10 );
		R[ n++ ] = invliftRestriction( f, &rootR, &countR, mode );
                int j = findspltable(R, n, &rootj, &countj);

	 num++;
        if ( nouvelle)  printf("\n#maps : %d cont=%d\n", num , countj);
     } 

     fclose(src);

     printf("\n#maps : %d cont=%d\n", num , countj);
     return 0;
}
