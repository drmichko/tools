
#define _XOPEN_SOURCE 500

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include "boolean.h"
#include "galois.h"
#include "nombre.h"
#include "agl.h"
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include "quad.h"
#include "mapping.h"
#include "code.h"
#include "invariant.h"



/* Boolean function are loaded sequentially
*  from a file, the degree and norme are printed.
*  dimension should be small!
*/
int xvalue[32];
int rvalue[32];

int Xvalue[32];
int Rvalue[32];
int tfr[ 64 ];
int cross[ 64 ];

int stab = 0, optval = 0;
int   optalpha = 0, oplin = 0, optx=0, optr=0, optyp=0, optdeg=0, optbal=0, opres=0, opwt = 0, optindex = 0, optnum=0, opt2=0;
int triphase =0;
int optX = 0, optR = 0;
int optmod = 0;
int degmin, degmax, optsha=0, shamin, shamax;
int opdiv=0, divmin = 0, divmax=0;
int linmin=0, linmax=0;
int optz = 0;
int zmax = 1;
nombre agsize, fixsize, total = 0;
float alphamin, alphamax;

int valcross[ 65]={0};
int valtfr[ 65 ]={0};

int optl = 2;
int valmin=0, valmax=0;

int modulo( int x, int r )
{
if ( x >=0 ) return x % r;
x = ( -x) % r;
return ( r -x ) % r;
}
void test( boole f )
{
int x, y, t, a, b;
int* F= calloc( 2*ffsize, sizeof(int) );
for( x = 0; x < ffsize; x++ )
        for( y = 0; y < ffsize; y++ ) {
                a = x ^ y;
                b = f[x] ^ f[y];
                F[ (a << 1) + b ]++;
        }
for( x = 0; x < 2*ffsize; x++ )
        for( y = x+1 ; y < 2*ffsize; y++ )
                if ( F[x] > F[y] )
                {
                        t = F[x];
                        F[x]= F[y];
                        F[y] = t;
                }
x = 0;
printf("\ndis:");
while ( x < 2*ffsize ){
        t = x;
        while ( t < 2*ffsize && F[x] == F[t] ) t++;
        printf(" %2d [%d]",  t - x , F[x] );
        x = t;
}

}

int indextwo( galois * f )
{
galois u, v;
for( u = 1; u < ffsize; u++ )
	for( v = u+1; v < ffsize; v++ )
		if ( ( f[u] ^ f[v] ) ==  f[ u ^ v ] ) return 1;
return 0;
}

void  doit( boole f )
{
galois x, y;
for( x = 0; x < ffsize; x++ )
	tfr[ x ] = 1 - 2* f[x];

for( x = 0; x < ffsize; x++)
	cross[x] = 0;
for( x = 0; x < ffsize; x++)
   for( y  = 0; y < ffsize; y++)
	cross[ x^y ] += tfr[x] * tfr[y]; 

Fourier( tfr, ffsize);
}

galois* myanfloadboole(FILE *src, int mode )
//-get the current boolean function
{
  boole res;
  char line[1024];
  while (! feof (src) )
    {
      line[0] = '\0';
      fgets (line, 1024, src);
      switch (*line)
	{
	case '#':
	  break;
	case 'a' :
	case 'h' :
	  res = strxtoboole ( line );
          if ( mode ) return res;
	  break;
        case 'f' :
	  sscanf( line, "fix=%Ld", &fixsize );
          return res;
	case 's' :
	  sscanf( line, "size=%Ld", &fixsize );
          return res;
	case 'c' :
	  sscanf( line, "card=%Ld", &fixsize );
          return res;
	default:;
	}
    }
    return NULL;
}


int poweroftwo( boole f )
{
int * fct;
galois x;
fct = booleChi(  f  );
Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++) {
	if ( ( fct[x] != 0 ) )  {
			int t = log( abs( fct[x] ) ) / log(2);
			if ( 1 << t !=  abs( fct[x]  )  ) 
			return 0;
	}
}
free( fct );
return 1;
}

int testphase( boole f )
{
galois x=0, y;
while ( tfr[x] == 0 ) x++;
for(  y= x;  y < ffsize; y++ )
	if ( tfr[y] &&  ( abs( tfr[y])  != abs(tfr[x]) )   ) return 0;
return 1;
}

int exclude( int t[] , int k, int v[] )
{
galois x;
for( int i = 0; i < k; i++ ) { 
	for( x = 0; x < ffsize; x++) 
		if ( abs( t[x] ) == v[ i ] ) 
			break;
	if (  x < ffsize ) { 
		return 0;
	}
}
return 1;
}

int allin( int t[], int k, int v[]  )
{
galois x;
for( x = 0; x < ffsize; x++) {
	int i;
	for( i = 0; i < k; i++ ) 
		if (  abs( t[x] ) ==  v[ i ] )  break;
	if ( i == k  )  
		return 0;
}
return 1;
}

void walshmode( void  )
{
int t[ 64  ]= {0};
galois x;
for( x = 0; x < ffsize; x++) 
	t[ modulo( tfr[x] , optmod )  ]++;
printf("mod:");
for( int i = 0 ; i < optmod; i++ )
	if ( t[i]) 	printf(" %d [%d]", t[i], i );
}
int balance( boole f )
{
int res = 0;
int * fct;
galois x;
fct = booleChi(  f  );
Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++)
	if ( 0 == fct[x] ) 
			res = 1;
free( fct );
return res;
}

void padic( galois *f )
{
long long int sum = 0, tmp;
int * fct;
galois x;
int r= 0, s, i;
fct = booleChi(  f  );
Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
    sum += tmp * tmp;
}
free( fct );
while ( ! (sum & 1) ) {
        r++;
        sum /=2;
}
printf( " padic:%d", r );
for( i = 0; i < 8; i++ ) {
	if ( sum & 1 ) printf(".1");
	else printf(".0");
	sum /=2;
}


printf( ": ");
}

int divisibility( boole f )
{
int res = 1024;
int * fct, v;
galois x;
fct = booleChi(  f  );
Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++)
	if ( fct[x] ) {
	     int r = 0;
	     v = abs( fct[x] );
	     while ( v  > 1 ) {
		     v /= 2;
		     r++;
	     }
	     if ( r < res ) res = r;
	}
free( fct );
return res;
}

int sha( boole f )
{
int res = 0, tmp;
galois x, y, z, t;
for( x = 0; x < ffsize; x++ )
for( y = x+1; y < ffsize; y++ )
for( z = y+1; z < ffsize; z++ ){
	t  = x ^y ^z;
	tmp = f[x]^f[y]^f[z]^f[t];
	res += (1 - 2 * tmp);
}
return res/24;
}


int zero( boole f )
{
int res = 0;
int * fct;
galois x;
fct = booleChi(  f  );
Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++)
	if (  fct[x] ) res++;

free( fct );
return res <=zmax;;
}
void lininfo( boole f )
{
int res = 0, count;
int * fct;
galois x;
fct = booleChi(  f  );
Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
	if ( abs( fct[x])  > res ){
			res = abs( fct[x] );
			count =0;
	}
	if ( abs( fct[x])  ==  res ) count++;
}
free( fct );
printf("linearity : %d [%d]", count, res );
}


void pspec( boole f )
{
int mult[ ffsize ];
int * fct = booleChi(  f  );
int x;

Fourier( fct, ffsize);

for( x = 0; x < ffsize; x++)
	mult[x]  = 0;


for( x = 0; x < ffsize; x++)
	mult[ abs(fct[x])   ]++;

int v = 0, w = 0;
int max = 0;
printf("spec : ");
for( x = 0; x <  ffsize  ; x++ )
	if ( mult[x] ) { 
	    printf(" %d [%d]", mult[ x ], x );
	    if ( x > max )  max = x;
	    if ( x % 8 == 4 ) v+= mult[x];
	    else w+=mult[x];
	}
printf(" type:%d %d walsh:%d", v, w, max);
free( fct );
}

int accept( galois *f , int optnum, int num )
{
int ok = 1;
int wt, x;
int tmp;
float t;

if (  stab )
	if ( fixsize % stab  ) 
		return  0;

if ( optnum )
	return optnum == num;

if ( oplin  ){
	int tmp =0;
	for( int x = 0; x < ffsize; x++ )
		if ( abs( tfr[x] ) >  tmp  ) tmp = abs( tfr[x] ) ;
        ok = ( linmin <= tmp ) && ( tmp <= linmax );
	
}
if ( opdiv  ){
	tmp  = divisibility ( f );
        ok = ( divmin <= tmp ) && ( tmp <= divmax );
	
}

if ( triphase  && ok ) {
   ok = testphase( f );
}

if ( optindex  && ok ) {
   ok = indextwo( f );
}

if ( opwt && ok ) {
   wt = 0;
   for( x = 0; x < ffsize; x++ )
	   wt += f[x];
   if ( wt != opwt ) ok = 0;
}

if ( ok && optbal  ) {
	int x;
	for(  x = 0; x < ffsize && tfr[x] !=  0 ; x++ );
        ok = x < ffsize;
}

if (  ok && optz ){
	ok = zero( f );
	
}
if ( ok && optalpha ){
	t = norme( f );
	ok = ( alphamin <= t ) && ( t <= alphamax );
}

if ( ok && opt2 ){
	ok  = poweroftwo( f );
	
}
if ( ok && optx ){
	ok  = exclude( tfr, optx, xvalue );
	
}
if ( ok && optr ){
	ok  = allin( tfr, optr, rvalue );
	
}
if ( ok && optX ){
	ok  = exclude( cross, optX, Xvalue );
	
}
if ( ok && optR ){
	ok  = allin( cross , optR, Rvalue);
	
}
if ( ok && optdeg  ){
	tmp  = degree( f );
	ok = ( degmin <= tmp ) && ( tmp <= degmax );
}
if ( ok && optsha  ){
	tmp  = sha( f );
	ok = ( shamin <= tmp ) && ( tmp <= shamax );
}
if ( ok && optval ){
	tmp  = valee( f );
	ok = ( valmin <= tmp ) && ( tmp <= valmax );
}
return ok;
}

void facteur( nombre z )
{
int p=2;
printf("\n#facteur=");
while ( z != 1 ) {
   if ( z % p ) p++;
   else {
	printf(".%d", p);
	while ( z % p == 0 ) z/=p;
   }
}
}

int64_t  Alfa( boole f , int t )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = 1;
     for ( int j = 1; j <= t; j++ )
	     tmp *=fct[x];
     sum += tmp;
}
sum /=  ffsize * ffsize;
free( fct );
return sum;
}

int intcmp( const void *p1, const void *p2 )
{
	return *( (int*) p1) - *( (int*) p2) ; 
}

void distribution( char* msg, int * f, int  n )
{

	printf("%s", msg );
	qsort( f, n, sizeof(int), intcmp );

	int i = 0, j;
	while ( i < n ) {
		j = i;
		while ( j < n && f[i] == f[j] ) j++;
		printf(" %d [ %d ]", j-i , f[ i ] );
		i = j;
	}
}

void bigFourier( int64_t *f, unsigned int n )
// Transformation de Fourier sur place.
{
int x,y;
int64_t z;
if (n>1) {
         bigFourier(f,n>>1);
         bigFourier(&(f[n>>1]),n>>1);
         for(x=0,y=(n>>1);  x< (n>>1); x++,y++)
                {
                z = f[x];
                f[x] = z + f[y];
                f[y] = z - f[y];
                }
         };
}


void  Vcor( boole f , int t )
{
int  * fct = calloc( ffsize, sizeof(int) );
int64_t   * big = calloc( ffsize, sizeof(int64_t) );
galois x;
for( x = 0; x < ffsize; x++ )
        big [ x ] = f[x] ? -1 : + 1;

bigFourier( big , ffsize);
for( x = 0; x < ffsize; x++){
     int64_t  tmp = 1;
     for ( int j = 1; j <= t; j++ )
	     tmp *=big[x];
     big[ x ] = tmp;
}
bigFourier( big, ffsize );

for( x = 0; x < ffsize; x++)
	fct[x] = big[ x ]  / ffsize;

for( x = 0; x < ffsize; x++){
	int tmp = fct[ x ];
	int  v = 0;
	if ( tmp < 0 ) tmp = -tmp;
	if ( tmp ) {
		while ( 0 ==  ( tmp & 1)  ){
				v++;
				tmp /= 2;
		}
	} else v = 1024;
	fct[ x ] = v;
}
fct[0] = 1024;
distribution("val:", fct, ffsize );
free( fct );
}


int   alfaX( boole f )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * tmp;
}
free( fct );
sum /= ffsize;
sum -= 3*ffsize*ffsize - 2*ffsize;
sum /= 24;
return sum;
}

void  correlation( boole f , int t )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
int64_t   * F = calloc( ffsize, sizeof(int64_t ) );
galois x;
for( x = 0; x < ffsize; x++ )
        F [ x ] = f[x] ? -1 : + 1;

bigFourier( F, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = 1;
     for ( int j = 1; j <= t; j++ )
	     tmp *=F[x];
     F[x] = tmp;
}
bigFourier( F , ffsize );
for( x = 0; x < ffsize; x++)
	fct[x] = F[x] / ffsize;

distribution( "cor:", fct, ffsize );
free( fct );
free( F );
}

int  Talfa( boole f , int p )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * fct[x];
}
free( fct );
if ( sum == 0 ) return -1;
int r = 0;
while ( (sum % p)   == 0 ) {
		sum/=p; r++;
		}
return r;
}

int Salfa( boole f, int p )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * fct[x];
}
free( fct );
if ( sum == 0 ) return -1;
int r = 0;
while ( (sum % p ) == 0 ) {sum/=p; r++;}
return sum % (p*p*p*p);
}
int  Ralfa( boole f, int p )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * tmp;
}
free( fct );
if ( sum == 0 ) return -1;
int r = 0;
while ( (sum % p ) == 0 ) {sum/=p; r++;}
return sum % (p*p*p*p);
}

int precision;
int  Falfa( boole f, int p, int* r)
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * tmp;
}
free( fct );
sum /= ffsize;
sum -=  3*ffsize*ffsize-2*ffsize;
assert ( sum % 24 == 0 );
sum /= 24;
*r = -1;
if ( sum == 0 )
	return -1;

int  v  = 0;
while ( (sum % p ) == 0 ) {
	sum  /=  p; 
	v++;
}
int i, m = 1;
for( i= 0; i < precision; i++ )
	m *= p;
sum   = sum % m;
*r   = (sum + m) % m;
return v;
}

int  Valfa( boole f, int p, int* r)
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * tmp;
}
free( fct );
if ( sum == 0 )
	return -1;
int  v  = 0;
while ( (sum % p ) == 0 ) {
	sum  /= p; 
	v++;
}

*r   = sum;

return v;
}

void  primes( boole f, int r)
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
	tmp = 1;
	for( int i = 0; i < r; i++ )
     		tmp *= fct[x];
         sum += tmp;
}
free( fct );
printf(" {" );
if ( sum == 0 )
	return; 
int p = 2;
sum = abs( sum );
while ( sum > 1 ) {
	int r = 0;
	while (  sum % p  == 0 ) { 
		r++; 
		sum/=p;
	}
	if ( r > 0 ) printf(" %d", p );
	p++;
}
printf(" }" );
}

int  Palfa( boole f, int p )
{
int64_t  sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * tmp;
}
free( fct );
return sum % p;
}

void  zeroes( boole f )
{
int  * fct = calloc( ffsize, sizeof(int) );
galois x;

for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++)
	fct[x] = fct[x]  * (1 - 2*f[x] ) ;

printf("raley:");
 pdistrib( fct , ffsize );
free(fct);

}
double   alfa( boole f )
{
long long int sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
     tmp = fct[x]* fct[x];
     sum += tmp * tmp;
}
free( fct );
return (double) sum  / ffsize / ffsize / ffsize;
}

void   vvvv( boole f, int p )
{
long long int sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++)
     fct[x] *=  fct[x]* fct[x]*fct[x];

for( x = 0; x < ffsize; x++)
        sum += fct[x];


int r = 0;
     	     while ( sum  % p == 0) {
		     sum /= p;
		     r++;
	     }

printf("\nvaluation %d-adic : %d  : ", p, r );
int i;
for( i = 0; i < 5; i++ ){
	printf("%2lld", sum % p );
	sum /= p;
}
}


void   vvv( boole f, int p )
{
long long int sum = 0, tmp;
int  * fct = calloc( ffsize, sizeof(int) );
galois x;
for( x = 0; x < ffsize; x++ )
        fct [ x ] = f[x] ? -1 : + 1;

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++)
     fct[x]= fct[x]* fct[x]*fct[x];

Fourier( fct, ffsize);
for( x = 0; x < ffsize; x++){
	fct[x] = abs( fct[x] );
     	if ( fct[x] ) {
	     int r = 0;
     	     while ( fct[x] % p == 0) {
		     fct[x] /= p;
		     r++;
	     }
	     fct[x] = r;
     }
     else fct[x] = 1024;
     
}
printf("\nval-3 at %d :", p );
pdistrib( fct, ffsize );
}


void pfboole( FILE *dst, char * format, boole f )
{
	int tmp;
	int64_t rop;
			    int v, r;
		precision=1;
		int t[64 ];
while ( *format ) {
        if ( *format == '%' ) {
                format++;
                switch( *format ) {
                 case 'x' : panf( dst, f ); break;
                 case 'd' : fprintf(dst, "deg=%d", degree(f)  ); break;
                 case 's' : pspec( f ); break;
                 case 'a' : fprintf(dst, "alpha=%.6f", alfa( f ) ); break;
		 case 'A' : format++; rop =  Alfa( f , *format - '0'); fprintf(dst, "m=%4ld %ld : %ld : %ld", rop, rop%3, rop%5, rop % 7   ); break;
		 case 'C' : format++; correlation( f, *format - '0' ); break;
		 case 'V' : format++; Vcor( f, *format - '0' ); break;
                 case 'S' : fprintf(dst, "size=%lld", fixsize); break;
                 case 'H' : fprintf(dst, "%d", alfaX( f ) ); break;
                 case 'T' : format++;
			    sscanf(format, "%d", &tmp );
			    format++;
			    fprintf(dst, "V=%d",Talfa( f, tmp ) ); 
			    break;
                 case 'F' : format++;
			    sscanf(format, "%d", &tmp );
			     do {      
                                    format++;
                            }  while ( isdigit(*format) );

			    r = 0;
			    v = Falfa( f, tmp, &r ) ; 
			    fprintf(dst, "[%d].%d", v, r); 
			    break;
                 case 'l' : format++;
			    lininfo( f  ) ; 
			    break;
                 case '2' : format++;
			    v = Valfa( f, 2, &r ) ; 
			    fprintf(dst, "[%d].%d", v , r); 
			    break;
                 case 'p' : format++;
			    int cpt = 0;
			    for( int i = 0; i < ffsize; i++ )
				    if ( cross[i] < 0 ) cpt++;
			    printf(" neg=%d", cpt );
			    break;
                 case 'R' : format++;
			    sscanf(format, "%d", &tmp );
			    format++;
			    fprintf(dst, "R=%d", Ralfa( f, tmp ) ); 
			    break;
                 case 'P' : format++;
			    sscanf(format, "%d", &precision );
			    do {
				    format++;
			    }  while ( isdigit(*format) );
			    break;
                 case 'w' : if ( optmod ) walshmode( ); else pwalshf( stdout, f ) ; break;
		 case 'c' : printf("cross: ");
			    if ( optmod) 
				    for( int i = 0; i < ffsize; i++ ){
						t[i] = modulo( cross[ i ], optmod ) ;
			    } else for( int i = 0; i < ffsize; i++ )
						t[i] = cross[ i ];
			    pdistrib( &t[1] , ffsize - 1   ) ; break;
                 case 'W' : pWalshf( stdout, f ) ; break;
                 case 'v' : format++;
		    if ( ! sscanf(format, "%d", &tmp ) ) tmp=2;
			    vvv(f, tmp); 
			    break;
                 case 'n' : fprintf(dst, "\n"); break;
                 case 'z' : zeroes( f ); break;
                 default  : fprintf(dst, "???"); break;
                }
          } else fprintf(dst, "%c", *format);
         format++;
}
}

int main( int argc, char*argv[] )

{
  FILE *src;
  boole f;
  int num = 0, count = 0;
  int opt;
   extern char *optarg;
       extern int optind, opterr, optopt;

  int dim=6;
  char *fn = NULL , *format="%x%n%a%n";
  printf("\n#command line : ");
  for( opt = 0; opt < argc; opt++ )
	  printf(" %s", argv[opt] );
 int optM = 0;
  while ((opt = getopt(argc, argv, "-A:a:x:r:bt:d:D:i:m:f:w:p:P:l:n:s:v:z:MS:23R:X:%:")) != -1) {
               switch (opt) {
               case 'a':
                   optalpha = 1; 
		   if ( 1 == sscanf(optarg,"%f:%f", &alphamin, &alphamax ) )
			   alphamax = alphamin;
                   break;
               case 'A':
                   optalpha = 1; 
		   alphamin = atof( optarg );
		   alphamax = alphamin+0.001;
		   alphamin = alphamin-0.001;
                   break;
               case 'x':
                   xvalue[ optx++] = atoi( optarg );
                   break;
               case 'r':
                   rvalue[ optr++] = atoi( optarg );
                   break;
               case 'X':
                   Xvalue[ optX++] = atoi( optarg );
                   break;
               case 'R':
                   Rvalue[ optR++] = atoi( optarg );
                   break;

               case '2':
                   opt2=1;
                   break;
               case '3':
                   triphase=1;
                   break;
               case 'b':
                   optbal=1;
                   break;
               case 'M':
                   optM=1;
                   break;
               case 'P':
                   optl = atoi( optarg );
                   break;
               case 'm':
                   dim = atoi( optarg );
                   break;
               case 's':
                   stab = atoi( optarg );
                   break;
		case 'n':
                   optnum = atoi( optarg );
                   break;
               case 'i':
                   optindex = atoi( optarg );
                   break;
               case 'f':
                   fn = strdup( optarg);
                   break;
               case 'p':
                   format   = strdup( optarg  );
                   break;
               case 't':
                   optyp   = atoi( optarg );
                   break;
	        case 'D':
                   opdiv   = atoi( optarg );
		   if ( 1 == sscanf(optarg,"%d:%d", &divmin, &divmax) )
			   divmax = divmin;
                   break;
	        case 'l':
                   oplin   = atoi( optarg );
		   if ( 1 == sscanf(optarg,"%d:%d", &linmin, &linmax) )
			   linmax = linmin;
                   break;
               case 'S':
                   optsha = 1;
		   if ( 1 == sscanf(optarg,"%d:%d", &shamin, &shamax) )
			   shamax = shamin;
                   break;
               case 'd':
                   optdeg = 1;
		   if ( 1 == sscanf(optarg,"%d:%d", &degmin, &degmax) )
			   degmax = degmin;
                   break;
               case 'v':
                   optval = 1;
		   if ( 1 == sscanf(optarg,"%d:%d", &valmin, &valmax) )
			   valmax = valmin;
                   break;
               case 'w': opwt = atoi( optarg );
                   break;
               case '%': optmod = atoi( optarg );
                   break;
               case 'z': optz = 1;
			 zmax = atoi( optarg);
                   break;
               default: /* '?' */
                   fprintf(stderr, "%s : check usage!\n", argv[0]);
                   exit(EXIT_FAILURE);
               }
           }
  if ( ! fn ){
	  fn = malloc( 64 );
	  sprintf(fn, "../../data/class-2-%d.txt", dim );
  } 
  src = fopen( fn , "r");
  if ( ! src ) {
	perror(fn);
	exit(1);
  }
  initfield( dim ,  0);


  agsize = aglgrpsize( ffdim );
  printf("#AG size = %Ld\n", agsize );

  while ( (f = myanfloadboole( src , optM ) ) ){
	doit( f ); 
	if ( accept( f, optnum, num  ) )  {
                pfboole(stdout, format, f );
		for( int i = 0; i < ffsize; i++ )
			valtfr[ abs( tfr[i] ) ] = 1;
		for( int i = 0; i < ffsize; i++ )
			valcross[ abs( cross[i] ) ] = 1;
		count++;
		if ( fixsize )
			total += agsize / fixsize ;
        }
	num++;
  }
  fclose( src );
  printf("\n# %Ld Boolean functions  in %d classes among %d\n", total, count, num);

	printf("\n  cross:");
  	for( int i = 0; i <=  ffsize; i++ )
	  if ( valcross[i] ) printf(" %d", i );

	printf("\nfourier:");
  	for( int i = 0; i <=  ffsize; i++ )
	  if ( valtfr[i] ) printf(" %d", i );
	printf("\n");
return 0;
}
