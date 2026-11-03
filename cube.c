#include <stdio.h>
#include <stdlib.h>
#include "galois.h"
#include "suites.h"
#include "boolean.h"
#include "basis.h"

	int
main (int argc, char *argv[])
{

  int m, s;
  galois a, b , x;
  m = atoi (argv[1]);
  printf ("\n---\nm=%d ", m);
  initfield (m, 0);
  printf("\n");
  boole f = getboole( );;
  

  for( a = 0; a < ffsize; a++)
      f[  a   ] = 1; 
  
  printf("\ndegree=%d\n", degree( f ) ) ;

  return 0;
}

