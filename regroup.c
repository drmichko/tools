#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "boolean.h"
#include "math.h"
#include "mapping.h"



int main(int argc, char *argv[])
{
    FILE *src = NULL;
    int opt;
    while ((opt = getopt(argc, argv, "a:m:f:r:hi:vs")) != -1) {
	switch (opt) {
	case 'm':
	    initboole(atoi(optarg));
	    break;
	}
    }



    src = fopen( argv[optind++], "r");
    if (!src) {
		perror(optarg);
		return 1;
	    }

    initagldim(ffdimen);

  boole  f;
  int num = 0, value;
  boole table[ 16000 ];
  while (   ( f =  loadBooleValue( src , &value ) ) )  {
	    table[num++] = f;
	    }
  fclose(src);


  printf( "\nnum=%d\n", num );

    src = fopen( argv[optind], "r");
    if (!src) {
		perror(optarg);
		return 1;
	    }

  uint64_t size;
  aglGroup g = NULL;
  int item = 0, class=0;
  while (   (  f =  loadaglboolesize( src , &g , &size ) )  )  {
	    int i = 0;
	    for ( i = 0; i < num && boolecmp( f, table[i] ) != 1; i++ );;
	    if ( i != num ) {
		panf( stdout, f );
		    paglGroup( stdout, g );
		    printf("\nstabSize=%ld", size );
		    item++;
	    }
	    free( f );
	    aglfreeGroup(  g );
	    g = NULL;
	    class++;
	    }
  fclose(src);
  printf("\n#class = %d, item : %d / %d\n", class, item, num );

  return 0;
}
