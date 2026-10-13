// ============================================================
// Cover set de B(3,3,7) sous Stab(h)
// Convention : T[i] = 0 (non visité), 1 (cover set), 2 (couvert)
// ============================================================


#define _GNU_SOURCE
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#define _GNU_SOURCE
#include <unistd.h>
#include <sys/syscall.h>

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "orbitData.h"

#include "boolean.h"
#include "orbitools.h"

// ============================================================
// Configuration
// ============================================================

#define N_ELEMENTS   (1ULL << 35)   // 2^35 = 34 359 738 368
#define NB_THREADS   72
#define FILE_TAILLE  (1 << 20)      // 1M éléments par file locale

// Valeurs de T
#define NON_VISITE   0
#define COVER_SET    2
#define COUVERT      1

// ============================================================
// Structures globales
// ============================================================

static atomic_uint_fast64_t curseur;        // Curseur global atomique


basis_t base;
aglVectorGroup lda;


// ============================================================
// Fonctions utilitaires
// ============================================================

// CAS sur un byte
static inline bool cas_byte(uint8_t *ptr, uint8_t expected, uint8_t new_val) {
    return atomic_compare_exchange_strong(
        (_Atomic(uint8_t)*)ptr, &expected, new_val
    );
}

// Cherche un indice i tel que T[i] == 0
// Retourne UINT64_MAX si aucun trouvé
static inline uint64_t trouver_non_visite(void) {
    while (true) {
        uint64_t i = atomic_fetch_add(&curseur, 1);
        if (i >= base.size ) {
            return UINT64_MAX;
        }
        if (  base.table[i] == NON_VISITE) {
            return i;
        }
    }
}

// ============================================================
// Exploration d'une orbite
// ============================================================

// Explore l'orbite de x (déjà marqué COVER_SET)
// Marque tous les éléments atteints avec COUVERT
static void explorer_orbite(uint64_t x) {
    // File locale (BFS)
    uint64_t *file = malloc(FILE_TAILLE * sizeof(uint64_t));
    if (!file) {
        fprintf(stderr, "Erreur malloc file\n");
        return;
    }
    
    uint64_t debut = 0, fin = 0;
    file[fin++] = x;
   
    aglVectorGroup aux;
    uint64_t count = 0;
    while (debut < fin) {
        uint64_t u = file[debut++];
        aux = lda;
        while (aux != NULL) {
            uint64_t v  = aglVectorImage( u , aux->act);
            // Tenter de marquer v comme couvert
            if (  base.table[ v ] ==  NON_VISITE  ) {
                // v est nouveau : l'ajouter à la file
		base.table[ v ] = COUVERT;
		count++;
                if (fin < FILE_TAILLE) {
                    file[fin++] = v;
                } else {
                    // abandon
                }
            }
            // Sinon : déjà marqué (COVER_SET ou COUVERT), ne rien faire
            aux = aux->next;
        }
        
        // Réinitialiser la file si elle est vide
        if (debut == fin) {
            debut = 0;
            fin = 0;
        }
    }
    
	if ( count )
		printf("\n %2d : %ld\n", gettid() , count );
    free(file);
}

// ============================================================
// Thread worker
// ============================================================

static void *thread_worker(void *arg) {
    (void)arg;
    
    while (true) {
        // 1. Chercher un élément non visité
        uint64_t x = trouver_non_visite();
        if (x == UINT64_MAX) {
            break;  // Plus rien à explorer
        }
        
        // 2. Tenter de le marquer comme élément du cover set
        if (! cas_byte( & base.table[ x ], NON_VISITE, COVER_SET)) {
            continue;  // Un autre thread l'a pris
        }
        
        // 3. Explorer l'orbite de x
        explorer_orbite(x);
    }
    
    return NULL;
}

// ============================================================
// Extraction du cover set
// ============================================================

static void extraire_cover_set(void) {
    uint64_t nb_cover = 0;
    uint64_t nb_couvert = 0;
    uint64_t nb_non_visite = 0;
    
    for (uint64_t i = 0; i < base.size ; i++) {
        switch ( base.table[ i ]) {
            case COVER_SET:   nb_cover++;      break;
            case COUVERT:     nb_couvert++;    break;
            case NON_VISITE:  nb_non_visite++; break;
        }
    }
    
    printf("=== Résultats ===\n");
    printf("Cover set      : %llu éléments\n", (unsigned long long)nb_cover);
    printf("Couverts       : %llu éléments\n", (unsigned long long)nb_couvert);
    printf("Non visités    : %llu éléments\n", (unsigned long long)nb_non_visite);
    printf("Total          : %llu éléments\n", (unsigned long long)N_ELEMENTS);
}

// ============================================================
// Initialisation des générateurs (exemple factice)
// ============================================================




size_t coverset(  int nb ) {
    // Initialiser le curseur
    atomic_init(&curseur, 0);
    
    // 4. Lancer les threads
    printf("\nLancement de %d threads...\n", nb );
    pthread_t threads[NB_THREADS];
    
    for (int i = 0; i < nb ; i++) {
        if (pthread_create(&threads[i], NULL, thread_worker, NULL) != 0) {
            fprintf(stderr, "Erreur pthread_create\n");
            return 1;
        }
    }
    
    // 5. Attendre la fin
    for (int i = 0; i < nb ; i++) {
        pthread_join(threads[i], NULL);
    }
    
    // 6. Extraire le cover set
    extraire_cover_set();
    
    return 0;
}

int main(int argc, char *argv[])
{  

    char *fn = "../boole/data/B-4-5-7.dat";


    aglGroup g = NULL;
    int opt; 
    int dimen = 7;
    int target = 0;
    int nbp = 2;

    while ((opt = getopt(argc, argv, "i:f:c:wb:d:m:r:t:p:")) != -1) {
	switch (opt) {
	case 'w':
	    break;
	case 'd':
	   break;
	case 'c':
	    break;
	case 'p': 
	    nbp = atoi( optarg );
	break;
	case 'm': 
		dimen  = atoi( optarg);
	break;
	case 'f': 
		fn = strdup( optarg );
                break;
	case 'r':
	    break;
	case 't':
	    target=atoi( optarg );
	    break;
	default:
	    exit(0);
	}

    }


    FILE *src = fopen( fn, "r" );
    initboole(  dimen );
    initagldim( dimen );

    boole f;
    int item = 0;
    uint64_t size = 0;
    
    base = monomialBasis( 3, 3, 7 );
    base.size  = ( (size_t) 1 ) <<   base.mat.nbl;
    base.table = calloc( base.size, sizeof(char));

    while ((f = loadaglboolesize( src , &g, &size ))) {
		item++; 
    	        if ( ! target )
			printf( "item=%d degree=%d order=%ld\n", item, degree(f), size );
                if ( item == target ) break;
		aglfreeGroup( g );
		free(f);
    }
    if ( item != target ){
	    exit(0);
    }
    fclose(src);

    char name[64];
    sprintf( name, "rho/rho-%d.out", target );

    FILE* dst  = stdout; // = fopen( name, "w");
    if ( access( fn , F_OK) == 0) {
	    puts("using stdout");
      }   else dst   = fopen( name, "w");
    clock_t debut = clock();
 
    fprintf( dst, "item=%d degree=%d order=%ld\n", item, degree(f), size );
    panf( dst , f );
    assert( degree( f ) > 4  );

    lda    = aglBoundaryGroupAction( f, g, &base );
    size_t count = coverset( nbp );
		
    fprintf( dst, "\n#orbits=%ld", count );
    clock_t fin = clock();
    double temps = (double)(fin - debut) / CLOCKS_PER_SEC;
    fprintf(dst, "\nruntime : %.6f secondes\n", temps);
    char bname[64];
    sprintf( bname, "%d", item );
    saveBasename( bname , base );
    aglfreeGroup( g );
    aglVectorGroupFree(  lda );
    free( base.table );
    
    free( f );
    fclose(dst);
    return 0;
}
