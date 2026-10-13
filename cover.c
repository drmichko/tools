// ============================================================
// Cover set de B(3,3,7) sous Stab(h)
// Convention : T[i] = 0 (non visité), 1 (cover set), 2 (couvert)
// ============================================================

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>

// ============================================================
// Configuration
// ============================================================

#define N_ELEMENTS   (1ULL << 35)   // 2^35 = 34 359 738 368
#define NB_THREADS   72
#define FILE_TAILLE  (1 << 20)      // 1M éléments par file locale

// Valeurs de T
#define NON_VISITE   0
#define COVER_SET    1
#define COUVERT      2

// ============================================================
// Structures globales
// ============================================================

static uint8_t *T;                          // Table de 2^35 bytes
static atomic_uint_fast64_t curseur;        // Curseur global atomique


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
        if (i >= N_ELEMENTS) {
            return UINT64_MAX;
        }
        if (T[i] == NON_VISITE) {
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
    
    while (debut < fin) {
        uint64_t u = file[debut++];
        
        for (int g = 0; g < NB_GENERATEURS; g++) {
            uint64_t v = apply_generator(g, u);
            
            // Tenter de marquer v comme couvert
            if (cas_byte(&T[v], NON_VISITE, COUVERT)) {
                // v est nouveau : l'ajouter à la file
                if (fin < FILE_TAILLE) {
                    file[fin++] = v;
                } else {
                    // File pleine : traiter v immédiatement
                    // (ou augmenter la taille de la file)
                    // Ici, on abandonne pour simplifier
                    // En pratique, utiliser une file dynamique
                }
            }
            // Sinon : déjà marqué (COVER_SET ou COUVERT), ne rien faire
        }
        
        // Réinitialiser la file si elle est vide
        if (debut == fin) {
            debut = 0;
            fin = 0;
        }
    }
    
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
        if (!cas_byte(&T[x], NON_VISITE, COVER_SET)) {
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
    
    for (uint64_t i = 0; i < N_ELEMENTS; i++) {
        switch (T[i]) {
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

static void charger_generateurs(void) {
    // Ici, il faut charger les générateurs de Stab(h)
    // Pour l'exemple, on met un générateur trivial (identité)
    // À remplacer par le vrai calcul
    
    NB_GENERATEURS = 1;
    A_g = calloc(NB_GENERATEURS * 35, sizeof(uint64_t));
    b_g = calloc(NB_GENERATEURS, sizeof(uint64_t));
    
    // Identité : A[i] = 1ULL << i, b = 0
    for (int i = 0; i < 35; i++) {
        A_g[i] = 1ULL << i;
    }
    b_g[0] = 0;
}

// ============================================================
// Programme principal
// ============================================================

int main(void) {
    // 1. Allouer la table T (32 Go)
    printf("Allocation de T (%llu bytes = %.2f Go)...\n",
           (unsigned long long)N_ELEMENTS,
           (double)N_ELEMENTS / (1024.0 * 1024.0 * 1024.0));
    
    T = aligned_alloc(64, N_ELEMENTS);
    if (!T) {
        fprintf(stderr, "Erreur d'allocation de T\n");
        return 1;
    }
    memset(T, NON_VISITE, N_ELEMENTS);
    
    // 2. Initialiser le curseur
    atomic_init(&curseur, 0);
    
    // 3. Charger les générateurs
    charger_generateurs();
    
    // 4. Lancer les threads
    printf("Lancement de %d threads...\n", NB_THREADS);
    pthread_t threads[NB_THREADS];
    
    for (int i = 0; i < NB_THREADS; i++) {
        if (pthread_create(&threads[i], NULL, thread_worker, NULL) != 0) {
            fprintf(stderr, "Erreur pthread_create\n");
            return 1;
        }
    }
    
    // 5. Attendre la fin
    for (int i = 0; i < NB_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    
    // 6. Extraire le cover set
    extraire_cover_set();
    
    // 7. Libérer
    free(T);
    free(A_g);
    free(b_g);
    
    return 0;
}
