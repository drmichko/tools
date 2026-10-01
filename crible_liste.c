// crible_liste.c : pour chaque f0 d'une liste (degre 6, 6 variables, poids impair),
// decide s'il existe f1 tel que f = f0 + x7*(f0+f1) soit de degre 6 en 7 variables
// et de linearite <= LIN.  Methode : f1 = (p,q), p sur x6=0, q sur x6=1.
//   1) crible de p parmi 2^32 (code de Gray) : |P(a')| <= (beta(a',0)+beta(a',1))/2
//   2) jointure par seaux sur les 2 points a' les plus contraints, puis test complet
//      |P+Q| <= beta(a',0), |P-Q| <= beta(a',1), wt(p)+wt(q) impair.
// Compilation : gcc -O3 -march=native -fopenmp -o crible_liste crible_liste.c
// Usage       : ./crible_liste fichier_f0.txt [LIN=16] [ALL=0]
//   fichier : un f0 par ligne, hexa 64 bits (bit x = f0(x), x1 = bit 0)
//   ALL=1   : liste toutes les solutions de chaque f0 (sinon arret a la premiere)
// Sortie (une ligne par f0) :  <f0> OUI f1=<f1> | <f0> NON | <f0> TROP (indecis)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef struct { uint32_t f; int8_t w[32]; } item;
static int8_t R[32][32];
#ifndef MAXN
#define MAXN 30000000UL
#endif

static void traite(uint64_t f0, int LIN, int ALL) {
  int beta[64];
  for (int a = 0; a < 64; a++) {
    int s = 0;
    for (int x = 0; x < 64; x++) s += ((__builtin_popcount(a & x) + (f0 >> x)) & 1) ? -1 : 1;
    beta[a] = LIN - abs(s);
  }
  int8_t B[32];
  for (int a = 0; a < 32; a++) {
    if (beta[a] < 0 || beta[a + 32] < 0) {
      #pragma omp critical
      { printf("%016llx NON (linearite f0 > %d)\n", (unsigned long long)f0, LIN); fflush(stdout); }
      return;
    }
    B[a] = (beta[a] + beta[a + 32]) / 2;
  }
  size_t cap = 1 << 16, n = 0; item *L = malloc(cap * sizeof(item));
  int8_t w[32] = {0}; w[0] = 32; uint32_t p = 0; int trop = 0;
  double t0 = omp_get_wtime();
  for (uint64_t i = 0; ; i++) {
    if (i > 0) {
      int x = __builtin_ctzll(i);
      int s = ((p >> x) & 1) ? -1 : 1;
      for (int a = 0; a < 32; a++) w[a] -= s * R[x][a];
      p ^= 1u << x;
    }
    int bad = 0;
    for (int a = 0; a < 32; a++) bad |= (abs(w[a]) > B[a]);
    if (!bad) {
      if (n == MAXN) { trop = 1; break; }
      if (n == cap) { cap *= 2; L = realloc(L, cap * sizeof(item)); }
      L[n].f = p; memcpy(L[n].w, w, 32); n++;
    }
    if (i == 0xFFFFFFFFull) break;
  }
  double t1 = omp_get_wtime();
  if (trop) {
    #pragma omp critical
    { printf("%016llx TROP (>%lu survivants)\n", (unsigned long long)f0, MAXN); fflush(stdout); }
    free(L); return;
  }
  // points les plus contraints
  int a0 = 0, a1 = 1;
  for (int a = 0; a < 32; a++) if (B[a] < B[a0]) a0 = a;
  for (int a = 0; a < 32; a++) if (a != a0 && (a1 == a0 || B[a] < B[a1])) a1 = a;
  // seaux (P(a0),P(a1)) : index = ((P0+32)/2)*33 + (P1+32)/2
  size_t *start = calloc(33 * 33 + 1, sizeof(size_t)), *pos;
  uint32_t *idx = malloc((n ? n : 1) * sizeof(uint32_t));
  for (size_t i = 0; i < n; i++) start[((L[i].w[a0] + 32) / 2) * 33 + (L[i].w[a1] + 32) / 2 + 1]++;
  for (int k = 0; k < 33 * 33; k++) start[k + 1] += start[k];
  pos = malloc(33 * 33 * sizeof(size_t)); memcpy(pos, start, 33 * 33 * sizeof(size_t));
  for (size_t i = 0; i < n; i++) idx[pos[((L[i].w[a0] + 32) / 2) * 33 + (L[i].w[a1] + 32) / 2]++] = (uint32_t)i;
  int ord[32]; for (int a = 0; a < 32; a++) ord[a] = a;
  for (int i = 0; i < 32; i++) for (int j = i + 1; j < 32; j++)
    if (B[ord[j]] < B[ord[i]]) { int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }

  long nsol = 0; uint64_t first = 0; int stop = 0;
  for (size_t i = 0; i < n && !stop; i++) {
    int P0 = L[i].w[a0], P1 = L[i].w[a1], pi = __builtin_popcount(L[i].f) & 1;
    int lo0 = -beta[a0] - P0, hi0 = beta[a0] - P0; if (P0 - beta[a0 + 32] > lo0) lo0 = P0 - beta[a0 + 32]; if (P0 + beta[a0 + 32] < hi0) hi0 = P0 + beta[a0 + 32];
    int lo1 = -beta[a1] - P1, hi1 = beta[a1] - P1; if (P1 - beta[a1 + 32] > lo1) lo1 = P1 - beta[a1 + 32]; if (P1 + beta[a1 + 32] < hi1) hi1 = P1 + beta[a1 + 32];
    if (lo0 < -32) lo0 = -32; if (hi0 > 32) hi0 = 32; if (lo1 < -32) lo1 = -32; if (hi1 > 32) hi1 = 32;
    for (int q0 = lo0; q0 <= hi0 && !stop; q0++) {
      if (q0 & 1) continue;
      for (int q1 = lo1; q1 <= hi1 && !stop; q1++) {
        if (q1 & 1) continue;
        int b = ((q0 + 32) / 2) * 33 + (q1 + 32) / 2;
        for (size_t k = start[b]; k < start[b + 1]; k++) {
          size_t j = idx[k];
          if (((__builtin_popcount(L[j].f) & 1) ^ pi) == 0) continue;
          int ok = 1;
          for (int m = 0; m < 32 && ok; m++) {
            int a = ord[m], P = L[i].w[a], Q = L[j].w[a];
            if (abs(P + Q) > beta[a] || abs(P - Q) > beta[a + 32]) ok = 0;
          }
          if (ok) {
            uint64_t f1 = (uint64_t)L[i].f | ((uint64_t)L[j].f << 32);
            nsol++;
            if (ALL) {
              #pragma omp critical
              { printf("%016llx SOL f1=%016llx\n", (unsigned long long)f0, (unsigned long long)f1); fflush(stdout); }
            } else { first = f1; stop = 1; break; }
          }
        }
      }
    }
  }
  double t2 = omp_get_wtime();
  #pragma omp critical
  {
    if (ALL) printf("%016llx %s (%ld solutions) [crible %.0fs, %zu surv, jointure %.0fs]\n", (unsigned long long)f0, nsol ? "OUI" : "NON", nsol, t1 - t0, n, t2 - t1);
    else if (nsol) printf("%016llx OUI f1=%016llx [crible %.0fs, %zu surv, jointure %.0fs]\n", (unsigned long long)f0, (unsigned long long)first, t1 - t0, n, t2 - t1);
    else printf("%016llx NON [crible %.0fs, %zu surv, jointure %.0fs]\n", (unsigned long long)f0, t1 - t0, n, t2 - t1);
    fflush(stdout);
  }
  free(L); free(start); free(pos); free(idx);
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s fichier [LIN=16] [ALL=0]\n", argv[0]); return 1; }
  int LIN = argc > 2 ? atoi(argv[2]) : 16, ALL = argc > 3 ? atoi(argv[3]) : 0;
  for (int x = 0; x < 32; x++) for (int a = 0; a < 32; a++) R[x][a] = (__builtin_popcount(a & x) & 1) ? -2 : 2;
  FILE *fp = fopen(argv[1], "r"); if (!fp) { perror("fichier"); return 1; }
  uint64_t *F = NULL; size_t nf = 0, cf = 0; char line[256];
  while (fgets(line, sizeof line, fp)) {
    char *s = line; while (*s == ' ' || *s == '\t') s++;
    if (!*s || *s == '\n' || *s == '#') continue;
    uint64_t f = strtoull(s, 0, 16);
    if (__builtin_popcountll(f) % 2 == 0) { fprintf(stderr, "ignore (poids pair): %s", line); continue; }
    if (nf == cf) { cf = cf ? 2 * cf : 1024; F = realloc(F, cf * sizeof *F); }
    F[nf++] = f;
  }
  fclose(fp);
  fprintf(stderr, "%zu fonctions, %d threads\n", nf, omp_get_max_threads());
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < nf; i++) traite(F[i], LIN, ALL);
  return 0;
}
