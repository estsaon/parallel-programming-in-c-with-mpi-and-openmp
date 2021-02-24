/*
 *   This program produces the "dislike" matrix used by the simulated
 *   annealing program that solves the roommate assignment problem.
 *
 *   Programmed by Michael J. Quinn
 */

#include <stdlib.h>

#define N 20

int main (int argc, char *argv[]) {
   int i, j;
   unsigned short xi[3];
   float a[N][N];

   xi[0] = 1;
   xi[1] = 2;
   xi[2] = 3;

   for (i = 0; i < N; i++) {
      a[i][i] = 0.0;
      for (j = i+1; j < N; j++) {
         a[i][j] = erand48(xi)*10.0;
         a[j][i] = a[i][j];
      }
   }
   for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
         printf (" %4.2f,", a[i][j]);
      printf ("\n");
   }
}
