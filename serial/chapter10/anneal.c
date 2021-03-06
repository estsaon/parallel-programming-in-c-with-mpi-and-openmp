/*
 *   This C program uses simulated annealing to solve the roommate assignment
 *   problem of Section 10.5.4.
 *
 *   Programmed by Michael J. Quinn
 */

#include <stdlib.h>
#include <math.h>

#define N 20
#define C 0.999
#define T_START 10.0

double dislike[N][N] = {
 0.00, 4.42, 2.63, 6.54, 4.21, 0.23, 7.65, 9.92, 0.53, 1.17,
 4.38, 6.95, 8.10, 1.56, 5.45, 8.55, 2.60, 2.53, 8.03, 5.16,
 4.42, 0.00, 9.14, 9.47, 9.16, 6.98, 0.74, 1.92, 0.70, 5.15,
 3.01, 0.72, 6.38, 7.12, 9.24, 5.70, 1.11, 4.14, 3.62, 1.42,
 2.63, 9.14, 0.00, 2.50, 5.38, 0.55, 6.74, 0.75, 8.18, 5.46,
 1.46, 6.86, 5.10, 7.69, 9.47, 1.41, 6.69, 7.54, 7.66, 0.32,
 6.54, 9.47, 2.50, 0.00, 2.80, 3.44, 4.76, 6.57, 1.35, 8.74,
 6.67, 3.72, 9.58, 9.65, 9.78, 0.68, 1.85, 4.83, 0.76, 4.97,
 4.21, 9.16, 5.38, 2.80, 0.00, 7.84, 4.21, 8.69, 3.43, 3.41,
 7.90, 0.30, 2.86, 1.79, 1.55, 2.64, 5.81, 8.76, 7.66, 4.58,
 0.23, 6.98, 0.55, 3.44, 7.84, 0.00, 7.27, 9.87, 4.80, 3.54,
 0.92, 0.32, 2.63, 8.71, 8.60, 6.38, 4.43, 4.54, 5.23, 6.64,
 7.65, 0.74, 6.74, 4.76, 4.21, 7.27, 0.00, 6.23, 8.14, 4.34,
 9.47, 8.50, 1.66, 5.62, 0.96, 9.07, 3.74, 6.74, 2.78, 4.36,
 9.92, 1.92, 0.75, 6.57, 8.69, 9.87, 6.23, 0.00, 5.24, 2.99,
 8.60, 4.27, 0.05, 3.85, 8.33, 0.22, 6.84, 6.08, 4.26, 6.90,
 0.53, 0.70, 8.18, 1.35, 3.43, 4.80, 8.14, 5.24, 0.00, 6.08,
 5.23, 8.95, 7.78, 7.13, 7.96, 0.05, 2.29, 3.55, 5.21, 9.86,
 1.17, 5.15, 5.46, 8.74, 3.41, 3.54, 4.34, 2.99, 6.08, 0.00,
 3.65, 3.01, 3.41, 9.64, 0.19, 9.85, 8.86, 9.79, 8.63, 7.63,
 4.38, 3.01, 1.46, 6.67, 7.90, 0.92, 9.47, 8.60, 5.23, 3.65,
 0.00, 7.35, 4.76, 3.83, 3.98, 0.08, 4.33, 8.08, 5.21, 3.40,
 6.95, 0.72, 6.86, 3.72, 0.30, 0.32, 8.50, 4.27, 8.95, 3.01,
 7.35, 0.00, 5.11, 3.21, 1.83, 7.06, 5.54, 6.32, 3.86, 6.91,
 8.10, 6.38, 5.10, 9.58, 2.86, 2.63, 1.66, 0.05, 7.78, 3.41,
 4.76, 5.11, 0.00, 4.30, 1.82, 3.80, 7.64, 4.74, 0.48, 3.27,
 1.56, 7.12, 7.69, 9.65, 1.79, 8.71, 5.62, 3.85, 7.13, 9.64,
 3.83, 3.21, 4.30, 0.00, 8.40, 0.10, 4.91, 1.45, 5.98, 7.40,
 5.45, 9.24, 9.47, 9.78, 1.55, 8.60, 0.96, 8.33, 7.96, 0.19,
 3.98, 1.83, 1.82, 8.40, 0.00, 6.44, 1.06, 2.25, 9.04, 5.05,
 8.55, 5.70, 1.41, 0.68, 2.64, 6.38, 9.07, 0.22, 0.05, 9.85,
 0.08, 7.06, 3.80, 0.10, 6.44, 0.00, 2.13, 2.21, 2.37, 0.91,
 2.60, 1.11, 6.69, 1.85, 5.81, 4.43, 3.74, 6.84, 2.29, 8.86,
 4.33, 5.54, 7.64, 4.91, 1.06, 2.13, 0.00, 1.04, 9.82, 2.25,
 2.53, 4.14, 7.54, 4.83, 8.76, 4.54, 6.74, 6.08, 3.55, 9.79,
 8.08, 6.32, 4.74, 1.45, 2.25, 2.21, 1.04, 0.00, 7.95, 9.86,
 8.03, 3.62, 7.66, 0.76, 7.66, 5.23, 2.78, 4.26, 5.21, 8.63,
 5.21, 3.86, 0.48, 5.98, 9.04, 2.37, 9.82, 7.95, 0.00, 9.17,
 5.16, 1.42, 0.32, 4.97, 4.58, 6.64, 4.36, 6.90, 9.86, 7.63,
 3.40, 6.91, 3.27, 7.40, 5.05, 0.91, 2.25, 9.86, 9.17, 0.00
};

int assignment[N];
int try[N];

int main (int argc, char *argv[]) {
   int i, j, k;
   int flip;
   double total_dislike;
   double T;
   double new_dislike;
   int cand1, cand2;
   unsigned short xi[3];
   int tmp;
   int no_change;

   xi[0] = atoi(argv[1]);
   xi[1] = atoi(argv[2]);
   xi[2] = atoi(argv[3]);

   for (i = 0; i < N; i++) assignment[i] = i/2;

   for (i = 0; i < N-1; i++) {
      flip = erand48(xi)*(N-i);
      tmp = assignment[flip];
      assignment[flip] = assignment[N-i-1];
      assignment[N-i-1] = tmp;
   }
   total_dislike = 0.0;
   for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
         if (assignment[i] == assignment[j]) total_dislike += dislike[i][j];
   T = T_START;
   no_change = 0;
   for (i = 0; i < 10000; i++) {
      if (++no_change >= 1000) break;
      for (j = 0; j < N; j++) try[j] = assignment[j];
      do {
         cand1 = erand48(xi) * N;
         cand2 = erand48(xi) * N;
      } while (assignment[cand1] == assignment[cand2]);
      tmp = try[cand1];
      try[cand1] = try[cand2];
      try[cand2] = tmp;
      new_dislike = 0.0;
      for (k = 0; k < N; k++)
         for (j = 0; j < N; j++)
            if (try[k] == try[j]) new_dislike += dislike[k][j];
      if ((new_dislike < total_dislike) ||
          (erand48(xi) <= exp((double)((total_dislike-new_dislike)/T)))) {
         for (j = 0; j < N; j++) assignment[j] = try[j];
         total_dislike = new_dislike;
         no_change = 0;
      }
      T = C * T;
   }
   for (i = 0; i < N; i++)
      for (j = i+1; j < N; j++)
         if (assignment[i] == assignment[j])
            printf ("Persons %2d and %2d assigned to the same room %6.3f\n",
            i, j, dislike[i][j]);
}
