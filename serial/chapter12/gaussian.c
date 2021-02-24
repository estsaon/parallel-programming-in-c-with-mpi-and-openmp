/*
 *   This C program solves a dense system of linear equations Ax = b using
 *   gaussian elimination with partial pivoting followed by back
 *   substitution. It inputs the system of equations from a file. The
 *   contents of the file are in this order:
 *
 *   A    A    ... A      b  A    A   ... A      b  ... A        b
 *    0,0  0,1      0,n-1  0  1,0  1,1     1,n-1  1      n-1,n-1  n-1
 *
 *   Written by Michael J. Quinn
 *
 *   Last modification: 23 February 2002
 */

#include <stdio.h>
#include <math.h>

#define N 9

int main (int argc, char *argv[]) {
   double a[N][N+1];        /* Stores matrix A in first n columns and vector b in
                         last column */
   int     i, j, k;
   int     loc[N];       /* loc[i] = r means during iteration i row r was pivot */
   double  magnitude; /* Largest value found so far in column i */
   int     marked[N];    /* Indicates which rows have been pivot rows */
   int     n;         /* Size of system */
   int     picked;    /* Row chosen as pivot row */
   int     pivot[N];     /* pivot[r] = i means row r was pivot row during
                         iteration i */
   double  t;         /* Temporary */
   double x[N];         /* Solution vector */
   int row, col;
   void print(double a[N][N+1]);

   n = N;
   for (i = 0; i < N; i++)
      for (j = 0; j < N+1; j++)
         a[i][j] = 0.0;
   for (i = 0; i < N; i++) {
      row = i / 3;
      col = i % 3;
      if (col > 0) a[i][i-1] = -0.25;
      if (col < 2) a[i][i+1] = -0.25;
      if (row > 0) a[i][i-3] = -0.25;
      if (row < 2) a[i][i+3] = -0.25;
      a[i][i] = 1.0;
      if (i < 3) a[i][N] = 25.0;
   }

   /* Initially no rows have been used as pivot rows */

   for (i = 0; i < n; i++) marked[i] = 0;

   /* Gaussian elimination. Transform matrix A from a dense matrix into
      an upper triangular matrix. */

   for (i = 0; i < n; i++) {

print(a);

      /* Find the pivot row */

      magnitude = 0.0;
      for (j = 0; j < n; j++) {
         if ((!marked[j]) && (fabs(a[j][i]) > magnitude)) {
            magnitude = fabs(a[j][i]);
            picked = j;
         }
      }

      /* We don't actually move the pivot row up to row "i". We just
         mark it where it is and remember its location by making
         appropriate entries in arrays "pivot" and "loc". */

      marked[picked] = 1;
      pivot[picked] = i;
      loc[i] = picked;

      /* Reduce all unmarked rows */

      for (j = 0; j < n; j++)
         if (!marked[j]) {
            t = a[j][i] / a[picked][i];
            for (k = i; k < n+1; k++)
               a[j][k] -= (a[picked][k] * t);
         }
   }

   /* Back substitution. Reduce matrix A from an upper triangular matrix to
      a diagonal matrix. */

   for (i = n-1; i >= 0; i--) {

print (a);
      /* Find value of x[i] */

      x[i] = a[loc[i]][n] / a[loc[i]][i];

      /* Now reduce all rows "above" row i. */
         
      for (j = 0; j < n; j++)
         if (pivot[j] < i) {
            a[j][n] -= (x[i] * a[j][i]);
            a[j][i] = 0.0;
         }
   }
   
   /* Print solution */

print (a);

   for (i = 0; i < n; i++)
      printf ("x[%d] = %10.4f\n", i, x[i]);
   return 0;
}

void print(double a[N][N+1]) {
   int i;
   int j;
   for (i = 0; i < N; i++) {
      for (j = 0; j < N+1; j++)
         if (a[i][j] == 0.0) printf (" . ");
         else printf (" X ");
/*
         printf (" %6.3f", a[i][j]);
*/
      printf ("\n");
   }
   printf ("\n");
}
