/*
 *   This program uses the Jacobi method to solve a sparse linear system.
 */

#include <stdio.h>
#include <math.h>

#define N 9
#define SQUARE_SIZE 10
#define SMALL_GAP 5
#define BIG_GAP 30
#define TEXT_OFFSET 15

int main (int argc, char *argv[]) {
   double a[N][N+1];        /* Stores matrix A in first n columns and vector b in
                         last column */
   int     i, j, k;
   double  sum; /* Largest value found so far in column i */
   double x[N];         /* Solution vector */
   double newx[N];
   int row, col;

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
      else a[i][N] = 0.0;
   }

   for (i = 0; i < N; i++) x[i] = 0.0;

   /* Jacobi iteration */

   for (i = 0; i < 100; i++) {
      printf ("Iteration %d:\n", i);
      for (j = 0; j < N; j++)
         printf ("   x[%d] = %7.3f\n", j, x[j]);
      printf ("\n");
      for (j = 0; j < N; j++) {
         sum = 0.0;
         for (k = 0; k < N; k++)
            if (k != j)
            sum += a[j][k]*x[k];
         newx[j] = (1.0/a[j][j])*(a[j][N] - sum);
/*
         x[j] = (1.0/a[j][j])*(a[j][N] - sum);
*/
      }
      for (j = 0; j < N; j++)
         x[j] = newx[j];
   }
   return 0;
}
