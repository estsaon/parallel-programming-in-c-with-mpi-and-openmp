/*
 *   This program generates a PostScript file showing how gaussian
 *   elimination fills in a sparse system.
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
   void print(int i, double a[N][N+1]);

   printf ("\%\%BoundingBox: %d %d %d %d\n", 0, 0,
      3*(9*SQUARE_SIZE + 8*SMALL_GAP) + 2*BIG_GAP,
      3*(9*SQUARE_SIZE + 8*SMALL_GAP) + 2*BIG_GAP);
   printf ("/blacksquare {2 copy moveto\n");
   printf ("2 copy 2 1 roll %d add 2 1 roll lineto\n", SQUARE_SIZE);
   printf ("2 copy %d add 2 1 roll %d add 2 1 roll lineto\n", SQUARE_SIZE, SQUARE_SIZE);
   printf ("%d add lineto closepath fill } def\n", SQUARE_SIZE);
   printf ("/whitesquare {2 copy moveto\n");
   printf ("2 copy 2 1 roll %d add 2 1 roll lineto\n", SQUARE_SIZE);
   printf ("2 copy %d add 2 1 roll %d add 2 1 roll lineto\n", SQUARE_SIZE, SQUARE_SIZE);
   printf ("%d add lineto closepath stroke } def\n", SQUARE_SIZE);
   printf ("/Times-Roman findfont 20 scalefont setfont\n");
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

print(i, a);

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


   return 0;
}

void print(int it, double a[N][N+1]) {
   int i;
   int j;
   int xoffset, yoffset;
   int block_size;
   int top_left_x, top_left_y, ylocal, xlocal;

   block_size = 9 * SQUARE_SIZE + 8 * SMALL_GAP;
   top_left_x = 0;
   top_left_y = 3 * block_size + 2 * BIG_GAP;
   xoffset = (it%3)*(block_size + BIG_GAP);
   yoffset = (it/3)*(block_size + BIG_GAP);
   printf ("%d %d moveto\n", top_left_x + xoffset,
      top_left_y - yoffset + TEXT_OFFSET);
   printf ("(%d) show\n", it);
   for (i = 0; i < N; i++) {
      ylocal = top_left_y - yoffset - i*(SQUARE_SIZE + SMALL_GAP);
      for (j = 0; j < N; j++) {
         xlocal = top_left_x + xoffset + j*(SQUARE_SIZE + SMALL_GAP);
         if (a[i][j] == 0.0) printf ("%d %d whitesquare\n", xlocal, ylocal);
         else printf ("%d %d blacksquare\n", xlocal, ylocal);
      }
   }
}
