/*
 *   Row pivot version of gaussian elimination.
 *
 *   This C program solves a dense system of linear equations Ax = b using
 *   gaussian elimination with partial pivoting of rows followed by back
 *   substitution. It inputs the system of equations from a file. The
 *   contents of the file are in this order:
 *
 *   A    A    ... A      b  A    A   ... A      b  ... A        b
 *    0,0  0,1      0,n-1  0  1,0  1,1     1,n-1  1      n-1,n-1  n-1
 *
 *   Written by Michael J. Quinn
 *
 *   Last modification: 2 February 2003
 */

#include <stdio.h>
#include <math.h>

int main (int argc, char *argv[]) {
   double **a;        /* Stores matrix A in first n columns and vector b in
                         last column */
   int     i, j, k;
   int    *loc;       /* loc[i] = r means during iteration i row r was pivot */
   double  magnitude; /* Largest value found so far in column i */
   int     n;         /* Size of system */
   int     picked;    /* Row chosen as pivot row */
   double  t;         /* Temporary */
   double *x;         /* Solution vector */
   int     tmp;
   int it;

   void read_linear_system (char *, double ***, int *);

   /* Make sure user put file name on command line */

   if (argc != 2) {
      printf ("Command line syntax: %s <file-name>\n", argv[0]);
      exit (1);
   }

   /* Read the system of equations from a file */

   read_linear_system (argv[1], &a, &n);

   /* Set aside memory for other useful arrays */

   loc = (int *) malloc (n * sizeof (int));
   for (i = 0; i < n; i++) loc[i] = i;
/*
printf ("loc: "); for (it = 0; it < n; it++) printf ("%d ", loc[it]);putchar('\n');
*/
   x = (double *) malloc (n * sizeof (double));

   /* Gaussian elimination. Transform matrix A from a dense matrix into
      an upper triangular matrix. */

   for (i = 0; i < n; i++) {

      /* Find the pivot col in row i */

      magnitude = 0.0;
      for (j = i; j < n; j++) {
         if (fabs(a[loc[j]][i]) > magnitude) {
            magnitude = fabs(a[loc[j]][i]);
            picked = j;
         }
      }
      tmp = loc[i];
      loc[i] = loc[picked];
      loc[picked] = tmp;

      /* We don't actually switch the pivot col with col "i". We just
         mark it where it is and remember its location by making
         appropriate entries in arrays "pivot" and "loc". */

      /* Reduce all unmarked rows */

      for (j = i+1; j < n; j++) {
         t = a[loc[j]][i] / a[loc[i]][i];
         for (k = i; k < n+1; k++) {
            a[loc[j]][k] -= (a[loc[i]][k] * t);
         }
      }
   }

   /* Back substitution. Reduce matrix A from an upper triangular matrix to
      a diagonal matrix. */

   for (i = n-1; i >= 0; i--) {

      /* Find value of x[i] */

      x[i] = a[loc[i]][n] / a[loc[i]][i];

      /* Now reduce all rows "above" row i. */
         
      for (j = 0; j < n; j++)
         a[loc[j]][n] -= (x[i] * a[loc[j]][i]);
   }
   
   /* Print solution */

   for (i = 0; i < n; i++)
      printf ("x[%d] = %10.4f\n", i, x[i]);
   return 0;
}

/*
 *    Function 'read_linear_system' inputs the system of 'n' linear equations
 *    from file 'fname' and returns it through the matrix pointer 'a'.
 *
 *    The file consists of two integers (r and c) followed by r * c doubles.
 *    This program expects a file of dimension r x (r+1), so if c != r+1
 *    the program terminates. The matrix elements are assumed to be in row
 *    major order. For this application, the last element of row i is
 *    actually the ith element of b.
 */

void read_linear_system (char *fname, double ***a, int *n) {

   int     c;         /* Columns */
   FILE   *fp;        /* File pointer */
   int     i;
   int     r;         /* Rows */
   double *storage;   /* Heap location where matrix elements stored */

   if ((fp = fopen (fname, "r")) == NULL) {
      printf ("Cannot open file '%s'\n", fname);
      exit (1);
   } 

   /* Get dimensions of matrix and make sure they meet expectations */

   fread (&r, sizeof(int), 1, fp);
   fread (&c, sizeof(int), 1, fp);
   if (c != (r+1)) {
      printf ("Incorrect matrix dimensions: %d x %d\n", r, c);
      fclose (fp);
      exit (1);
   }
   *n = r;

   /* Read matrix from file */

   if ((storage = (double *) malloc (r * c * sizeof(double))) == NULL) {
      printf ("Allocation of matrix storage failed\n");
      exit (1);
   }
   fread (storage, sizeof(double), r*c, fp);
   fclose (fp);

   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &(storage[i*c]);
}
