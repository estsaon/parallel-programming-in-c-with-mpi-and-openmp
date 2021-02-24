/*
 *   Conjugate Gradient Method
 *
 *   This hybrid MPI/OpenMP program uses the conjugate gradient method to
 *   solve a postive definite system of linear equations.
 *
 *   Programmer: Michael J. Quinn
 *
 *   Last modification: 28 September 2006
 */

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#include "MyMPI.h"

main (int argc, char *argv[])
{
   double **a;                /* Solving Ax = b for x */
   double *astorage;          /* Holds elements of A */
   double *b;                 /* Constant vector */
   double *x;                 /* Solution vector */
   int     p;                 /* MPI Processes */
   int     id;                /* Process rank */
   int     m;                 /* Rows in A */
   int     n;                 /* Columns in A */
   int     n1;                /* Elements in b */

   /* Initialize a and b so that solution is x[i] = i */

   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   read_row_striped_matrix (argv[1], (void *) &a, (void *) &astorage,
      MPI_DOUBLE, &m, &n, MPI_COMM_WORLD);
   read_replicated_vector (argv[2], (void **) &b,
      MPI_DOUBLE, &n1, MPI_COMM_WORLD);
   if ((m != n) || (n != n1)) {
      if (!id)
         printf ("Incompatible dimensions (%d x %d) x (%d)\n", m, n, n1);
   } else {
      x = (double *) malloc (n * sizeof(double));
      cg (p, id, a, b, x, n);
      print_replicated_vector (x, MPI_DOUBLE, n, MPI_COMM_WORLD);
   }
   MPI_Finalize();
}

#define EPSILON 1.0e-10       /* Convergence criterion */

double *piece;                /* Temporary space used in mv product */

/* Conjugate gradient method solves ax = b for x */

cg (int p, int id, double **a, double *b, double *x, int n )
{
   int     i, it;                /* Loop indices */
   double *d;
   double *g;                    /* Gradient vector */
   double  denom1, denom2, num1,
           num2, s, *tmpvec;     /* Temporaries */

   double dot_product (double *, double *, int);
   void matrix_vector_product (int, int, int, double **,
                               double *, double *);

   /* Initialize solution and gradient vectors */

   d = (double *) malloc (n * sizeof(double));
   g = (double *) malloc (n * sizeof(double));
   tmpvec = (double *) malloc (n * sizeof(double));
   piece = (double *) malloc (BLOCK_SIZE(id,p,n) *
                              sizeof(double));
   for (i = 0; i < n; i++) {
      d[i] = x[i] = 0.0;
      g[i] = -b[i];
   }

   /* Algorithm converges in n or fewer iterations */

   for (it = 0; it < n; it++) {
      denom1 = dot_product (g, g, n);
      matrix_vector_product (id, p, n, a, x, g);
      for (i = 0; i < n; i++)
         g[i] -= b[i];
      num1 = dot_product (g, g, n);

      /* When g is sufficiently close to 0, it is time to halt */
      if (num1 < EPSILON) break;

      for (i = 0; i < n; i++)
         d[i] = -g[i] + (num1/denom1) * d[i];
      num2 = dot_product (d, g, n);
      matrix_vector_product (id, p, n, a, d, tmpvec);
      denom2 = dot_product (d, tmpvec, n);
      s = -num2 / denom2;
      for (i = 0; i < n; i++) x[i] += s * d[i];
   }
}

/*
 *   Return the dot product of two vectors
 */

double dot_product (double *a, double *b, int n)
{
   int i;
   double answer;

   answer = 0.0;
   for (i = 0; i < n; i++)
      answer += a[i] * b[i];
   return answer;
}

/*
 *   Compute the product of matrix a and vector b and
 *   store the result in vector c
 */

void matrix_vector_product (int id, int p, int n,
                            double **a, double *b, double *c)
{
   int    i, j;
   double tmp;       /* Accumulates sum */

   for (i = 0; i < BLOCK_SIZE(id,p,n); i++) {
      tmp = 0.0;
      for (j = 0; j < n; j++)
         tmp += a[i][j] * b[j];
      piece[i] = tmp;
   }
   replicate_block_vector (piece, n, (void *) c, MPI_DOUBLE, MPI_COMM_WORLD);
}
