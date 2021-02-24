/*
 *   Matrix-vector multiplication, Version 2
 *
 *   This program multiplies a matrix and a vector.
 *   The matrix and vector are input from files.
 *   The answer is printed to standard output.
 *
 *   Data distribution of matrix: columnwise block striped
 *   Data distribution of vector: block
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 9 September 2002
 */

#include <stdio.h>
#include <mpi.h>
#include "../MyMPI.h"

/* Change these two definitions when the matrix and vector
   element types change */

typedef double dtype;
#define mpitype MPI_DOUBLE

int main (int argc, char *argv[]) {
   dtype **a;             /* The first factor, a matrix */
   dtype  *b;             /* The second factor, a vector */
   dtype  *c;             /* The product, a vector */
   dtype  *c_part_out;    /* Partial sums, sent */
   dtype  *c_part_in;     /* Partial sums, received */
   int    *cnt_out;       /* Elements sent to each proc */
   int    *cnt_in;        /* Elements received per proc */
   int    *disp_out;      /* Indices of sent elements */
   int    *disp_in;       /* Indices of received elements */
   int     i, j;          /* Loop indices */
   int     id;            /* Process ID number */
   int     local_els;     /* Cols of 'a' and elements of 'b'
                             held by this process */
   int     m;             /* Rows in the matrix */
   int     n;             /* Columns in the matrix */
   int     nprime;        /* Size of the vector */
   int     p;             /* Number of processes */
   double    max_seconds;
   double    seconds;     /* Elapsed time */
   dtype  *storage;       /* This process's portion of 'a' */

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

   read_col_striped_matrix (argv[1], (void ***) &a,
      (void **) &storage, mpitype, &m, &n, MPI_COMM_WORLD);
   print_col_striped_matrix ((void **) a, mpitype, m, n, 
      MPI_COMM_WORLD);
   read_block_vector (argv[2], (void **) &b, mpitype,
      &nprime, MPI_COMM_WORLD);
   print_block_vector ((void *) b, mpitype, nprime,
      MPI_COMM_WORLD);

   /* Each process multiplies its columns of 'a' and vector
      'b', resulting in a partial sum of product 'c'. */

   c_part_out = (dtype *) my_malloc (id, n * sizeof(dtype));
   local_els = BLOCK_SIZE(id,p,n);
   MPI_Barrier (MPI_COMM_WORLD);
   seconds = -MPI_Wtime();
   for (i = 0; i < n; i++) {
      c_part_out[i] = 0.0;
      for (j = 0; j < local_els; j++)
         c_part_out[i] += a[i][j] * b[j];
   }

   create_mixed_xfer_arrays (id, p, n, &cnt_out, &disp_out);
   create_uniform_xfer_arrays (id, p, n, &cnt_in, &disp_in);
   c_part_in =
      (dtype*) my_malloc (id, p*local_els*sizeof(dtype));
   MPI_Alltoallv (c_part_out, cnt_out, disp_out, mpitype,
      c_part_in, cnt_in, disp_in, mpitype, MPI_COMM_WORLD);

   c = (dtype*) my_malloc (id, local_els * sizeof(dtype));
   for (i = 0; i < local_els; i++) {
      c[i] = 0.0;
      for (j = 0; j < p; j++)
         c[i] += c_part_in[i + j*local_els];
   }
   MPI_Barrier (MPI_COMM_WORLD);
   seconds += MPI_Wtime();
   MPI_Allreduce (&seconds, &max_seconds, 1, mpitype, MPI_MAX,
      MPI_COMM_WORLD);
   if (!id) {
      printf ("MV3) N = %d, Processes = %d, Time = %12.6f sec,",
         n, p, max_seconds);
      printf ("Mflop = %6.2f\n", 2*n*n/(1000000.0*max_seconds));
   }
   print_block_vector ((void *) c, mpitype, n, MPI_COMM_WORLD);
   MPI_Finalize();
   return 0;
}
