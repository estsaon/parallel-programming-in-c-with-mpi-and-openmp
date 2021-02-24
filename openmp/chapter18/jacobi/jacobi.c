/*
 *   MPI/OpenMP program that solves the steady-state temperature
 *   distribution problem using the Jacobi method.
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 17 July 2003
 */

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "../../MyMPI.h"

#define N 10                        /* Grid size */
#define EPSILON 0.001               /* Termination condition */

int main (int argc, char *argv[])
{
   double **u;                      /* Previous temperatures */
   double **w;                      /* New temperatures */
   int     my_rows;                 /* Rows controlled by this process */
   int     p;                       /* Number of processes */
   int     id;                      /* Process rank */
   int     its;                     /* Iterations to converge */
   double elapsed;                  /* Execution time */

   void initialize_arrays (int, int, int *, double ***,
                           double ***);
   int find_steady_state (int, int, int, double **, double **);
   void print_solution (int, int, double **);

   MPI_Init (&argc, &argv);
   elapsed = -MPI_Wtime();
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   omp_set_num_threads (atoi(argv[1]));

   initialize_arrays (p, id, &my_rows, &u, &w);

   its = find_steady_state (p, id, my_rows, u, w);
   elapsed += MPI_Wtime();
   print_solution (p, id, u);
   if (!id) {
      printf ("Converged after %d iterations\n", its);
      printf ("Elapsed time = %8.6f\n", elapsed);
   }
   MPI_Finalize();
}

/* Allocate two-dimensional array. See Section 6.3 of the book for
   more details. */

void allocate_2d_array (int r, int c, double ***a)
{
   double *storage;
   int     i;

   storage = (double *) malloc (r * c * sizeof(double));
   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &storage[i * c];
}

void initialize_arrays (int p, int id, int *my_rows,
                        double ***u, double ***w)
{
   int    i, j;
   double mean;

   if (p == 1)
      *my_rows = BLOCK_SIZE(id,p,N);
   else if ((id == 0) || (id == p-1))
      *my_rows = BLOCK_SIZE(id,p,N) + 1;
   else
      *my_rows = BLOCK_SIZE(id,p,N) + 2;
   allocate_2d_array (*my_rows, N, u);
   allocate_2d_array (*my_rows, N, w);

   /* Set boundary conditions */

   if (id == 0)
      for (j = 0; j < N; j++)
         (*u)[0][j] = 100.0;
   if (id == p-1)
      for (j = 0; j < N; j++)
         (*u)[*my_rows-1][j] = 0.0;
   for (i = 0; i < *my_rows; i++) {
      (*u)[i][0] = (*u)[i][N-1] = 100.0;
   }

   /* Set initial inner values */

   mean = 75.0;
   for (i = 1; i < *my_rows-1; i++)
      for (j = 1; j < N-1; j++) (*u)[i][j] = mean;
}

void print_solution (int p, int id, double **u)
{
   int i, j, k, kblocks;
   MPI_Status status;
   double **t;

   if (p == 1) {
      for (i = 0; i < BLOCK_SIZE(0,1,N); i++) {
         for (j = 0; j < N; j++)
            printf ("%6.2f ", u[i][j]);
         printf ("\n");
      }
      printf ("\n");
   } else if (!id) {
      allocate_2d_array (BLOCK_SIZE(p-1,p,N), N, &t);
      for (i = 0; i < BLOCK_SIZE(0,p,N); i++) {
         for (j = 0; j < N; j++)
            printf ("%6.2f ", u[i][j]);
         printf ("\n");
      }
      for (k = 1; k < p; k++) {
         kblocks = BLOCK_SIZE(k,p,N);
         MPI_Recv (t[0], kblocks*N, MPI_DOUBLE, k, 1, MPI_COMM_WORLD,
            &status);
         for (i = 0; i < kblocks; i++) {
            for (j = 0; j < N; j++)
               printf ("%6.2f ", t[i][j]);
            printf ("\n");
         }
      }
      printf ("\n");
   } else {
      MPI_Send (u[1], BLOCK_SIZE(id,p,N)*N, MPI_DOUBLE, 0, 1,
         MPI_COMM_WORLD);
   }
}

int find_steady_state (int p, int id, int my_rows,
                       double **u, double **w)
{
   double     diff;             /* Maximum difference on this process */
   double     global_diff;      /* Globally maximum temperature difference */
   int        i, j;
   int        its;              /* Iteration count */
   MPI_Status status;           /* Result of receive */
   double     tdiff;            /* Maximum difference on this thread */

   its = 0;
   for (;;) {

      /* Exchange rows for ghost buffers */

      if (id > 0)
         MPI_Send (u[1], N, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD);
      if (id < p-1) {
         MPI_Send (u[my_rows-2], N, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD);
         MPI_Recv (u[my_rows-1], N, MPI_DOUBLE, id+1, 0, MPI_COMM_WORLD,
            &status);
      }
      if (id > 0)
         MPI_Recv (u[0], N, MPI_DOUBLE, id-1, 0, MPI_COMM_WORLD, &status);

      /* Update temperatures */

      diff = 0.0;
#pragma omp parallel private (i, j, tdiff)
{
      tdiff = 0.0;
#pragma omp for
      for (i = 1; i < my_rows-1; i++)
         for (j = 1; j < N-1; j++) {
            w[i][j] = (u[i-1][j] + u[i+1][j] +
                       u[i][j-1] + u[i][j+1])/4.0;
            if (fabs(w[i][j] - u[i][j]) > tdiff)
               tdiff = fabs(w[i][j] - u[i][j]);
         }
#pragma omp for nowait
      for (i = 1; i < my_rows-1; i++)
         for (j = 1; j < N-1; j++)
            u[i][j] = w[i][j];
#pragma omp critical
      if (tdiff > diff) diff = tdiff;
}
      MPI_Allreduce (&diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX,
         MPI_COMM_WORLD);

      /* Terminate if temperatures have converged */

      if (global_diff <= EPSILON) break;

      its++;
   }
   return its;
}
