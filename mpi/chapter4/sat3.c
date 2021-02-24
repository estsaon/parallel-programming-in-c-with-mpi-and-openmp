/*
 *   Circuit Satisfiability, Version 3
 *
 *   This version of the program prints the total number of
 *   solutions and the execution time.
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 3 September 2002
 */

#include "mpi.h"
#include <stdio.h>

int main (int argc, char *argv[]) {
   int count;            /* Solutions found by this proc */
   double elapsed_time;  /* Time to find, count solutions */
   int global_count;     /* Total number of solutions */
   int i;
   int id;               /* Process rank */
   int p;                /* Number of processes */
   int check_circuit (int, int);

   MPI_Init (&argc, &argv);

   /* Start timer */
   MPI_Barrier (MPI_COMM_WORLD);
   elapsed_time = - MPI_Wtime();

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

   count = 0;
   for (i = id; i < 65536; i += p)
      count += check_circuit (id, i);

   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0,
      MPI_COMM_WORLD); 

   /* Stop timer */
   elapsed_time += MPI_Wtime();

   if (!id) {
      printf ("Execution time %8.6f\n", elapsed_time);
      fflush (stdout);
   }
   MPI_Finalize();
   if (!id) printf ("There are %d different solutions\n",
      global_count);
   return 0;
}

/* Return 1 if 'i'th bit of 'n' is 1; 0 otherwise */
#define EXTRACT_BIT(n,i) ((n&(1<<i))?1:0)

int check_circuit (int id, int z) {
   int v[16];        /* Each element is a bit of z */
   int i;

   for (i = 0; i < 16; i++) v[i] = EXTRACT_BIT(z,i);
   if ((v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
      && (!v[3] || !v[4]) && (v[4] || !v[5])
      && (v[5] || !v[6]) && (v[5] || v[6])
      && (v[6] || !v[15]) && (v[7] || !v[8])
      && (!v[7] || !v[13]) && (v[8] || v[9])
      && (v[8] || !v[9]) && (!v[9] || !v[10])
      && (v[9] || v[11]) && (v[10] || v[11])
      && (v[12] || v[13]) && (v[13] || !v[14])
      && (v[14] || v[15])) {
      printf ("%d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n", id,
         v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
         v[10],v[11],v[12],v[13],v[14],v[15]);
      fflush (stdout);
      return 1;
   } else return 0;
}
