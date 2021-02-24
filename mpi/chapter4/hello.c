#include "mpi.h"
#include <stdio.h>

int main (int argc, char *argv[]) {
   int id;               /* Process rank */

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);

   printf ("hello, world, from process %d\n", id);
   fflush (stdout);
   MPI_Finalize();
   return 0;
}
