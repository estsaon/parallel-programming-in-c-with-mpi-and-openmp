/*
 *   Sieve of Eratosthenes
 *
 *   This MPI program computes the number of prime numbers
 *   less than N, where N is a command-line argument.
 *
 *   Enhancements:
 *      Only odd integers are represented
 *      Each process finds its own prime sieve values: there
 *         is no broadcast step
 *
 *   Programmer: Michael J. Quinn
 *
 *   Last modification: 6 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))

main (int argc, char *argv[])
{
   int    current_prime;     /* Sieve by this prime */
   double elapsed_time;      /* Stopwatch */
   int    id_num;            /* Process rank */
   int    n;                 /* Top integer to check */
   int    p;                 /* Number of processors */
   int    sqrt_n;            /* Top value for sieve primes */
   char  *small_primes;      /* 1's show primes to sqrt(n) */
   char  *primes;            /* Share of values 3..n */
   int    small_prime_count; /* Number of sieve primes */
   int   *small_prime_values;/* Array of sieving primes */
   int    index;             /* Sieving location */
   int    i, j, k;
   int    size;              /* Size of array 'primes' */
   int    prime_count;       /* Primes on this proc */
   int    low_proc_value;    /* Smallest int on this proc */
   int    high_proc_value;   /* Highest int on this proc */
   int    smaller_size;      /* Smaller block size */
   int    larger_size;
   int    num_larger_blocks;
   int blocks;
   int low;
   int high;
   int global_count;
   int small_prime_array_size;
   int els;
   

   MPI_Init (&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();
   n = atoi (argv[1]);
   MPI_Comm_rank (MPI_COMM_WORLD, &id_num);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   sqrt_n = (int) sqrt((double) n);
   small_prime_array_size = (sqrt_n - 1)/2;
   small_primes = (char *) malloc (small_prime_array_size);
   small_prime_count = seq_sieve (small_primes, small_prime_array_size, sqrt_n);
   small_prime_values = (int *) malloc (small_prime_count * sizeof(int));
   index = 0;
   for (i = 0; i < small_prime_array_size; i++)
      if (small_primes[i]) {
         small_prime_values[index++] = 2*i+3;
      }
   els = (n-1) / 2;
   smaller_size = els / p;
   larger_size = smaller_size + 1;
   num_larger_blocks = els % p;
   if (id_num < num_larger_blocks) size = larger_size;
   else size = smaller_size;
   low_proc_value =
      2*(id_num*smaller_size + MIN(id_num, num_larger_blocks)) + 3;
   high_proc_value = low_proc_value + 2 * (size-1);
   primes = (char *) malloc (size);


   for (j = 0; j < size; j++) primes[j] = 1;

   for (j = 0; j < small_prime_count; j++) {
      current_prime = small_prime_values[j];
      if (current_prime * current_prime > low_proc_value)
         index = (current_prime * current_prime - low_proc_value)/2;
      else {
         int r = low_proc_value % current_prime;
         if (!r) index = 0;
         else if ((current_prime - r) & 1)
            index = (2*current_prime - r)/2;
         else index = (current_prime - r)/2;
      }
      if (index >= size) break;
      for (k = index; k < size; k+= current_prime)
         primes[k] = 0;
   }

   prime_count = 0;
   for (j = 0; j < size; j++)
      if (primes[j]) prime_count++;
   MPI_Reduce (&prime_count, &global_count, 1, MPI_INT, MPI_SUM, 0,
      MPI_COMM_WORLD);
   if (!id_num) global_count++;   /* To account for only even prime, 2 */
   elapsed_time += MPI_Wtime();
   if (!id_num) {
      printf ("Total prime count is %d\n", global_count);
      printf ("SIEVE_ODD_NO_BCAST (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize();
   return 0;
}

seq_sieve (char *small_primes, int small_prime_array_size, int sqrt_n)
{
   int i, j;
   int prime_index;
   int prime_value;
   int count;

   /* small_primes[i] represents integer 2i+3 */

   for (i = 0; i < small_prime_array_size; i++) small_primes[i] = 1;
   prime_index = 0;
   prime_value = 3;
   while (prime_value * prime_value <= sqrt_n) {
      j = prime_value * prime_value / 2 - 1;
      while (j < small_prime_array_size) {
         small_primes[j] = 0;
         j += prime_value;
      }
      while (small_primes[++prime_index] == 0);
      prime_value = 2*prime_index + 3;
   }
   count = 0;
   for (i = 0; i < small_prime_array_size; i++)
      if (small_primes[i] == 1) {
         count++;
      }
   return count;
}
