/*
 *   Generate files needed for conjuage gradient algorithm
 *
 *   Programmer: Michael J. Quinn
 *
 *   Last modification: 5 February 1999
 */

#include <stdlib.h>
#include <stdio.h>

main (int argc, char *argv[])
{
   int i, j;                   /* Loop indices */
   double **a, *astorage;
   double *b;
   unsigned short xi[3];
   int n;
   FILE *fp;

   /* Initialize a and b so that solution is x[i] = i */

   n = atoi(argv[1]);
   printf ("Will create a system of size %d\n", n);
   astorage = (double *) malloc (n * n * sizeof(double));
   a = (double **) malloc (n * sizeof (double *));
   for (i = 0; i < n; i++)
      a[i] = &astorage[i*n];
   b = (double *) malloc (n * sizeof(double));

   xi[0] = 4;
   xi[1] = 153;
   xi[2] = 17;
   for (i = 0; i < n; i++) {
      a[i][i] = 0.0;
      for (j = 0; j < n; j++) {
         if (i != j) a[i][j] = erand48(xi);
         a[i][i] += 2.0 * a[i][j];
      }
      b[i] = 0.0;
      for (j = 0; j < n; j++) {
         b[i] += a[i][j] * j;
      }
   }
   printf ("Will put matrix in file '%s'\n", argv[2]);
   fp = fopen (argv[2], "w");
   fwrite (&n, sizeof(int), 1, fp);
   fwrite (&n, sizeof(int), 1, fp);
   fwrite (astorage, sizeof(double), n*n, fp);
   fclose (fp);
   printf ("Will put vector in file '%s'\n", argv[3]);
   fp = fopen (argv[3], "w");
   fwrite (&n, sizeof(int), 1, fp);
   fwrite (b, sizeof(double), n, fp);
   fclose (fp);
}
