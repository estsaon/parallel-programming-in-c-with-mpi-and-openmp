/* Generate an n-element double vector and write it to a file */

#include <stdio.h>
main (int argc, char * argv[]) {
   int i, j;
   int n;
   FILE *foutptr;
   double *a;
   double *ptr;

   printf ("argv[0] is '%s'\n", argv[0]);
   printf ("argv[1] is '%s'\n", argv[1]);
   printf ("argv[2] is '%s'\n", argv[2]);
   n = atoi (argv[1]);
   a = (double *) malloc (n * sizeof(double));
   ptr = a;
   for (i = 0; i < n; i++) {
      *(ptr++) = (double) (i) / (double) (n);
   }
   foutptr = fopen (argv[2], "w");
   fwrite (&n, sizeof(int), 1, foutptr);
   fwrite (a, sizeof(double), n, foutptr);
   fclose (foutptr);
}
