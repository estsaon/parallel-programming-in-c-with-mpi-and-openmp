/* Generate a square matrix of integers and write the values
   to a file. */
#include <stdio.h>
main (int argc, char * argv[]) {
   int i, j;
   int n;
   FILE *foutptr;
   int *a;
   int *ptr;

   printf ("argv[0] is '%s'\n", argv[0]);
   printf ("argv[1] is '%s'\n", argv[1]);
   printf ("argv[2] is '%s'\n", argv[2]);
   n = atoi (argv[1]);
   a = (int *) malloc ((n * n + 1)* sizeof(int));
   ptr = a;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) 
         if (i == j) *(ptr++) = 0;
         else *(ptr++) = (i + j) % 7;
   }
   foutptr = fopen (argv[2], "w");
   fwrite (&n, sizeof(int), 1, foutptr);
   fwrite (&n, sizeof(int), 1, foutptr);
   fwrite (a, sizeof(int), n*n, foutptr);
   fclose (foutptr);
}
