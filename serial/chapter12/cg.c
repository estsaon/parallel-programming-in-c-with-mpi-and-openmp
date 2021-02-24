/*
 *   Conjugate Gradient Method
 */

#include <stdlib.h>
#include <stdio.h>

#define EPSILON 1.0e-10            /* Convergence criterion */
#define N       2                /* Size of linear system */

main (int argc, char *argv[])
{
   int i, j, it;                   /* Loop indices */
   double a[N][N];                 /* Solving Ax = b for x */
   double b[N];
   double x[N];
   double d[N];
   double g[N];                    /* Gradient */
   double denom1, denom2, num1, num2, s, tmpvec[N];  /* Temp variables */

   double dot_product (double a[N], double b[N]);
   void matrix_vector_product (double a[N][N], double b[N], double c[N]);

   /* Initialize a and b so that solution is x[i] = i */

   a[0][0] = 2.0;
   a[0][1] = 1.0;
   a[1][0] = 1.0;
   a[1][1] = 3.0;
   b[0] = 7.0;
   b[1] = 11.0;

   /* Initialize solution and gradient vectors */

   for (i = 0; i < N; i++) {
      d[i] = 0.0;
      x[i] = 0.0;
      g[i] = -b[i];
   }

   /* Algorithm converges in N or fewer iterations */

   for (it = 0; it < N; it++) {
printf ("x[0] = %6.3f   x[1] = %6.3f\n", x[0], x[1]);
      denom1 = dot_product (g, g);
      matrix_vector_product (a, x, g);
      for (i = 0; i < N; i++)
         g[i] -= b[i];
      num1 = dot_product (g, g);

      /* When g is sufficiently close to 0, time to halt */

      if (num1 < EPSILON) break;

      for (i = 0; i < N; i++)
         d[i] = -g[i] + (num1/denom1) * d[i];
      num2 = dot_product (d, g);
      matrix_vector_product (a, d, tmpvec);
      denom2 = dot_product (d, tmpvec);
      s = -num2 / denom2;
      for (i = 0; i < N; i++)
         x[i] += s * d[i];
   }

   printf ("Terminated after %d iterations\n", it+1);

   for (i = 0; i < N; i++)
      printf ("x[%d] = %8.4f\n", i, x[i]);
}

/*
 *   This function returns the dot product (inner product) of two vectors
 *    of length N
 */

double dot_product (double a[N], double b[N])
{
   int i;
   double answer;

   answer = 0.0;
   for (i = 0; i < N; i++)
      answer += a[i] * b[i];
   return answer;
}

/*
 *   This function computes the product of matrix a and vector b and
 *   stores the result in vector c
 */

void matrix_vector_product (double a[N][N], double b[N], double c[N])
{
   int i, j;
   double tmp;

   for (i = 0; i < N; i++) {
      tmp = 0.0;
      for (j = 0; j < N; j++)
         tmp += a[i][j] * b[j];
      c[i] = tmp;
   }
}
