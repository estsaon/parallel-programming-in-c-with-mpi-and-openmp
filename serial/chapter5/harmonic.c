/*
 *   C program that solves Exercise 5.11.
 *
 *   Programmed by Michael J. Quinn
 */

int main (int argc, char *argv[]) {

   int d;
   int i, j;
   int num;
   int carry;
   int n;
   int *v;
   double answer;

   n = atoi(argv[1]);
   d = atoi(argv[2]);
   printf ("Computing S%d to %d digits of precision\n", n, d);
   answer = 0.0;
   for (i = 1; i <= n; i++)
      answer += (1.0 / (double) i);
   printf ("Estimate of answer: %f\n\n", answer);
   v = (int *) malloc ((d+1) * sizeof(int));
   v[0] = 1;
   for (i = 1; i <= d; i++) v[i] = 0;
   for (i = 2; i <= n; i++) {
      num = 10;
      for (j = 1; j <= d; j++) {
         v[j] += (num / i);
         num = (num - i*(num/i))*10;
      }
      if (num/i >= 5) v[d]++;
   }
   carry = 0;
   for (i = d; i >= 1; i--) {
      v[i] += carry;
      carry = v[i] / 10;
      v[i] = v[i] % 10;
   }
   v[0] += carry;
   printf ("%d.", v[0]);
   for (i = 1; i <= d; i++) printf ("%d", v[i]);
   printf ("\n");
}
