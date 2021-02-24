/*
 *   C program that uses the Monte Carlo method to solve the parking
 *   garage problem posed in Section 10.5.5.
 *
 *   Programmed by Michael J. Quinn
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define S 80
#define A 3.0

unsigned short xi[3];       /* random number seed */
int active = 0;
double g2;

double box_muller (void)
{
   double f, g1, r, v1, v2;   /* See pseudocode in book */

   if (active) {
      active = 0;
      return g2;
   } else {
      do {
         v1 = 2.0 * erand48(xi) - 1.0;
         v2 = 2.0 * erand48(xi) - 1.0;
         r = v1 * v1 + v2 * v2;
      } while ((r <= 0.0) || (r >= 1.0));
      f = sqrt(-2.0 * log(r)/r);
      g1 = f * v1;
      g2 = f * v2;
      active = 1;
      return g1;
   }
}

main (int argc, char *argv[])
{
   double limit;
   double t;
   double G[S];
   double u;
   unsigned short xi[3];
   int occupied;
   int turned_away;
   double area;
   double stay;
   int found;
   int i;
   int total;

   int parked;

   limit = (double) atoi(argv[1]);

   
   xi[0] = 0;
   xi[1] = 1;
   xi[2] = 2;
   turned_away = 0;
   total = 0;
   area = 0.0;
   for (i = 0; i < S; i++) G[i] = 0.0;
   t = 0.0;
   parked = 0;
   do {
      total++;
      found = 0;
      for (i = 0; i < S; i++) {
         if (G[i] <= t) {
            parked++;
            stay = 240.0 + 60.0 * box_muller();
            G[i] = t + stay;
            area += stay;
            found = 1;
            break;
         }      
      }
      if (!found) turned_away++;
      u = erand48(xi);
      t -= A * log(u);
   } while (t < limit);
   printf ("   Average capacity = %6.3f\n", area/t);
   printf ("   Average length of stay = %6.3f\n", area/parked);
   printf ("   Pct cars parked = %6.2f\n", (100.0 * parked)/total);
   printf ("   Pct cars turned away = %6.2f\n", (100.0*turned_away)/total);
}
