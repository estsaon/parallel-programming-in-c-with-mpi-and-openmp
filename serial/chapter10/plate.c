/*
 *   C program that uses the Monte Carlo method to solve the neutron
 *   transport problem posed in Section 10.5.1.
 *
 *   Programmed by Michael J. Quinn
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define C  2.0
#define CS 1.5
#define CC (C - CS)
#define PI 3.1415926

int main (int argc, char *argv[])
{
   double h;
   double l;
   int i, its;
   double x;
   double d;
   int active;
   unsigned short xi[3];
   int reflected, absorbed, transmitted;

   xi[0] = 0;
   xi[1] = 1;
   xi[2] = 2;

   its = 10000000; /* Book says 10 million tests */
   h = 1.0;
   while (h < 10.1) {
      reflected = 0;
      absorbed = 0;
      transmitted = 0;
      for (i = 1; i <= its; i++) {
         d = 0.0;
         x = 0.0;
         active = 1;
         while (active) {
            x += -(1.0/C)*log(erand48(xi))* cos(d);
            if (x < 0.0) {
               reflected++;
               active = 0;
            } else if (x >= h) {
               transmitted++;
               active = 0;
            } else if (erand48(xi) <  (CC/C)) {
               absorbed++;
               active = 0;
            }
            d = erand48(xi) * PI;
         }
      }
      printf ("H = %d, Reflected     Absorbed      Transmitted\n",
         (int) (h+0.1));
      printf ("       %14.2f%14.2f%14.2f\n\n",
         (float)(100.0*reflected)/(float)its,
         (float)(100.0*absorbed)/(float)its,
         (float)(100.0*transmitted)/(float)its);
      h += 1.0;
   }
}
