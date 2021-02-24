/*
 *   C program solving the traffic circle problem posed in Section 10.5.6.
 *
 *   Programmed by Michael J. Quinn
 */

#include <stdlib.h>
#include <stdio.h>

#define CIRCLE_SIZE 16
#define ENTRANCES   4
#define N_ENTRANCE  0
#define W_ENTRANCE  4
#define S_ENTRANCE  8
#define E_ENTRANCE  12

#define EMPTY       -1

int offset[4] = { 0, 4, 8, 12 };

/* Mean time between car arrivals at each entrance: N W S E */
float f[4] = { 3.0, 3.0, 4.0, 2.0 };

/* Probability that car entering at i will leave at j: N W S E */
float d[4][4] = { 
   0.1, 0.2, 0.5, 0.2,
   0.2, 0.1, 0.3, 0.4,
   0.5, 0.1, 0.1, 0.3,
   0.3, 0.4, 0.2, 0.1,
};

unsigned short xi[3];
int enter[ENTRANCES];
int wait[ENTRANCES];
int entrance_count[ENTRANCES];
int entrance[ENTRANCES];

int main (int argc, char *argv[]) {

   int circle[CIRCLE_SIZE];
   int clicks;
   int i, j;

   void initialize_circle (int *, int *);
   void cars_arrive (int *);
   void cars_enter_circle (int *, int *);
   void move_cars_in_circle (int *);
   void print_circle (int *);
   
   clicks = atoi(argv[1]);
   xi[0] = atoi(argv[2]);
   xi[1] = atoi(argv[3]);
   xi[2] = atoi(argv[4]);

   initialize_circle (circle, entrance);

   for (j = 0; j < ENTRANCES; j++) {
      enter[j] = 1;
      wait[j] = 0;
      entrance_count[j] = 0;
   }
   
   for (i = 0; i < clicks; i++) {
      cars_arrive (entrance);
      move_cars_in_circle (circle);
      cars_enter_circle (circle, entrance);
      for (j = 0; j < ENTRANCES; j++) entrance_count[j] += entrance[j];
      if (i % 1000 == 0) {
         printf ("Time = %d\n", i+1);
         print_circle (circle);
         printf ("ITERATION %5d)   N     W     S     E\n", i+1);
         printf ("  Cars entered: %4d %5d %5d %5d\n",
            enter[0], enter[1], enter[2], enter[3]);
         printf ("  Wait prob:   %5.2f %5.2f %5.2f %5.2f\n",
            (100.0*wait[0])/enter[0], (100.0*wait[1])/enter[1],
            (100.0*wait[2])/enter[2], (100.0*wait[3])/enter[3]);
         printf ("  Ave queue:   %5.2f %5.2f %5.2f %5.2f\n\n",
            (float) entrance_count[0] / (i+1),
            (float) entrance_count[1] / (i+1),
            (float) entrance_count[2] / (i+1),
            (float) entrance_count[3] / (i+1));
      }
   }
}

void initialize_circle (int *circle, int *entrance) {
   int i;
   for (i = 0; i < CIRCLE_SIZE; i++)
      circle[i] = EMPTY;
   for (i = 0; i < ENTRANCES; i++)
      entrance[i] = 0;
}

void cars_arrive (int *entrance) {
   int i;
   
   for (i = 0; i < ENTRANCES; i++) {
      /* Is a car arriving at entrance i? */
      if (erand48(xi) <= 1.0/f[i]) {
         enter[i]++;
         if (entrance[i] > 0) wait[i]++;
         entrance[i]++;
      }
   }
}

void move_cars_in_circle (int *circle) {
   int i, j;
   void print_circle (int *);

   int new_circle[CIRCLE_SIZE];
   for (i = 0; i < CIRCLE_SIZE; i++) {
      j = (i + 1) % CIRCLE_SIZE;
      if ((circle[i] == EMPTY) || (circle[i] == j)) new_circle[j] = EMPTY;
      else new_circle[j] = circle[i];
   }
   for (i = 0; i < CIRCLE_SIZE; i++)
      circle[i] = new_circle[i];
}

void print_circle (int *circle) {
   int i;
   printf ("%2d %7d %7d %7d\n", enter[0], enter[1],
      enter[2], enter[3]);
   printf ("%2d %7d %7d %7d\n", wait[0], wait[1], wait[2], wait[3]);
   printf ("%2d %7d %7d %7d\n", entrance_count[0], entrance_count[1],
      entrance_count[2], entrance_count[3]);
   printf ("%2d %7d %7d %7d\n", entrance[0], entrance[1], entrance[2],
      entrance[3]);
   printf (" N       W       S       E\n");
   for (i = 0; i < CIRCLE_SIZE; i++)
      switch (circle[i]) {
         case EMPTY: printf (" -");
                     break;
         case N_ENTRANCE: printf (" N");
                          break;
         case W_ENTRANCE: printf (" W");
                          break;
         case S_ENTRANCE: printf (" S");
                          break;
         case E_ENTRANCE: printf (" E");
                          break;
      }
   printf ("\n\n");
}

void cars_enter_circle (int *circle, int *entrance) {
   int i, j;
   double dest;
   float cum;

   for (i = 0; i < ENTRANCES; i++) {
      if ((entrance[i] > 0) && (circle[offset[i]] == EMPTY)) {
         entrance[i]--;
         dest = erand48(xi);
         cum = 0.0;
         for (j = 0; j < ENTRANCES; j++) {
            cum += d[i][j];
            if (dest <= cum) {
               switch (j) {
                  case 0: circle[offset[i]] = N_ENTRANCE;
                          break;
                  case 1: circle[offset[i]] = W_ENTRANCE;
                          break;
                  case 2: circle[offset[i]] = S_ENTRANCE;
                          break;
                  case 3: circle[offset[i]] = E_ENTRANCE;
                          break;
               }
               break;
            }
            if (circle[offset[i]] != EMPTY) break;
         }
      }
   }
}
