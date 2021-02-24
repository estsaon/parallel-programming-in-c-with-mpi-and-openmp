/*
 *   Solution to Exercise 3-7
 */

int main (int argc, char *argv[])
{
   int i;
   int id;  /* Task ID */
   int p;   /* Number of tasks */
   int pow; /* ID diff between sending/receiving procs */

   if (argc != 3) {
      printf ("Command line syntax: %s <p> <id>\n", argv[0]);
      exit (-1);
   }
   p = atoi(argv[1]);
   id = atoi(argv[2]);
   if ((id < 0) || (id >= p)) {
      printf ("Task id must be in range 0..p-1\n");
      exit (-1);
   }
   if (p == 1) {
      printf ("No messages sent or received when only one task\n");
      exit (-1);
   }
   pow = 1;
   while(2*pow < p) pow *= 2;
   while (pow > 0) {
      if ((id < pow) && (id + pow <= p))
         printf ("Message received from task %d\n", id + pow);
      else if (id >= pow) {
         printf ("Message sent to task %d\n", id - pow);
         break;
      }
      pow /= 2;
   }
}
