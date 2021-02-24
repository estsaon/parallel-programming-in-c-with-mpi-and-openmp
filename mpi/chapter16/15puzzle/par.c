/*   Parallel best-first branch-and-bound algorithm to solve
 *   Loyd's 15-puzzle.
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last revision: 6 November 2002
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define POSITIONS     16
#define HOLE          (POSITIONS-1)
#define DIMENSION     4
#define MAX_MOVES     60
#define MAX_HEAP_SIZE 5000000
#define MAX_NUM_PROCS 16
#define COMM_INTERVAL 0.01

#define WHITE 1
#define BLACK -1

#define MAX(a,b)  ((a)>(b)?(a):(b))
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define PARENT(i) (((i)-1)/2)
#define LCHILD(i) (2*(i)+1)
#define RCHILD(i) (2*((i)+1))
#define SUCCESSOR(i) (((i)+1)%p)

struct puzzle {
   int moves_made;
   int lower_bound;
   int hole;
   char val[POSITIONS];
   char move[MAX_MOVES];
};
typedef struct puzzle puzzle;

int possible_moves[POSITIONS] = { 2, 3, 3, 2,
                                  3, 4, 4, 3,
                                  3, 4, 4, 3,
                                  2, 3, 3, 2 };

int new_hole[POSITIONS][4] = {
   1, 4, -1, -1,
   0, 2,  5, -1,
   1, 3,  6, -1,
   2, 7, -1, -1,
   0, 5,  8, -1,
   1, 4,  6,  9,
   2, 5,  7, 10,
   3, 6, 11, -1,
   4, 9, 12, -1,
   5, 8, 10, 13,
   6, 9, 11, 14,
   7, 10, 15, -1,
   8, 13, -1, -1,
   9, 12, 14, -1,
   10, 13, 15, -1,
   11, 14, -1, -1 
};
  

puzzle heap[MAX_HEAP_SIZE];
int heap_size;

void initialize_heap()
{
   heap_size = 0;
}

int superior (int i, int j)
{
   if (heap[i].lower_bound < heap[j].lower_bound) return 1;
   if (heap[j].lower_bound < heap[i].lower_bound) return 0;
   if (heap[i].moves_made <= heap[j].moves_made) return 1;
   return 0;
}

void insert_heap(puzzle b)
{
   puzzle tmp;
   int    index;

   if (heap_size == MAX_HEAP_SIZE) {
      printf ("Heap overflow\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
      exit(-1);
   }
   heap[heap_size] = b;
   index = heap_size;
   heap_size++;
   while ((index > 0) && superior(index, PARENT(index))) {
      tmp = heap[index];
      heap[index] = heap[PARENT(index)];
      heap[PARENT(index)] = tmp;
      index = PARENT(index);
   }
}


puzzle delete_heap (void)
{
   puzzle return_val;
   puzzle tmp;
   int favored_child;
   int index;

   if (heap_size == 0) {
      printf ("Heap underflow\n");
      exit (-1);
   }
   return_val = heap[0];
   heap[0] = heap[heap_size-1];
   heap_size--;
   index = 0;
   while (LCHILD(index) < heap_size) {
      if (RCHILD(index) >= heap_size) favored_child = LCHILD(index);
      else if (superior(RCHILD(index),LCHILD(index)))
         favored_child = RCHILD(index);
      else favored_child = LCHILD(index);
      if (superior(index, favored_child)) break;
      tmp = heap[index];
      heap[index] = heap[favored_child];
      heap[favored_child] = tmp;
      index = favored_child;
   }
   return return_val;
}

int solved (puzzle b)
{
   int i;

   for (i = 0; i < POSITIONS; i++)
      if (b.val[i] != i) return 0;
   return 1;
}

int lower_bound (puzzle b)
{
   int i;
   int sum;
   int diff;
   
   sum = 0;
   for (i = 0; i < POSITIONS; i++) {
      if (b.val[i] != HOLE) {
         diff = MAX(b.val[i]-i, i-b.val[i]);
         sum += (diff / DIMENSION) + (diff % DIMENSION);
      }
   }
   return sum;
}

void get_puzzle (puzzle *b)
{
   int    i;
   int    d;

   for (i = 0; i < POSITIONS; i++) {
      scanf ("%d", &d);
      if (d > 0) b->val[i] = (char) (d-1);
      else {
         b->val[i] = HOLE;
         b->hole = i;
      }
   }
   b->lower_bound = lower_bound(*b);
   b->moves_made = 0;
}

puzzle make_move (puzzle b, int i)
{
   b.val[b.hole] = b.val[new_hole[b.hole][i]];
   b.move[b.moves_made] = b.val[b.hole]+1;
   b.val[new_hole[b.hole][i]] = HOLE;
   b.hole = new_hole[b.hole][i];
   b.moves_made++;
   b.lower_bound = b.moves_made + lower_bound(b);
   return b;
}

print_puzzle (puzzle b)
{
   int i, j;
   
   for (i = 0; i < POSITIONS; i++)
   printf ("%d ", (int) b.val[i]);
   printf ("\n\n");
   for (i = 0; i < DIMENSION; i++) {
      for (j = 0; j < DIMENSION; j++)
         if (b.val[DIMENSION*i+j] != HOLE)
            printf ("%3d", (int) (b.val[DIMENSION*i+j])+1);
         else
            printf (" --");
      if (i == 0) printf ("   Hole: %d\n", b.hole);
      else if (i == 1) printf ("   Lower bound: %d\n", b.lower_bound);
      else if (i == 2) printf ("   Moves made: %d\n", b.moves_made);
      else {
         printf ("   ");
         for (j = 0; j < b.moves_made; j++)
            printf ("%d-", b.move[j]);
         printf ("\n");
      }
   }
   printf ("\n");
}

print_solution (puzzle b)
{
   int i;
   printf ("Here is the solved puzzle:\n");
   print_puzzle (b);
   printf ("Here are the moves needed to solve the puzzle:\n");
   for (i = 0; i < b.moves_made; i++)
      printf ("%d-", b.move[i]);
   printf ("\n");
}

unsigned short xi[3];

int partner (int p, int id)
{
   int rv;
   do {
      rv = (int) (p * erand48(xi));
   } while (rv == id);
   return rv;
}

#define UNEXAMINED_SUBPROBLEM_TAG 1
#define STATUS_CHECK_TAG          2
#define TERMINATION_TAG           5

struct status_check {
   puzzle global_s;
   int    global_c;
   int    color;
   int    count;
} ;

int msg_count;
int color;
int best_soln;

void BandB_Communication (int p, int id, int local_c, puzzle local_s,
 int *latest_best)
{
   int        dummy1[MAX_NUM_PROCS], dummy2[MAX_NUM_PROCS];
   int        flag;
   int        i;
   int        min;
   MPI_Status status;
   puzzle     u, v;
   struct     status_check scr;
   int        dest;
   int        dummy;

   MPI_Iprobe (MPI_ANY_SOURCE, TERMINATION_TAG, MPI_COMM_WORLD, &flag,
      &status);
   if (flag) {
         MPI_Finalize ();
         exit (0);
   }
   MPI_Iprobe (MPI_ANY_SOURCE, STATUS_CHECK_TAG, MPI_COMM_WORLD, &flag,
      &status);
   if (flag) {
      MPI_Recv (&scr, sizeof(struct status_check), MPI_CHAR, MPI_ANY_SOURCE,
         STATUS_CHECK_TAG, MPI_COMM_WORLD, &status);
      if (local_c < scr.global_c) {
         scr.global_c = local_c;
         scr.global_s = local_s;
      }
      if (scr.global_c <= lower_bound(heap[0])) heap_size = 0;
      best_soln = scr.global_c;
      if (!id) {
         if ((color == WHITE) && (scr.color == WHITE) &&
             ((scr.count + msg_count) == 0)) {
            for (i = 1; i < p; i++)
               MPI_Send (NULL, 0, MPI_CHAR, i, TERMINATION_TAG,
                  MPI_COMM_WORLD);
            printf ("Solution reported by process %d:\n", id);
            print_puzzle (scr.global_s);
            MPI_Finalize();
            exit(0);
         }
         color = WHITE;
         scr.color = WHITE;
         scr.count = 0;
         MPI_Send (&scr, sizeof(struct status_check), MPI_CHAR,
                  SUCCESSOR(id), STATUS_CHECK_TAG, MPI_COMM_WORLD);
      } else {
         if ((color == BLACK) && (scr.color == WHITE)) scr.color = BLACK;
         color = WHITE;
         scr.count += msg_count;
         MPI_Send (&scr, sizeof(struct status_check), MPI_CHAR,
            SUCCESSOR(id), STATUS_CHECK_TAG, MPI_COMM_WORLD);
      }
   }
   for (;;) {
      MPI_Iprobe (MPI_ANY_SOURCE, UNEXAMINED_SUBPROBLEM_TAG,
         MPI_COMM_WORLD, &flag, &status);
      if (!flag) break;
      MPI_Recv (&u, sizeof(puzzle), MPI_CHAR, MPI_ANY_SOURCE,
         UNEXAMINED_SUBPROBLEM_TAG, MPI_COMM_WORLD, &status);
      msg_count--;
      color = BLACK;
      if (lower_bound(u) < best_soln) insert_heap (u);
    }
   if (heap_size > 1) {
      u = delete_heap();
      v = delete_heap();
      dest = partner (p, id);
      MPI_Send (&v, sizeof(puzzle), MPI_CHAR, dest,
         UNEXAMINED_SUBPROBLEM_TAG, MPI_COMM_WORLD);
      msg_count++;
      color = BLACK;
      insert_heap (u);
   }
}

int main (int argc, char *argv[])
{
   puzzle local_s;
   int    local_c;
   int    i;
   puzzle s, t, u;
   int    p;
   int    id;
   struct status_check scr;
   int    latest_best;
   double last_comm;
   double     start_time;
   double     elapsed;
   double     comm_time;
   int        cnt;

   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   if (p < 2) {
      if (!id) printf ("Need at least two processes\n");
      MPI_Finalize();
      exit (0);
   }
   xi[0] = 0;
   xi[1] = 1;
   xi[2] = id;
   initialize_heap();
   if (id == 0) {
      get_puzzle(&s);
      print_puzzle(s);
      insert_heap(s);
      scr.global_c = 999;
      scr.color = WHITE;
      scr.count = 0;
      MPI_Send (&scr, sizeof(struct status_check), MPI_CHAR,
         SUCCESSOR(id), STATUS_CHECK_TAG, MPI_COMM_WORLD);
   }
   local_c = best_soln = 999;
   last_comm = MPI_Wtime();
   msg_count = 0;
   for (;;) {
      if ((heap_size == 0) || ((MPI_Wtime()-last_comm) > COMM_INTERVAL)) {
         comm_time -= MPI_Wtime();
         BandB_Communication(p, id, local_c, local_s, &latest_best);
         last_comm = MPI_Wtime();
/*
         comm_time += MPI_Wtime();
         elapsed = MPI_Wtime()-start_time;
         if (!id) {
            cnt++;
            if (cnt == 100) {
            printf ("Process %d is spending %5.2f percent "
                    "of its time communicating\n",
                    id, 100.0*comm_time/elapsed);
            fflush (stdout);
            cnt = 0;
            }
         }
*/
      }
      if (heap_size > 0) {
         u = delete_heap();
         if (u.lower_bound < best_soln) {
            color = BLACK;
            if (solved(u)) {
               local_s = u;
               local_c = u.lower_bound;
            } else {
               for (i = 0; i < possible_moves[u.hole]; i++) {
                  t = make_move (u, i);
                  if (lower_bound(t) < best_soln)
                     insert_heap (t);
               }
            }
         }
      }
   }
}
