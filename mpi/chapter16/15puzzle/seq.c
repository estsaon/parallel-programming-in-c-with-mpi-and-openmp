/*
 *   C program to solve the 15-puzzle
 *
 *   Programmed by Michael J. Quinn
 *
 */

#define POSITIONS     16
#define HOLE          (POSITIONS-1)
#define DIMENSION     4
#define MAX_MOVES     60
#define MAX_HEAP_SIZE 5000000


#define MAX(a,b)  ((a)>(b)?(a):(b))
#define PARENT(i) (((i)-1)/2)
#define LCHILD(i) (2*(i)+1)
#define RCHILD(i) (2*((i)+1))

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

int main (int argc, char *argv[])
{
   puzzle best_soln;
   int    best_cost;
   int    i;
   puzzle s, t, u;
   int    lbsf;

   initialize_heap();
   get_puzzle(&s);
   print_puzzle(s);
   insert_heap(s);
   best_cost = 999;
   lbsf = 0;
   while(heap_size > 0) {
      u = delete_heap();
      if (u.lower_bound > lbsf) {
         lbsf = u.lower_bound;
      }
      if (u.lower_bound >= best_cost) break;
      if (solved(u)) {
         if (u.lower_bound < best_cost) {
            s = u;
            best_cost = u.lower_bound;
         }
      } else {
         for (i = 0; i < possible_moves[u.hole]; i++) {
            t = make_move (u, i);
            insert_heap (t);
         }
      }
   }
   print_solution (s);
}
