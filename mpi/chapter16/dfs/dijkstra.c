/*   Parallel depth-first search that incorporates Dijkstra's
 *   distributed termination detection algorithm.
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 9 November 2002
 */
   
#include <stdio.h>
#include <mpi.h>

#define CUTOFF_DEPTH 4
#define MAX_DEPTH    8

#define SOLUTION_TAG 1
#define STATUS_CHECK_TAG 2
#define TERMINATION_TAG 3
#define WHITE 0
#define BLACK 1

struct status_check {
   int count;
   int color;
} ;

int moves[MAX_DEPTH];
int p;
int id;
int cutoff_count;
int color;
int msg_count;

struct position {
   long i;
};

void terminate_gracefully (int, int);


int main (int argc, char *argv[]) {
   struct position board;

   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);

   get_initial_position (&board);
   msg_count = 0;
   color = WHITE;
   par_dfs (&board, 0);
   terminate_gracefully(p, id);
}

par_dfs (struct position *board, int level)
{
   int i;
   int possible_moves;
   int flag;
   MPI_Status status;

   if (level == MAX_DEPTH) {
      if (solution(board)) {
         if (id > 0) {
            MPI_Send (moves, MAX_DEPTH, MPI_INT, 0,
               SOLUTION_TAG, MPI_COMM_WORLD);
            msg_count++;
            color = BLACK;
         }
         terminate_gracefully (p, id);
      }
   } else if ((level == CUTOFF_DEPTH) &&
            (cutoff_count++ % p != id)) return;
   else {
      if (level == CUTOFF_DEPTH) {
         MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
            &flag, &status);
         if (flag == 1) {
            terminate_gracefully(p, id);
         }
      }
      possible_moves = count_moves (board);
      for (i = 0; i < possible_moves; i++) {
         make_move(board, i);
         moves[level] = i;
         par_dfs (board, level+1);
         unmake_move(board, i);
      }
   }
}

count_moves (struct position *board)
{
return 2;
}

make_move (struct position *board, int move)
{
}

unmake_move (struct position *board, int i)
{
}

get_initial_position (struct position *board)
{
}

solution (struct position *board)
{
   if (cutoff_count == 3) return 1;
   else return 0;
}

void terminate_gracefully (int p, int id)
{
   MPI_Status status;
   int flag;
   int started_status_check;
   struct status_check scr;
   int i;

   if (!id) {
      started_status_check = 0;
      for (;;) {
         for (;;) {
            MPI_Iprobe (MPI_ANY_SOURCE, SOLUTION_TAG, MPI_COMM_WORLD,
               &flag, &status);
            if (!flag) break;
            MPI_Recv (moves, MAX_DEPTH, MPI_INT, MPI_ANY_SOURCE, 
               SOLUTION_TAG, MPI_COMM_WORLD, &status);
            color = BLACK;
            msg_count--;
         }
         if (!started_status_check) {
            scr.count = 0;
            scr.color = WHITE;
            MPI_Send (&scr, sizeof(struct status_check), MPI_CHAR,
               (id+1)%p, STATUS_CHECK_TAG, MPI_COMM_WORLD);
            color = WHITE;
            started_status_check = 1;
         }
         MPI_Iprobe (MPI_ANY_SOURCE, STATUS_CHECK_TAG, MPI_COMM_WORLD,
            &flag, &status);
         if (flag) {
            MPI_Recv (&scr, sizeof(struct status_check), MPI_CHAR,
               MPI_ANY_SOURCE, STATUS_CHECK_TAG, MPI_COMM_WORLD, &status);
            if ((color == WHITE) && (scr.color == WHITE) &&
             ((scr.count + msg_count) == 0)) {
               for (i = 1; i < p; i++)
                  MPI_Send (NULL, 0, MPI_CHAR, i, TERMINATION_TAG,
                     MPI_COMM_WORLD);
               printf ("Process 0 prints a solution:\n");
               for (i = 0; i < MAX_DEPTH; i++) printf ("%d-", moves[i]);
               printf ("\n");
               fflush (stdout);
               MPI_Finalize();
               exit(0);
            } else {
               scr.count = 0;
               scr.color = WHITE;
               MPI_Send (&scr, sizeof(struct status_check), MPI_CHAR,
               (id + 1)%p, STATUS_CHECK_TAG, MPI_COMM_WORLD);
               color = WHITE;
            }
         }
      }
   } else {
      for (;;) {
         MPI_Iprobe (0, TERMINATION_TAG, MPI_COMM_WORLD, &flag, &status);
         if (flag) {
            MPI_Finalize();
            exit(0);
         }
         MPI_Iprobe (MPI_ANY_SOURCE, STATUS_CHECK_TAG, MPI_COMM_WORLD,
            &flag, &status);
         if (flag) {
            MPI_Recv (&scr, sizeof(struct status_check), MPI_CHAR,
               MPI_ANY_SOURCE, STATUS_CHECK_TAG, MPI_COMM_WORLD, &status);
            if (color == BLACK) scr.color = BLACK;
            scr.count += msg_count;
            MPI_Send (&scr, sizeof(struct status_check), MPI_CHAR,
               (id + 1)%p, STATUS_CHECK_TAG, MPI_COMM_WORLD);
            color = WHITE;
         }
      }
   }
}
