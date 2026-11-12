/*
 *   Sequential program that solves the steady-state temperature
 *   distribution problem using the Jacobi method.
 *
 *   gcc -O3 -o jacobi_seq jacobi_seq.c -lm
 *
 *   Execute:
 *   ./jacobi_seq [optional: <rows> <cols>]
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>              /* For timing */
#include <sys/time.h>            /* For timing */
#include <sys/resource.h>
#include <string.h>

/****************Global****************************/

#define MAX(a,b) ((a)>(b)?(a):(b))
#define EPSILON 0.001            /* Termination condition */

char filename[PATH_MAX];                  /* File name of output file */

double** u;                   /* Previous temperatures */
double** w;                   /* New temperatures after each iteration */

/* Grid size */
int M = 200;                     /* Number of rows */
int N = 200;                     /* Number of cols */
long max_its = 1000000;          /* Maximum iterations, a safe bound to avoid infinite loop */
double final_diff;               /* Temperature difference between iterations at the end */

/**************************************************/

int main (int argc, char *argv[])
{

   int      its;                 /* Iterations to converge */
   double   elapsed;             /* Execution time */
   struct timeval stime, etime;  /* Start and end times */
   struct rusage usage;

   void allocate_2d_array (int, int, double ***);
   void initialize_array (double ***);
   void print_solution (char *, double **);
   int  find_steady_state ();

   /* For convenience of other problem size testing */
   if ((argc == 1) || (argc == 3)) {
      if (argc == 3) {
         M = atoi(argv[1]);
         N = atoi(argv[2]);
      } // Otherwise use default grid size
   } else {
     printf("Usage: %s [ <rows> <cols> ]\n", argv[0]);
     exit(-1);
   }

   printf("Problem size: M=%d, N=%d\n", M, N);
   snprintf(filename, sizeof(filename), "%s_%d_%d.dat", argv[0], M, N);

   allocate_2d_array (M, N, &u);
   allocate_2d_array (M, N, &w);
   initialize_array (&u);
   initialize_array (&w);

   gettimeofday (&stime, NULL);
   its = find_steady_state ();
   gettimeofday (&etime, NULL);

   elapsed = ((etime.tv_sec*1000000+etime.tv_usec)-(stime.tv_sec*1000000+stime.tv_usec))/1000000.0;

   printf("%s after %d iterations with error: %.6f.\n",
      final_diff <= EPSILON ? "Converged" : "Iteration limit reached", its, final_diff);
   printf("Elapsed time = %.4f sec.\n", elapsed);

   /* get the resource usage of the whole process */
   getrusage(RUSAGE_SELF, &usage);
   printf("Program completed - user: %.4f s, system: %.4f s\n",
      (usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.0),
    (usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.0));
   printf("no. of context switches: vol %ld, invol %ld\n\n", 
 		  usage.ru_nvcsw, usage.ru_nivcsw);

   print_solution (filename, u);
}

/* Allocate two-dimensional array. */
void allocate_2d_array (int r, int c, double ***a)
{
   double *storage;
   int     i;
   storage = (double *) malloc (r * c * sizeof(double));
   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &storage[i * c];
}

/* Set initial and boundary conditions */
void initialize_array (double ***u)
{
   int i, j;

   /* Set initial values and boundary conditions */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         (*u)[i][j] = 25.0;      /* Room temperature */
      (*u)[i][0] = 0.0;
      (*u)[i][N-1] = 0.0;
   }

   for (j = 0; j < N; j++) {
      (*u)[0][j] = 0.0;
      (*u)[M-1][j] = 1000.0;     /* Heat source */
   }
}

/* Print solution to standard output or a file */
void print_solution (char *filename, double **u)
{
   int i, j;
   char sep;
   FILE *outfile;

   if (!filename) { /* if no filename specified, print on screen */
      sep = '\t';   /* tab added for easier view */
      outfile = stdout;
   } else {
      sep = '\n';   /* for gnuplot format */
      outfile = fopen(filename,"w");
      if (outfile == NULL) {
         printf("Can't open output file.");
         exit(-1);
      }
   }

   /* Print the solution array */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         fprintf (outfile, "%6.2f%c", u[i][j], sep);
      fprintf(outfile, "\n"); /* Empty line for gnuplot */
   }
   if (outfile != stdout)
      fclose(outfile);

}

int find_steady_state ()
{
   double diff;            /* Maximum temperature difference */
   int    its;             /* Iteration count */
   int    i, j;
   double **temp;

   for (its = 1; its <= max_its; its++) {
      diff = 0.0;
      for (i = 1; i < M-1; i++) {
         for (j = 1; j < N-1; j++) {
            w[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
            if (fabs(w[i][j] - u[i][j]) > diff)
               diff = fabs(w[i][j] - u[i][j]);
         }
      }

      /* Swap matrix u, w by exchanging the pointers */
      temp = u;
      u = w;
      w = temp;

      /* Terminate if temperatures have converged */
      if (diff <= EPSILON)
         break;
   }

   final_diff = diff;  //return the diff value via global variable
   return its > max_its ? max_its : its;  //return completed iterations
}
