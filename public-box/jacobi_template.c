
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>              /* For timing */
#include <sys/time.h>            /* For timing */
#include <sys/resource.h>
#include <string.h>
#include <pthread.h>
#include <stdbool.h>
#include <semaphore.h>


/****************Global****************************/

#define MAX(a,b) ((a)>(b)?(a):(b))
#define EPSILON 0.001            /* Termination condition */

char filename[PATH_MAX];                  /* File name of output file */

/* Grid size */
int M = 200;                     /* Number of rows */
int N = 200;                     /* Number of cols */
long max_its = 1000000;          /* Maximum iterations, a safe bound to avoid infinite loop */
double final_diff;               /* Temperature difference between iterations at the end */

/* Thread count */
int thr_count = 2;

/* shared variables between threads */
/*************************************************************/
double** u;                   /* Previous temperatures */
double** w;                   /* New temperatures */


// (1) Add your variables here

/* Information assigned to each worker thread */
typedef struct {
    int id;                 // Worker thread ID
    int start_row;          // First row assigned to this worker
    int end_row;            // Row after the last assigned row
    double diff;             // Maximum temperature difference
    struct rusage usage;    // Execution statistics
} WorkerArgs;

/* Thread management */
pthread_t *threads;
WorkerArgs *worker_args;

/* Synchronization variables */
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_cond_t start_cond = PTHREAD_COND_INITIALIZER;
pthread_cond_t done_cond = PTHREAD_COND_INITIALIZER;

/* Shared state */
int generation = 0;       // Current iteration number
int completed = 0;        // Number of workers finished
int terminate = 0;        // Whether workers should terminate




/**************************************************************/

int main (int argc, char *argv[])
{
   int      its;                 /* Iterations to converge */
   double   elapsed;             /* Execution time */
   struct timeval stime, etime;  /* Start and end times */
   struct rusage usage;

   void allocate_2d_array (int, int, double ***);
   void initialize_array (double ***);
   void print_solution (char *, double **);
   int  find_steady_state (void);

   /* For convenience of other problem size testing */
   if ((argc == 1) || (argc == 4)) {
      if (argc == 4) {
         M = atoi(argv[1]);
         N = atoi(argv[2]);
         thr_count = atoi(argv[3]);
      } // Otherwise use default grid and thread size
   } else {
     printf("Usage: %s [ <rows> <cols> <threads ]>\n", argv[0]);
     exit(-1);
   }

   printf("Problem size: M=%d, N=%d\nThread count: T=%d\n", M, N, thr_count);

   /* Create the output file */
   snprintf(filename, sizeof(filename), "%s_%d_%d_%d.dat", argv[0], M, N, thr_count);

   allocate_2d_array (M, N, &u);
   allocate_2d_array (M, N, &w);
   initialize_array (&u);
   initialize_array (&w);

   gettimeofday (&stime, NULL);
   its = find_steady_state();
   gettimeofday (&etime, NULL);

   elapsed = ((etime.tv_sec*1000000+etime.tv_usec)-(stime.tv_sec*1000000+stime.tv_usec))/1000000.0;

   printf("%s after %d iterations with error: %8.6f.\n",
      final_diff <= EPSILON ? "Converged" : "Iteration limit reached", its, final_diff);
   printf("Elapsed time = %8.4f sec.\n", elapsed);

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

/* Entry function of the worker threads */
void *thr_func(void *arg) {

// (2) Add the worker's logic here
/* Entry function of the worker threads */

    WorkerArgs *args = (WorkerArgs *)arg;

    int i, j;
    int my_generation = 0;

    double local_diff;

    while (1) {

        /* =====================================
           STEP 1: Wait for the master to start
           ===================================== */

        pthread_mutex_lock(&mutex);

        while (generation == my_generation && !terminate) {
            pthread_cond_wait(&start_cond, &mutex);
        }

        /* Check whether we should terminate */

        if (terminate) {
            pthread_mutex_unlock(&mutex);
            break;
        }

        /* Remember which iteration we are doing */

        my_generation = generation;

        pthread_mutex_unlock(&mutex);


        /* =====================================
           STEP 2: Calculate assigned rows
           ===================================== */

        local_diff = 0.0;

        for (i = args->start_row; i < args->end_row; i++) {

            for (j = 1; j < N - 1; j++) {

                w[i][j] = 0.25 * (
                              u[i-1][j] +
                              u[i+1][j] +
                              u[i][j-1] +
                              u[i][j+1]
                          );

                double change = fabs(w[i][j] - u[i][j]);

                if (change > local_diff) {
                    local_diff = change;
                }
            }
        }

        /* Save this worker's maximum difference */

        args->diff = local_diff;


        /* =====================================
           STEP 3: Report completion
           ===================================== */

        pthread_mutex_lock(&mutex);

        completed++;

        /* If everyone has finished, notify master */

        if (completed == thr_count) {
            pthread_cond_signal(&done_cond);
        }

        pthread_mutex_unlock(&mutex);

    }


    /* =====================================
       STEP 4: Collect thread execution stats
       ===================================== */

    getrusage(RUSAGE_THREAD, &args->usage);

    return NULL;


}


int find_steady_state (void)
{

// (3) Implement the thread creation and the main control logic here

    int its;
    int i;
    double diff;
    double **temp;

    struct rusage master_usage;

    /* =====================================
       STEP 1: Allocate thread information
       ===================================== */

    threads = malloc(thr_count * sizeof(pthread_t));

    worker_args = malloc(thr_count * sizeof(WorkerArgs));

    if (threads == NULL || worker_args == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(EXIT_FAILURE);
    }


    /* =====================================
       STEP 2: Create worker threads
       ===================================== */

    for (i = 0; i < thr_count; i++) {

        worker_args[i].id = i;

        /* Divide the interior rows between workers */

        worker_args[i].start_row =
            1 + i * (M - 2) / thr_count;

        worker_args[i].end_row =
            1 + (i + 1) * (M - 2) / thr_count;

        worker_args[i].diff = 0.0;

        int rc = pthread_create(
            &threads[i],
            NULL,
            thr_func,
            &worker_args[i]
        );

        if (rc != 0) {
            fprintf(stderr, "pthread_create failed: %s\n",
                    strerror(rc));
            exit(EXIT_FAILURE);
        }
    }


    /* =====================================
       STEP 3: Perform Jacobi iterations
       ===================================== */

    for (its = 1; its <= max_its; its++) {

        /* Wake up all workers */

        pthread_mutex_lock(&mutex);

        completed = 0;

        generation++;

        pthread_cond_broadcast(&start_cond);

        pthread_mutex_unlock(&mutex);


        /* Wait until all workers finish */

        pthread_mutex_lock(&mutex);

        while (completed < thr_count) {
            pthread_cond_wait(&done_cond, &mutex);
        }

        pthread_mutex_unlock(&mutex);


        /* =====================================
           STEP 4: Find global maximum difference
           ===================================== */

        diff = 0.0;

        for (i = 0; i < thr_count; i++) {

            if (worker_args[i].diff > diff) {
                diff = worker_args[i].diff;
            }
        }


        /* =====================================
           STEP 5: Swap the two matrices
           ===================================== */

        temp = u;
        u = w;
        w = temp;


        /* =====================================
           STEP 6: Check termination condition
           ===================================== */

        if (diff <= EPSILON) {
            break;
        }
    }


    /* Store the final temperature difference */

    final_diff = diff;


    /* =====================================
       STEP 7: Tell all workers to terminate
       ===================================== */

    pthread_mutex_lock(&mutex);

    terminate = 1;

    pthread_cond_broadcast(&start_cond);

    pthread_mutex_unlock(&mutex);


    /* =====================================
       STEP 8: Wait for all workers to exit
       ===================================== */

    for (i = 0; i < thr_count; i++) {

        pthread_join(threads[i], NULL);

        /* Print worker execution statistics */

        printf("Thread %d has completed - user: %.4f s, system: %.4f s\n",
               i,
               worker_args[i].usage.ru_utime.tv_sec +
                   worker_args[i].usage.ru_utime.tv_usec / 1000000.0,
               worker_args[i].usage.ru_stime.tv_sec +
                   worker_args[i].usage.ru_stime.tv_usec / 1000000.0
        );
    }


    /* =====================================
       STEP 9: Get master thread statistics
       ===================================== */

    getrusage(RUSAGE_THREAD, &master_usage);

    printf("find_steady_state - user: %.4f s, system: %.4f s\n",
           master_usage.ru_utime.tv_sec +
               master_usage.ru_utime.tv_usec / 1000000.0,
           master_usage.ru_stime.tv_sec +
               master_usage.ru_stime.tv_usec / 1000000.0
    );


    /* =====================================
       STEP 10: Clean up
       ===================================== */

    free(threads);
    free(worker_args);

    pthread_mutex_destroy(&mutex);

    pthread_cond_destroy(&start_cond);
    pthread_cond_destroy(&done_cond);


    /* Return completed iteration count */

    return its > max_its ? max_its : its;
}