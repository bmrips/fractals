#include <complex.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>
#include <sys/param.h>
#include <sys/sysinfo.h>
#include <time.h>

#define PI 3.14159265358979323846264338327950288

// Mutex for item_done list
pthread_mutex_t item_done_mutex;

// Arrays for communication between computing and writing threads
unsigned char ** roots; // The matrix that holds the results of roots
unsigned char ** iters; // The matrix that holds the results of iterations
char * item_done; // To mark a row as done

// The exact roots
double * exact_roots_real;
double * exact_roots_imag;

double step_size; // Step size for starting points of newton method

int n_rows = 1000; // The size of the image (number of rows and columns)
int n_threads; // The number of threads

#define TOLERANCE 1e-6 // The tolerance
#define CANCEL_THRESHOLD 1e10 // The thresold at which to cancel the iteration
int max_iters = 100; // The iterations value cap in the convergence image

int exp_d = 0; // The power d

/**
 * Parse the arguments and save the results in the global variables.
 */
void
parse_args (int argc,
            char * argv[])
{// {{{
    char * rem; // The remainder of the conversion from string to integer

    // Iterate over the input arguments
    for (size_t i = 1; i < argc; i++) {

        // Option -h: help
        if (strncmp ("-h", argv[i], 2) == 0) {
            printf (
                "Usage: newton [OPTION]... <d>\n"
                "Prints the newton fractal (attractors and iterations) of the complex function\n"
                "f(z) = z^d-1.\n\n"
                "Various options may be given:\n"
                "  -m<n>  Cap the number of iterations in the convergence image at n. Must be\n"
                "         less than 256. Defaults to 100 iterations.\n"
                "  -h     Print this help.\n"
                "  -s<n>  Make an image with a resolution of n×n. Defaults to 1000×1000.\n"
                "  -t<n>  Use <n> threads. By default all available cores are used.\n");
            exit (EX_OK);
        }

        // Option -m: The iterations value cap in the convergence image
        else if (strncmp ("-m", argv[i], 2) == 0) {
            max_iters = (int) strtol (argv[i]+2, &rem, 0);
            if (max_iters == 0 || rem[0] != '\0') {
                printf ("Error: option -m requires a number as argument! Found: %s\n", argv[i]+2);
                exit (EX_USAGE);
            } else if (max_iters <= 0 || max_iters >= 256) {
                printf ("Error: iterations cap passed to the -m flag has to be in the interval [1,255]!\n");
                exit (EX_USAGE);
            }
        }

        // Option -s: size of the image (number of rows and columns)
        else if (strncmp ("-s", argv[i], 2) == 0) {
            n_rows = (int) strtol (argv[i]+2, &rem, 0);
            if (n_rows == 0 || rem[0] != '\0') {
                printf ("Error: option -s requires a number as argument! Found: %s\n", argv[i]+2);
                exit (EX_USAGE);
            } else if (n_rows <= 0) {
                printf ("Error: number of rows passed to the -s flag has to be positive!\n");
                exit (EX_USAGE);
            }
        }

        // Option -t: number of threads
        else if (strncmp ("-t", argv[i], 2) == 0) {
            n_threads = (int) strtol (argv[i]+2, &rem, 0);
            if (n_threads == 0 || rem[0] != '\0') {
                printf ("Error: option -t requires a number as argument! Found: %s\n", argv[i]+2);
                exit (EX_USAGE);
            } else if (n_threads <= 0) {
                printf ("Error: number of threads passed to the -t flag has to be positive!\n");
                exit (EX_USAGE);
            }
        }

        // Unknown option
        else if (strncmp("-", argv[i], 1) == 0) {
            printf ("Error: unkown option: %s\n", argv[i]);
            exit (EX_USAGE);
        }

        // Positional argument: exponent d
        else if (exp_d == 0) {
            exp_d = (int) strtol (argv[i], &rem, 0);
            if (exp_d == 0 || rem[0] != '\0') {
                printf ("Error: positional argument has to be a number! Found: %s\n", argv[i]);
                exit (EX_USAGE);
            } else if (exp_d <= 0) {
                printf ("Error: exponent has to be positive!\n");
                exit (EX_USAGE);
            }
        }

        // Too many arguments
        else {
            printf ("Error: too many arguments!\n");
            exit (EX_USAGE);
        }
    }

    // Check that positional arguments have been passsed
    if (exp_d == 0) {
        printf ("Error: exonent as positional argument missing!\n");
        exit (EX_USAGE);
    }
}// }}}

/**
 * Calculate z^d.
 */
inline double complex
power (double complex z,
       int d)
{// {{{
    double complex res = 1;
    while (d > 0) {
        // If d is odd, multiply the result with z
        if (d & 1)
            res = res * z;

        // d must be even now
        d = d >> 1;
        z = z * z;
    }
    return res;
}// }}}

/**
 * Calculate the next iteration of Newton's method: \Phi(z) = z - f(z)/f'(z)
 */
inline double complex
next_newton (double complex z,
             const int d)
{// {{{
    complex double z_raised = power (z, d-1);
    return (z*z_raised*(d-1) + 1) / (z_raised * d);
}// }}}

/**
 * Apply the Newton method to z
 */
inline void
newton_method (complex double z,
               unsigned char * root,
               unsigned char * iters)
{// {{{
    unsigned int n_iter = 0;
    double z_real, z_imag, z_real_d, z_imag_d;

    // Apply Newton's method until it converges or diverges
    bool cont = true;
    while (cont) {
        z_imag = cimag (z);
        z_real = creal (z);

        // Test on divergence
        if ( (z_real*z_real + z_imag*z_imag) < TOLERANCE
                ||  z_real > CANCEL_THRESHOLD
                || -z_real > CANCEL_THRESHOLD
                ||  z_imag > CANCEL_THRESHOLD
                || -z_imag > CANCEL_THRESHOLD) {
            *root = exp_d;
            *iters = 0;
            return;
        }

        // Test on convergence
        for (size_t mi = 0; mi < exp_d; mi++) {
            z_real_d = z_real - exact_roots_real[mi];
            z_imag_d = z_imag - exact_roots_imag[mi];
            if (((z_real_d * z_real_d) + (z_imag_d * z_imag_d)) < TOLERANCE) {
                *root = mi;
                cont = false;
            }
        }

        n_iter++;
        z = next_newton (z, exp_d); // Iterate z
    }

    // Set number of iterations, but maximal max_iters
    *iters = MIN (n_iter, max_iters);
}// }}}

/**
 * Compute Newton's method for the specifed values and save the results in the ROOTS and ITERS
 * matrices.
 */
void *
compute_main (void * args)
{// {{{
    size_t offset = *((size_t*)args);
    free (args);

    double complex z;
    size_t buffer_size = n_rows * sizeof(unsigned char);

    for (size_t ix = offset; ix < n_rows; ix += n_threads) {

        // Initialize a local buffer
        unsigned char * roots_row = malloc (buffer_size);
        unsigned char * iters_row = malloc (buffer_size);

        // Compute the root using the Newton method for each item in the row ix
        for (size_t jx = 0; jx < n_rows; jx++) {
            z = (-2 + jx*step_size) + I * (2 - ix*step_size);
            newton_method (z, roots_row + jx, iters_row + jx);
        }

        // Write the results to the matrices
        roots[ix] = roots_row;
        iters[ix] = iters_row;

        // Set row ix as done
        pthread_mutex_lock (&item_done_mutex);
        item_done[ix] = 1;
        pthread_mutex_unlock (&item_done_mutex);
    }

    return NULL;
}// }}}

/**
 * Writes the matrices to files.
 */
void *
write_main (void * args)
{// {{{
    // Initialize local item_done copy
    char * item_done_loc = calloc (n_rows, sizeof(char));

    // Define sleeping time in order to prevent deadlock
    struct timespec sleep_time;
    sleep_time.tv_sec = 0;
    sleep_time.tv_nsec = 50000;

    // Initialize files for saving roots and iterations
    char buf_attractors  [35];
    char buf_convergence [35];

    sprintf (buf_attractors,  "newton_attractors_x%02d.pgm",  exp_d);
    sprintf (buf_convergence, "newton_convergence_x%02d.pgm", exp_d);

    FILE * roots_file = fopen (buf_attractors, "wb");
    FILE * iters_file = fopen (buf_convergence, "wb");

    fprintf (roots_file, "P5\n%d %d\n%d\n", n_rows, n_rows, exp_d);
    fprintf (iters_file, "P5\n%d %d\n%d\n", n_rows, n_rows, max_iters);

    for (size_t ix = 0; ix < n_rows;) {

        // Copy item_done if the next item is done, otherwise sleep a short while and check again
        pthread_mutex_lock (&item_done_mutex);
        if (item_done[ix] != 0) {
            memcpy (item_done_loc, item_done, n_rows*sizeof(char));
            pthread_mutex_unlock (&item_done_mutex);
        } else {
            pthread_mutex_unlock (&item_done_mutex);
            nanosleep (&sleep_time, NULL);
            continue;
        }

        // Write to the files
        for (; ix < n_rows && item_done_loc[ix] != 0; ++ix) {
            fwrite (roots[ix], sizeof(unsigned char), n_rows, roots_file);
            fwrite (iters[ix], sizeof(unsigned char), n_rows, iters_file);
            free (roots[ix]);
            free (iters[ix]);
        }
    }

    // Close files and free memory
    fclose (roots_file);
    fclose (iters_file);
    free (item_done_loc);

    // Return a null pointer
    return NULL;
}// }}}

int
main (int argc,
      char * argv[])
{// {{{
    n_threads = get_nprocs() - 1; // Default the number of threads to #CPUs-1

    parse_args (argc, argv);

    step_size = 4.0 / (n_rows-1); // The step size when changing z.

    // Compute the exact roots
    exact_roots_real = malloc (exp_d * sizeof(double));
    exact_roots_imag = malloc (exp_d * sizeof(double));
    double roots_var;
    for (size_t ix = 0; ix < exp_d; ++ix) {
        roots_var = (2*PI*ix) / exp_d;
        exact_roots_real[ix] = cos (roots_var);
        exact_roots_imag[ix] = sin (roots_var);
    }

    // Allocating memory for the roots and iters matrices, and the item_done array
    roots = malloc (n_rows * sizeof(unsigned char *));
    iters = malloc (n_rows * sizeof(unsigned char *));
    item_done = calloc (n_rows, sizeof(char));

    // Threading
    pthread_mutex_init (&item_done_mutex, NULL);
    pthread_t * compute_threads = malloc (n_threads * sizeof(pthread_t));
    pthread_t write_thread;

    // Start threads
    for (size_t tx = 0; tx < n_threads; ++tx) {
        size_t * arg = malloc (sizeof(size_t));
        *arg = tx;
        if (pthread_create (compute_threads + tx, NULL, compute_main, (void*)arg)) {
            printf ("Error: cannot create the computation threads!\n");
            exit (EX_OSERR);
        }
    }
    if (pthread_create (&write_thread, NULL, write_main, NULL)) {
        printf ("Error: cannot create the writing thread!\n");
        exit (EX_OSERR);
    }

    // Close threads
    for (int tx = 0; tx < n_threads; ++tx) {
        if (pthread_join (compute_threads[tx], NULL)) {
            printf ("Error: cannot join the computation threads!\n");
            exit (EX_OSERR);
        }
    }
    if (pthread_join(write_thread, NULL)) {
        printf ("Error: cannot join the writing thread!\n");
        exit (EX_OSERR);
    }

    pthread_mutex_destroy (&item_done_mutex);
    free (compute_threads);

    free (roots);
    free (iters);
    free (item_done);
    free (exact_roots_real);
    free (exact_roots_imag);

    return EX_OK;
}// }}}
