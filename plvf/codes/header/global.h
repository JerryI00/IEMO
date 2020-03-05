/*
 * global.h:
 *  This is the header file for global variables and data structures.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

# ifndef Samaritan_GLOBAL_H
# define Samaritan_GLOBAL_H

# include <pthread.h>
# include <stdio.h>
# include <stdlib.h>
# include <stdarg.h>
# include <time.h>
# include <math.h>
# include <float.h>
# include <string.h>
# include <unistd.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <limits.h>
# include "vector.h"

# define PI  M_PI
# define INF 1.0e14
# define EPS 1.0e-7

# define BUFSIZE_S 64
# define BUFSIZE_M 128
# define BUFSIZE_L 256

/* common paramters */
extern int run_index;
extern int run_index_begin;
extern int run_index_end;
extern int max_evaluation;              // maximum number of evaluations (stopping criterion)
extern int evaluation_count;            // evaluation counter
extern int popsize;                     // population size
extern int number_variable;             // number of variables
extern int number_objective;            // number of objectives
extern double *ideal_point;             // ideal point
extern double *nadir_point;             // nadir point
extern double *variable_lowerbound;     // variable lower bound
extern double *variable_upperbound;     // variable upper bound
extern char dummy[BUFSIZE_S];
extern char problem_name[BUFSIZE_S];
extern char algorithm_name[BUFSIZE_S];
extern char analyse_stream[BUFSIZE_L];
extern char problem_param_stream[BUFSIZE_L];
/* crossover and mutation */
extern double eta_c;                    // eta_c in SBX
extern double eta_m;                    // eta_m in polynomial mutation
extern double pcross_real;              // crossover rate for real encoded
extern double pmut_real;                // mutation rate for real encoded
extern double CR;                       // CR in DE
extern double F;                        // F in DE
extern double K;

/* performance metrics */
extern int PF_size;                 // size of the true Pareto-optimal Front
extern double **PF_data;            // true Pareto-optimal front data
extern double *ref_point;           // reference point for Hypervolume calculation

/* MOEA/D variants */
extern int neighbor_size;                           // neighborhood length
extern int number_weight;                           // number of weight vectors
extern char* weight_file;
extern int function_type;                           // type of the aggregation function
extern int maximumNumberOfReplacedSolutions;        // the maximum replacement number of a superior offspring
extern double neighborhood_selection_probability;   // probability to replace in the neighborhood
extern double **lambda;                             // weight vectors
extern double **rbf_lambda;
extern int **neighborhood;                          // neighborhood structure
extern int *permutation;                            // subproblem index permutation
extern int *frequency;                              // subproblem usages counter arrary
extern double* utility;                             // subproblem utility array
extern struct int_vector *selected;
extern struct int_vector *candidate;

/* analysis platform */
extern int runtime_output;
extern int output_interval;
extern int analyse_list[BUFSIZE_S];
extern FILE *pythonplot;
extern pthread_t *plot_thread;

enum analyse_name{VAR, FUN, GD, IGD, HV, PLOT, END};
enum NeighborType{NEIGHBOR, POPULATION};
enum MoeadFunction{WS, N_WS, TCH, N_TCH, ITCH, N_ITCH, PBI, N_PBI};

// gaps setting, read from file?
static int weight_gaps_table[8][3] = {{0,   0, 0},
                                      {0,   0, 0},
                                      {299, 0, 0},  // 2 obj 100->100
                                      {23,  0, 0},  // 3 obj 12->91
                                      {10,  0, 0},  // 4 obj
                                      {6,   4, 0},  // 5 obj
                                      {5,   2, 0},  // 6 obj
                                      {4,   3, 0}}; // 7 obj

typedef struct
{
    int rank;
    double *xreal;
    double *obj;
    double fitness;
    double cv;
} individual_real;

typedef struct
{
    individual_real *ind;
} population_real;

typedef struct lists
{
    int index;
    int index2;
    struct lists *parent;
    struct lists *child;
} list;

typedef struct double_lists
{
    double value;
    struct double_lists *parent;
    struct double_lists *child;
} double_list;

void insert (list *node, int x);
list* del (list *node);
void append (list* node, list* to_list);
int length(list* node);
list* get_item(list* start, int n);
list* free_list (list* node);

# endif // Samaritan_GLOBAL_H
