/*
 * main.c:
 *  This is the main procedures of a general EMO algorithm (generational evolution model).
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

# include "header/rand.h"
# include "header/metaheuristics.h"
#include <unistd.h>

/* common paramters */
int run_index;
int run_index_begin;
int run_index_end;
int max_evaluation;              // maximum number of evaluations (stopping criterion)
int evaluation_count;            // evaluation counter
int popsize;                     // population size
int number_variable;             // number of variables
int number_objective;            // number of objectives
double* ideal_point;             // ideal point
double* nadir_point;             // nadir point
double* variable_lowerbound;     // variable lower bound
double* variable_upperbound;     // variable upper bound
char dummy[BUFSIZE_S];
char problem_name[BUFSIZE_S];
char algorithm_name[BUFSIZE_S];
char analyse_stream[BUFSIZE_L];
char problem_param_stream[BUFSIZE_L];
/* crossover and mutation */
double eta_c;                    // eta_c in SBX
double eta_m;                    // eta_m in polynomial mutation
double pcross_real;              // crossover rate for real encoded
double pmut_real;                // mutation rate for real encoded
double CR;                       // CR in DE
double F;                        // F in DE
double K;

/* performance metrics */
int PF_size;                 // size of the true Pareto-optimal Front
double **PF_data;            // true Pareto-optimal front data
double *ref_point;           // reference point for Hypervolume calculation

/* MOEA/D variants */
int neighbor_size;                           // neighborhood length
int number_weight;                           // number of weight vectors
char *weight_file;
int function_type;                           // type of the aggregation function
int maximumNumberOfReplacedSolutions;        // the maximum replacement number of a superior offspring
double neighborhood_selection_probability;   // probability to replace in the neighborhood
double **lambda;                             // weight vectors
int **neighborhood;                          // neighborhood structure
int *permutation;                            // subproblem index permutation
int *frequency;                              // subproblem usages counter arrary
double *utility;                             // subproblem utility array
struct int_vector *selected;
struct int_vector *candidate;

/* analysis platform */
int runtime_output;
int output_interval;
int analyse_list[BUFSIZE_S];
FILE *pythonplot;
pthread_t *plot_thread;
double *reference_point, *weights_obj;

int main(int argc, char *argv[])
{
    int i;
    // initialize parameter settings
    initialization_real (argc,argv);

    population_real *parent_pop;
    population_real *offspring_pop;
    population_real *mixed_pop;
    parent_pop    = (population_real *) malloc (sizeof(population_real));
    offspring_pop = (population_real *) malloc (sizeof(population_real));
    mixed_pop     = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (offspring_pop, popsize);
    allocate_memory_pop (mixed_pop, 2 * popsize);

    randomize ();

    /* parameter settings for preference-based methods */
    reference_point = (double *) malloc (sizeof(double) * number_objective);
    printf ("Please input a reference point:\n"); //reference point.
    for (i = 0; i < number_objective; i++)
        scanf ( "%lf", &(reference_point[i]));
    weights_obj = (double *) malloc (sizeof(double) * number_objective);
    for (i = 0; i < number_objective; i++)
        weights_obj[i] = 1 / (double) number_objective;  // weights_obj for different objectives

    double specificity = 0.001;     // the parameter of PBEA
    double sigma       = 0.2;       // the parameter of r-NSGA-II
    double epsilon     = 0.01;    // the parameter of R-NSGA-II

    // run experiments
    for (run_index = run_index_begin; run_index <= run_index_end; run_index++)
    {
        printf ("-----------------------------\n");
        printf ("|\tThe %d run\t|\t", run_index);
        if (!strcmp (algorithm_name, "NSGA2"))
            NSGA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD"))
            MOEAD (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD_DRA"))
            MOEAD_DRA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD_STM"))
            MOEAD_STM (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD_STM_DRA"))
            MOEAD_STM_DRA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "SMSEMOA"))
            SMSEMOA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "IBEA"))
            IBEA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "HYPE"))
            HypE (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "NSGA3"))
            NSGA3 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "SPEA2"))
            SPEA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "rNSGA2"))
            r2NSGA2 (parent_pop, offspring_pop, mixed_pop, sigma);
        else if(!strcmp (algorithm_name, "gNSGA2"))
            gNSGA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "RNSGA2"))
            RNSGA2 (parent_pop, offspring_pop, mixed_pop, reference_point, weights_obj, epsilon);
        else if(!strcmp (algorithm_name, "PBEA"))
            PBEA (parent_pop, offspring_pop, mixed_pop, reference_point, weights_obj, specificity);
        else
            print_error (1, 2, "UNKNOWN algorithm:", algorithm_name);

        printf ("\n");
    }
    printf ("-----------------------------\n");

    // free memory
    if (number_variable != 0)
    {
        free (variable_lowerbound);
        free (variable_upperbound);
    }
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (offspring_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2 * popsize);
    free (parent_pop);
    free (offspring_pop);
    free (mixed_pop);

    for (i = 0; i < PF_size; i++)
        free (PF_data[i]);
    free (PF_data);
    free(reference_point);
    free(weights_obj);

    return 0;
}
