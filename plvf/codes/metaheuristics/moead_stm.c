/*
 * moead_stm.c:
 *  This file contains the all procedures of the MOEAD-STM.
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

# include "../header/metaheuristics.h"

int *idx ;
int *next;
int *statusWoman;
double *nicheCount;
double **distMatrix;
double **fitnessMatrix;
struct double_with_index **solMatrix;
struct double_with_index **subpMatrix;

void MOEAD_STM (population_real *parent_pop, population_real *offspring_pop, population_real *mixed_pop)
{

    int i;
    int generation;
    int subproblem_id, neighbor_type;

    // initialize uniform weight vectors
    initialize_uniform_weight ();
    print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");

    idx           = malloc (sizeof(int) * number_weight);
    nicheCount    = malloc (sizeof(double)* number_weight);
    distMatrix    = malloc (sizeof(double *) * number_weight * 2);
    fitnessMatrix = malloc (sizeof(double *) * number_weight * 2);
    statusWoman   = malloc (sizeof(int) * number_weight * 2);
    next          = malloc (sizeof(int) * number_weight * 2);
    solMatrix     = malloc (sizeof(struct double_with_index *) * number_weight * 2);
    subpMatrix    = malloc (sizeof(struct double_with_index *) * number_weight);

    for (i = 0; i < 2 * number_weight; i++)
    {
        solMatrix[i]     = malloc (sizeof(struct double_with_index) * number_weight);
        distMatrix[i]    = malloc (sizeof(double) * number_weight);
        fitnessMatrix[i] = malloc (sizeof(double) * number_weight);
    }

    for (i = 0; i < number_weight; i++)
        subpMatrix[i] = malloc (sizeof(struct double_with_index) * number_weight * 2);

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize process
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);
    initialize_neighborhood ();
    initialize_idealpoint (parent_pop);
    initialize_nadirpoint(parent_pop);

    permutation = malloc (number_weight * sizeof(int));

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);

    while (evaluation_count < max_evaluation)
    {
        print_progress (generation);

        // offspring generation
        random_permutation (permutation,number_weight);
        for (i = 0; i < number_weight; i++)
        {
            subproblem_id = permutation[i];
            crossover_moead_real (parent_pop, &(offspring_pop->ind[i]), subproblem_id, &neighbor_type);
            mutation_ind (&(offspring_pop->ind[i]));
            evaluate_individual (&(offspring_pop->ind[i]));
            update_ideal_point (&(offspring_pop->ind[i]));
            update_nadir_point(&(offspring_pop->ind[i]));
        }
        merge (parent_pop, offspring_pop, mixed_pop);

        // environmental selection
        stm_selection (parent_pop, mixed_pop);
        generation++;

        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    // free all malloc memory
    for (i = 0; i < 2 * number_weight; i++)
    {
        free (solMatrix[i]);
        free (distMatrix[i]);
        free (fitnessMatrix[i]);
    }

    for (i = 0; i < number_weight; i++)
        free (subpMatrix[i]);
    free (permutation);
    free (idx) ;
    free (nicheCount);
    free (solMatrix);
    free (distMatrix);
    free (fitnessMatrix);
    free (statusWoman);
    free (next);
    moead_free ();
    return;
}