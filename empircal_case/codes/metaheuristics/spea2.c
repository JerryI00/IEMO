/*
 * spea2.c:
 *  This file contains the main procedures of the standard SPEA2.
 *
 * Authors:
 *  Qi Xu <qixu.student@gmail.com>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Qi Xu, Renzhi Chen, Ke Li
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

void SPEA2(population_real *parent_pop, population_real *archive, population_real *mixed_pop)
{
    // Here, the offspring is actually used as archive, so the name is changed for convenience.
    int archive_size = popsize;
    int total_size = popsize + archive_size;
    int k_min = (int)sqrt(total_size);
    int i, generation;

    int *dominated_Num;            // Number of an individual dominating others
    int *bedominated_Num;        // Number of an individual dominated by others
    int **dominated_Matrix;
    int *R_i;

    double **distance_Matrix;
    double *D_i;
    double *kth_distance;
    individual_real *temp_ind;
    population_real *temp_pop;

    temp_ind         = (individual_real *) malloc (sizeof(individual_real));
    temp_pop         = (population_real *) malloc (sizeof(population_real));
    dominated_Num    = (int *) malloc(total_size * sizeof(int));
    bedominated_Num  = (int *) malloc(total_size * sizeof(int));
    R_i              = (int *) malloc(total_size * sizeof(int));
    D_i              = (double *) malloc(total_size * sizeof(double));
    kth_distance     = (double *) malloc(total_size * sizeof(double));
    dominated_Matrix = (int **) malloc(total_size * sizeof(int *));
    distance_Matrix  = (double **) malloc(total_size * sizeof(double *));

    for (i = 0; i < total_size; i++)
    {
        dominated_Matrix[i] = (int *) malloc(total_size * sizeof(int));
        distance_Matrix[i]  = (double *) malloc(total_size * sizeof(double));
    }
    allocate_memory_ind (temp_ind);
    allocate_memory_pop (temp_pop, total_size);

    generation = 1;
    evaluation_count = 0;
    printf("Progress: 1%%");

    // initialize population
    initialize_population_real(parent_pop);
    evaluate_population(parent_pop);
    for (i = 0; i < popsize; i++)
    {
        copy_ind(&(parent_pop->ind[i]), &(archive->ind[i]));
    }

    // track the current evolutionary progress, including population and metrics
    track_evolution(archive, generation, 0);
    while (evaluation_count < max_evaluation)
    {
        generation++;
        print_progress();

        // reproduction
        crossover_spea2(archive, parent_pop); // Mating selection included here.
        mutation_real(parent_pop);
        evaluate_population(parent_pop);

        // environmental selection
        merge(parent_pop, archive, mixed_pop);
        fitness_spea2(mixed_pop,total_size,k_min,dominated_Num,bedominated_Num,dominated_Matrix,R_i,distance_Matrix,D_i,kth_distance);
        selection_spea2(mixed_pop,total_size,archive,archive_size,temp_ind,temp_pop,distance_Matrix);

        // track the current evolutionary progress, including population and metrics
        track_evolution(archive, generation, evaluation_count >= max_evaluation);
    }

    free(dominated_Num);
    free(bedominated_Num);
    free(R_i);
    free(D_i);
    free(kth_distance);
    for (i = 0; i < total_size; i++)
    {
        free(dominated_Matrix[i]);
        free(distance_Matrix[i]);
    }
    free(distance_Matrix);
    free(dominated_Matrix);
    deallocate_memory_ind (temp_ind);
    deallocate_memory_pop (temp_pop, total_size);
    free (temp_ind);
    free (temp_pop);

    return;
}















