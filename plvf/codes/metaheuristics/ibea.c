/*
 * nsga2.c:
 *  This file contains the main procedures of the standard IBEA.
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

void IBEA (population_real* parent_pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i;
    int generation;
    int *flag;
    double **Indicator_value;

    nadir_point = (double *) malloc (number_objective * sizeof(double));
    ideal_point = (double *) malloc (number_objective * sizeof(double));
    for(i=0;i < number_objective;i++)
        nadir_point[i]=-INF;
    for(i=0;i < number_objective;i++)
        ideal_point[i]=INF;
    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);

    flag    = (int *) malloc (2 * popsize * sizeof(int));
    Indicator_value = (double **) malloc (2 * popsize * sizeof(double*));
    for (i = 0; i < 2 * popsize; i++)
        Indicator_value[i] = (double *) malloc (2 * popsize * sizeof(double*));

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while (evaluation_count < max_evaluation)
    {
        generation++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_real (parent_pop, offspring_pop);
        mutation_real (offspring_pop);
        evaluate_population (offspring_pop);

        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);
        for(i=0;i < 2*popsize;i++)
            update_ideal_point(&(mixed_pop->ind[i]));
        for(i=0;i < 2*popsize;i++)
            update_nadir_point(&(mixed_pop->ind[i]));
        ibea_selection(mixed_pop, parent_pop, flag, Indicator_value);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    // release memory
    for (i = 0; i < 2 * popsize; i++)
        free (Indicator_value[i]);
    free (Indicator_value);
    free (flag);

    return;
}
