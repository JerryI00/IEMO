/*
 * R-MEAD2.c:
 *  This file contains the main procedures of R-MEAD2. It is based on the following reference
 *
 *  A. Mohammadi, M. N. Omidvar, X. Li, and K. Deb. Integrating user preferences and decomposition methods for
 *  many-objective optimiza-tion. In Evolutionary Computation (CEC), 2014 IEEE Congress on, pages 421–428. IEEE, 2014.
 *
 * Authors:
 *  Joe Billingsley <jb931@exeter.ac.uk>
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
#include "../header/global.h"

void RMEAD2 (population_real *parent_pop, population_real *offspring_pop, population_real *mixed_pop, double* reference_point, double radius)
{
    int generation;

    generation       = 1;
    evaluation_count = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", run_index);

    // initialization process
    lambda = (double **) malloc (popsize * sizeof(double *));

    for (int i = 0; i < popsize; i++)
        lambda[i] = (double *) malloc(number_objective * sizeof(double));

    uniform_random_weights_const(lambda, popsize, number_objective, 0, 1);

    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);
    initialize_idealpoint (parent_pop);

    track_evolution (parent_pop, generation, 0);

    while (evaluation_count < max_evaluation)
    {
        print_progress ();

        // crossover and mutation
        crossover_real (parent_pop, offspring_pop);
        mutation_real (offspring_pop);
        evaluate_population (offspring_pop);

        for(int i = 0; i < popsize; i++)
            update_ideal_point (&offspring_pop->ind[i]);

        for(int i = 0; i < popsize; i++) {
            double offspring_fitness = fitnessFunction(&offspring_pop->ind[i], lambda[i]);
            double parent_fitness = fitnessFunction(&parent_pop->ind[i], lambda[i]);

            if(offspring_fitness < parent_fitness)
                copy_ind(&offspring_pop->ind[i], &parent_pop->ind[i]);
        }

        int nearest_weight_idx = 0;
        double min_dist = INF;
        for(int i = 0; i < popsize; i++) {
            double dist = euclidian_distance(parent_pop->ind[i].obj, reference_point, number_objective);

            if(dist < min_dist) {
                nearest_weight_idx = i;
                min_dist = dist;
            }
        }

        double lower_bounds[number_objective], upper_bounds[number_objective];

        for(int i = 0; i < number_objective; i++) {
            lower_bounds[i] = lambda[nearest_weight_idx][i] - radius;
            upper_bounds[i] = lambda[nearest_weight_idx][i] + radius;
        }

        uniform_random_weights(lambda, popsize, number_objective, lower_bounds, upper_bounds);

        generation++;

        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    free (ideal_point);

    for (int i = 0; i < number_weight; i++)
        free (lambda[i]);
    free (lambda);

    ideal_point = NULL;
    lambda = NULL;
}
