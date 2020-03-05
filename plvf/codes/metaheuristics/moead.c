/*
 * moead.c:
 *  This file contains the main procedures of the standard MOEAD.
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

void MOEAD (population_real *pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i, j;
    int generation;
    int subproblem_id, neighbor_type;

    generation       = 1;
    evaluation_count = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", run_index);

    // initialization process
  // initialize_uniform_weight();

   nums_weight();
  //  read_uniform_weight("1.dat");
    print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_neighborhood ();
    initialize_population_real (pop);
    evaluate_population (pop);
    initialize_idealpoint (pop);

    track_evolution (pop, generation, 0);

    permutation = malloc (number_weight * sizeof(int));
    individual_real* offspring = &(offspring_pop->ind[0]);

    while (evaluation_count < max_evaluation)
    {
        print_progress ();

        random_permutation (permutation, number_weight);
        for (i = 0; i < number_weight; i++)
        {
            subproblem_id = permutation[i];

            // crossover and mutation
            crossover_moead_real (pop, offspring, subproblem_id, &neighbor_type);
            mutation_ind (offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point (offspring);

            // update subproblem
            update_subproblem (pop, offspring, subproblem_id, neighbor_type);
        }
        generation++;

        track_evolution (pop, generation, evaluation_count >= max_evaluation);
    }

    free (permutation);
    moead_free ();

    return;
}
