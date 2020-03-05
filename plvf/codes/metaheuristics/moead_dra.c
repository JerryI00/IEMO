/*
 * moead_dra.c:
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

void MOEAD_DRA (population_real* pop, population_real* offspring_pop, population_real* mixed_pop)
{
    int i, j;
    int selected_size;
    int generation_count;
    int subproblem_id, neighbor_type;
    population_real* saved_pop;
    individual_real *offspring;

    generation_count = 1;
    evaluation_count = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", run_index);

    // initialization process
    initialize_uniform_weight ();
    print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_neighborhood ();
    initialize_population_real (pop);
    evaluate_population (pop);
    initialize_idealpoint (pop);

    track_evolution (pop, generation_count, 0);

    utility   = malloc (number_weight * sizeof(double));
    frequency = malloc (number_weight * sizeof(int));
    saved_pop = malloc (number_weight * sizeof(population_real));
    allocate_memory_pop (saved_pop, number_weight);

    for (i = 0; i < number_weight; i++)
    {
        utility[i]   = 1.0;
        frequency[i] = 0;
    }

    offspring = &(offspring_pop->ind[0]);
    while (evaluation_count < max_evaluation)
    {
        //create empty head for selected and candidate
        selected  = malloc (sizeof(struct int_vector));
        candidate = malloc (sizeof(struct int_vector));
        selected->value  = INT_MIN;
        selected->next   = NULL;
        candidate->value = INT_MIN;
        candidate->next  = NULL;

        print_progress ();
        // select the current most active subproblems to evolve (based on utility)
        tour_selection_subproblem (neighbor_size);
        selected_size = int_vector_size (selected);
        for (i = 0; i < selected_size; i++)
        {
            j = int_vector_get (selected, i + 1);
            frequency[j]++;


            // crossover and mutation
            crossover_moead_real (pop, offspring, j, &neighbor_type);
            mutation_ind (offspring);
            evaluate_individual (offspring);

            // update the ideal point
            update_ideal_point (offspring);

            // update the subproblem
            update_subproblem (pop, offspring, j, neighbor_type);
        }

        generation_count++;
        if (generation_count % 30 == 0)
            comp_utility (pop, saved_pop);

        track_evolution (pop, generation_count, evaluation_count >= max_evaluation);

        int_vector_free (selected);
        int_vector_free (candidate);
    }

    moead_free ();
    free (utility);
    free (frequency);
    deallocate_memory_pop (saved_pop, number_weight);
    free (saved_pop);

    return;
}
