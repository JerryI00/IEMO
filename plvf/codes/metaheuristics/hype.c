/*
 * hype.c:
 *  This file implements the main procedures of HypE. It is based on the following reference:
 *
 * J. Bader and E. Zitzler, "HypE: An algorithm for fast hypervolume-based many-objective optimization",
 * Evol. Comput. 19(1): 45-76, 2011.
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

void HypE (population_real *parent_pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i, j;
    int generation;
    int maxdepth, maxStackSize;

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);
    initialize_idealpoint (parent_pop);
    initialize_nadirpoint (parent_pop);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while (evaluation_count < max_evaluation)
    {
        generation++;
        print_progress ();

        print_error (popsize % 2 != 0, 1, "Population size shall be even number in HypE!");
        for (i = 0; i < popsize; i += 2)
        {
            fflush (stdout);
            // reproduction (crossover and mutation)
            crossover_real_steadystate (parent_pop, &(offspring_pop->ind[i]), &(offspring_pop->ind[i + 1]));
            mutation_ind (&(offspring_pop->ind[i]));
            mutation_ind (&(offspring_pop->ind[i + 1]));
            evaluate_individual (&(offspring_pop->ind[i]));
            evaluate_individual (&(offspring_pop->ind[i + 1]));

            update_ideal_point(&(offspring_pop->ind[i]));
            update_ideal_point(&(offspring_pop->ind[i + 1]));
            update_nadir_point (&(offspring_pop->ind[i]));
            update_nadir_point (&(offspring_pop->ind[i + 1]));
        }
        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);
        fill_hype_sort (parent_pop, mixed_pop, 2 * popsize);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    return;
}
