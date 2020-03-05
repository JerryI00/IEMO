/*
 * crossover_nsga2.c:
 *  This file implements the crossover operations of NSGA-II (SBX).
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

# include "../../header/reproduction.h"

void crossover_real (population_real *parent_pop, population_real *offspring_pop)
{
    int i, temp, rand;
    int *a1, *a2;
    individual_real *parent1, *parent2;

    a1 = (int *) malloc (popsize * sizeof(int));
    a2 = (int *) malloc (popsize * sizeof(int));
    for (i = 0; i < popsize; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < popsize; i++)
    {
        rand     = rnd (i, popsize - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        rand     = rnd (i, popsize - 1);
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }
    for (i = 0; i < popsize; i += 4)
    {
        parent1 = tournament (&parent_pop->ind[a1[i]], &parent_pop->ind[a1[i + 1]]);
        parent2 = tournament (&parent_pop->ind[a1[i + 2]], &parent_pop->ind[a1[i + 3]]);
        sbx_crossover (parent1, parent2, &offspring_pop->ind[i], &offspring_pop->ind[i + 1]);
        parent1 = tournament (&parent_pop->ind[a2[i]], &parent_pop->ind[a2[i + 1]]);
        parent2 = tournament (&parent_pop->ind[a2[i + 2]], &parent_pop->ind[a2[i + 3]]);
        sbx_crossover (parent1, parent2, &offspring_pop->ind[i + 2], &offspring_pop->ind[i + 3]);
    }
    free (a1);
    free (a2);

    return;
}
