/*
 * crossover_nsga2.c:
 *  This file implements the crossover operations of NSGA-II (SBX).
 *
 * Authors:
 *  Qi Xu <qixu.student@gmail.com>
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

void crossover_spea2 (population_real *parent_pop, population_real *offspring_pop)
{
    int i;
    int rand1, rand2;
    individual_real *parent1, *parent2;

    for (i = 0; i < popsize; i += 2)
    {
        rand1 = rnd (0, popsize - 1);
        rand2 = rnd (0, popsize - 1);
        parent1 = tournament_min (&parent_pop->ind[rand1], &parent_pop->ind[rand2]);
        rand1 = rnd (0, popsize - 1);
        rand2 = rnd (0, popsize - 1);
        parent2 = tournament_min (&parent_pop->ind[rand1], &parent_pop->ind[rand2]);
        sbx_crossover (parent1, parent2, &offspring_pop->ind[i], &offspring_pop->ind[i + 1]);
    }
    return;
}
