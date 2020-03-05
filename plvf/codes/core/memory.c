/*
 * memory.c:
 *  This file contains the functions to perform memory operations.
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

#include "../header/global.h"
#include "../header/rand.h"
#include "../header/population.h"

/* Function to allocate memory to an individual */
void allocate_memory_ind (individual_real *ind)
{
    if (number_variable != 0)
    {
        ind->xreal = (double *)malloc(number_variable * sizeof(double));
    }
    ind->obj = (double *)malloc(number_objective * sizeof(double));
    return;
}

/* Function to allocate memory to a population */
void allocate_memory_pop (population_real *pop, int size)
{
    int i;

    pop->ind = (individual_real *)malloc(size * sizeof(individual_real));
    for (i = 0; i < size; i++)
    {
        allocate_memory_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to deallocate memory to an individual */
void deallocate_memory_ind (individual_real *ind)
{
    if (number_variable != 0)
    {
        free(ind->xreal);
    }

    free(ind->obj);

    return;
}

/* Function to deallocate memory to a population */
void deallocate_memory_pop (population_real *pop, int size)
{
    int i;

    for (i = 0; i < size; i++)
    {
        deallocate_memory_ind (&(pop->ind[i]));
    }

    free (pop->ind);

    return;
}
