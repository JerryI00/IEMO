/*
 * mutation_real.c:
 *  This file contains the functions to perform mutation operations for real-coded individuals.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Ke Li, Renzhi Chen
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

#include "../../header/reproduction.h"

void mutation_real (population_real *pop)
{
    int i;

    for (i = 0; i < popsize; i++)
        mutation_ind (&(pop->ind[i]));

    return;
}


/* Function to perform mutation of an individual */
void mutation_ind (individual_real *ind)
{
    if (number_variable != 0)
        polymut_ind (ind);

    return;
}
