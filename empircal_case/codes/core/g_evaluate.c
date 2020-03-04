/*
 * g_evaluate.c:
 *  This file contains the functions to implement g-dominance check.
 *
 * Authors:
 *  Minhui Liao <>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization for Learning and Adaptive System (COLA) Laboratory @ University of Exeter
 *
 * Copyright (c) 2018 Minhui Liao, Ke Li
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

# include "../header/population.h"
# include "../header/dominance.h"

void g_evaluate_population (population_real *mixed_pop)
{
   int i;
   int flag;

   for (i = 0; i < 2 * popsize; i++)
   {
       flag = check_g_dominance (&(mixed_pop->ind[i]));
       if (flag == 0)
       calc_gvalue (&(mixed_pop->ind[i])) ;
    }

    return;
}

void calc_gvalue (individual_real *ind)
{
    int i;

    for (i = 0; i < number_objective; i++)
        ind->obj[i] = ind->obj[i] + 5;

    return;
}
