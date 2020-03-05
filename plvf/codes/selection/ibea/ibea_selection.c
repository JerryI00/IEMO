/*
 * ibea_selection.c:
 *  This file contains the environmental selection function for IBEA.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization for Learning and Adaptive System (COLA) Laboratory @ University of Exeter
 *
 * Copyright (c) 2018 Renzhi Chen, Ke Li
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

# include "../../header/selection.h"
#include "../../header/global.h"

void environmental_selection (population_real *mixed_pop, population_real *new_pop, int *flag, double **Indicator_value, int size)
{
    int i, j;
    int worst, new_size;



    for (i = 0; i < size; i++)
        flag[i] = 0;

    for (i = size - popsize; i > 0; i--)
    {
        for (j = 0; j < size && flag[j] == 1; j++);

        worst = j;

        for (j = j + 1; j < size; j++)
        {
            if (flag[j] != 1)
            {
                if (mixed_pop->ind[j].fitness <
                    mixed_pop->ind[worst].fitness)
                    worst = j;
            }
        }

        for (j = 0; j < size; j++)
            if (flag[j] != 1)
                mixed_pop->ind[j].fitness -= Indicator_value[worst][j];

        flag[worst] = 1;
    }

    //Move remaining solutions.
    new_size = 0;
    for (i = 0; i < size; i++)
    {
        if (flag[i] != 1)
        {
            copy_ind (&mixed_pop->ind[i], &new_pop->ind[new_size]);
            new_size++;
        }
    }

    return;
}

void ibea_selection (population_real *mixed_pop, population_real *new_pop, int *flag, double **Indicator_value)
{
    int size;
    size = 2 * popsize;
    calcFitness (mixed_pop, Indicator_value, size);
    environmental_selection (mixed_pop, new_pop, flag, Indicator_value, size);
}

void pbea_selection (population_real *mixed_pop, population_real *new_pop, int *flag, double **Indicator_value, double* weights, double* reference_point, double specificity)
{
    int size;

    size = 2 * popsize;
    pbea_calcFitness (mixed_pop, Indicator_value, size, weights, reference_point, specificity);
    environmental_selection (mixed_pop, new_pop, flag, Indicator_value, size);
}
