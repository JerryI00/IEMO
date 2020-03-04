/*
 * stm_selection.c:
 *  This file contains the functions for stable matching based selection procedure.
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

# include "../../header/selection.h"
# include "../../header/metaheuristics.h"

/* Environmental selection based on stable matching */
void stm_selection (population_real *parent_pop, population_real *mixed_pop)
{
    int i, j;
    int min_index;

    // Calculate the preference values of solution matrix
    for (i = 0; i < 2 * number_weight; i++)
    {
        min_index = 0;
        for (j = 0; j < number_weight; j++) {
            fitnessMatrix[i][j] = fitnessFunction (&(mixed_pop->ind[i]), lambda[j]);
            distMatrix[i][j]  	= calculateDistance2 (&(mixed_pop->ind[i]), lambda[j]);
            if (distMatrix[i][j] < distMatrix[i][min_index])
                min_index = j;
        }
        nicheCount[min_index] = nicheCount[min_index] + 1;
    }

    // calculate the preference values of subproblem matrix and solution matrix
    for (i = 0; i < 2 * number_weight; i++)
    {
        for (j = 0; j < number_weight; j++) {
            subpMatrix[j][i].x   = fitnessFunction (&(mixed_pop->ind[i]), lambda[j]);
            subpMatrix[j][i].idx = i;
            solMatrix[i][j].x    = distMatrix[i][j] + nicheCount[j];
            solMatrix[i][j].idx  = j;
        }
    }

    // sort the preferences of subproblems and solutions
    for (i = 0; i < number_weight ; i++)
        qsort (subpMatrix[i], 2 * number_weight, sizeof(struct double_with_index), double_with_index_greater_cmp);
    for (i = 0; i < 2 * number_weight ; i++)
        qsort (solMatrix[i], number_weight, sizeof(struct double_with_index), double_with_index_greater_cmp);

    // find the stable matching between subproblems and solutions
    stableMatching (idx, statusWoman, next, subpMatrix, solMatrix, number_weight, number_weight * 2);

    for (i = 0; i < number_weight; i++)
        copy_ind (&(mixed_pop->ind[idx[i]]), &(parent_pop->ind[i]));

    return;
}

/* Environmental selection based on stable matching (dynamic resource allocation version) */
void stm_dra_selection (population_real *parent_pop, population_real *mixed_pop, int size)
{
    int i, j;
    int min_index;

    // Calculate the preference values of solution matrix
    for (i = 0; i < size; i++)
    {
        min_index = 0;
        for (j = 0; j < number_weight; j++) {
            fitnessMatrix[i][j] = fitnessFunction (&(mixed_pop->ind[i]), lambda[j]);
            distMatrix[i][j]  	= calculateDistance2 (&(mixed_pop->ind[i]), lambda[j]);
            if (distMatrix[i][j] < distMatrix[i][min_index])
                min_index = j;
        }
        nicheCount[min_index] = nicheCount[min_index] + 1;
    }

    // calculate the preference values of subproblem matrix and solution matrix
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < number_weight; j++) {
            subpMatrix[j][i].x   = fitnessFunction (&(mixed_pop->ind[i]), lambda[j]);
            subpMatrix[j][i].idx = i;
            solMatrix[i][j].x    = distMatrix[i][j] + nicheCount[j];
            solMatrix[i][j].idx  = j;
        }
    }

    // sort the preferences of subproblems and solutions
    for (i = 0; i < number_weight ; i++)
        qsort (subpMatrix[i], size, sizeof(struct double_with_index), double_with_index_greater_cmp);
    for (i = 0; i < size ; i++)
        qsort (solMatrix[i], number_weight, sizeof(struct double_with_index), double_with_index_greater_cmp);

    // find the stable matching between subproblems and solutions
    stableMatching(idx,statusWoman,next,subpMatrix, solMatrix, number_weight, size);

    for (i = 0; i < number_weight; i++)
        copy_ind (&(mixed_pop->ind[idx[i]]), &(parent_pop->ind[i]));

    return;
}
