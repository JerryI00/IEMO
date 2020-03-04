/*
 * update.c:
 *  This file contains the update subproblem procedure in MOEA/D.
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

/* Population update procedure of MOEA/D */
void update_subproblem (population_real *pop, individual_real *individual, int subProblem_id, int neighborType)
{
    int i, k;
    int time, size;
    double f1, f2;
    int *perm;

    if (neighborType == NEIGHBOR)
        size = neighbor_size;
    else
        size = number_weight;

    perm = malloc (sizeof(int) * size);
    random_permutation (perm, size);

    time = 0;
    for (i = 0; i < size; i++)
    {
        if (neighborType == NEIGHBOR)
            k = neighborhood[subProblem_id][perm[i]];
        else
            k = perm[i];

        f1 = fitnessFunction (&(pop->ind[k]), lambda[k]);
        f2 = fitnessFunction (individual, lambda[k]);

        if (f2 < f1 )
        {
            copy_ind (individual,&(pop->ind[k]));
            time++;
        }
        if (time >= maximumNumberOfReplacedSolutions)
            break;
    }

    free (perm);

    return;
}

/* Population update procedure of MOEA/D with constraint handling */
void update_subproblem_constraint (population_real *pop, individual_real *individual, int subProblem_id, int neighborType)
{
    int i, k;
    int time, size;
    double f1, f2, cv1, cv2;
    int *perm;

    if (neighborType == NEIGHBOR)
        size = neighbor_size;
    else
        size = number_weight;

    perm = malloc (sizeof(int) * size);
    random_permutation (perm, size);

    time = 0;
    for (i = 0; i < size; i++)
    {
        if (neighborType == NEIGHBOR)
            k = neighborhood[subProblem_id][perm[i]];
        else
            k = perm[i];

        f1  = fitnessFunction (&(pop->ind[k]), lambda[k]);
        f2  = fitnessFunction (individual, lambda[k]);
        cv1 = pop->ind[k].cv;
        cv2 = individual->cv;
        if ((cv1 > -EPS && cv2 > -EPS && f2 < f1 ) || (cv2 > cv1))
        {
            copy_ind (individual, &(pop->ind[k]));
            time++;
        }
        if (time >= maximumNumberOfReplacedSolutions)
            break;
    }

    free (perm);

    return;
}
