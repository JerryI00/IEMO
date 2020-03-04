/*
 * utility.c:
 *  This file contains the functions to facilitate some common usages for MOEA/D-STM variants.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
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

# include "../../header/utility.h"
# include "../../header/selection.h"

int prefers (int x, int y, struct double_with_index *womanPref, int size)
{
    int i;
    int pref;

    for (i = 0; i < size; i++)
    {
        pref = womanPref[i].idx;
        if (pref == x)
            return 1;
        if (pref == y)
            return 0;
    }

    return 0;
}

double calculateDistance2 (individual_real *solution, double *lambda)
{
    int i;
    double sum, distance;
    double *vecInd;
    double *normalized_obj;

    vecInd         = malloc (sizeof(double) * number_objective);
    normalized_obj = malloc (sizeof(double) * number_objective);

    sum = 0.0;
    for (i = 0; i < number_objective; i++)
        sum += solution->obj[i];
    for (i = 0; i < number_objective; i++)
        normalized_obj[i] = solution->obj[i] / sum;
    distance = euclidian_distance (normalized_obj, lambda, number_objective);

    free (vecInd);
    free (normalized_obj);

    return distance;
}

/* Find the stable matching between subproblems and solutions */
void stableMatching (int *statusMan, int *statusWoman, int *next, struct double_with_index **man_pref, struct double_with_index **woman_pref, int menSize, int womenSize)
{
    int i;
    int m1;
    int NOT_ENGAGED; // Indicates the mating status
    struct int_vector *freeMen;

    NOT_ENGAGED = -1;
    for (i = 0; i < womenSize; i++)
    {
        next[i]        = 0;
        statusWoman[i] = NOT_ENGAGED;
    }

    freeMen        = malloc (sizeof(struct int_vector));
    freeMen->value = INT_MIN;
    freeMen->next  = NULL;

    for (i = 0; i < menSize; i++)
        int_vector_pushback (freeMen, i);

    while (freeMen->next!=NULL)
    {
        int m = int_vector_pop (freeMen);
        int w = man_pref[m][next[m]].idx;
        next[m]++;
        if (statusWoman[w] == NOT_ENGAGED)
        {
            statusMan[m]   = w;
            statusWoman[w] = m;
        }
        else
        {
            m1 = statusWoman[w];
            if (prefers (m, m1, woman_pref[w], menSize))
            {
                statusMan[m]   = w;
                statusWoman[w] = m;
                int_vector_pushback (freeMen, m1);
            }
            else
            {
                int_vector_pushback (freeMen, m);
            }
        }
    }

    return;
}
