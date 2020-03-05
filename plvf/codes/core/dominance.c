/*
 * dominance.c:
 *  This file contains the functions to perform the dominance check.
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

# include "../header/global.h"
# include "../header/dominance.h"

double *weights_obj;
double *reference_point;

/* Check the Pareto dominance relationship between two solutions */
int check_dominance (individual_real *a, individual_real *b)
{
    int i;
    int flag1;
    int flag2;

    flag1 = flag2 = 0;
    for (i = 0; i < number_objective; i++)
    {
        if (a->obj[i] < b->obj[i])
            flag1 = 1;
        else
        {
            if (a->obj[i] > b->obj[i])
                flag2 = 1;
        }
    }
    if (flag1 == 1 && flag2 == 0)
        return (1);
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return (-1);
        else
            return (0);
    }
}

/* Check the g-dominance value */
int check_g_dominance (individual_real *a)
{
    int i;
    int flag1, flag2;
    
    flag1 = flag2 = 0;
    for (i = 0; i < number_objective; i++)
    {
        if (a->obj[i] < reference_point[i])
            flag1 = 1;
        else
        {
            if (a->obj[i] > reference_point[i])
                flag2 = 1;
        }
    }
    if ((flag1 == 1 && flag2 == 0) || (flag1 == 0 && flag2 == 1))
        return 1;
    else
        return 0;
    
}

/* Check the r-dominance relationship between two solutions */
int check_r_dominance (individual_real *a, individual_real *b, population_real *mixed_pop, double sigma)
{
    int flag1, flag2;

    flag1 = check_dominance (a, b);
    if (flag1 == 1)
        return 1;
    else
    {
        if(flag1 == 0)
        {
            flag2 = dist_measure (mixed_pop,a,b,sigma);
            return flag2;
        }
        else
            return -1;
    }
}

int dist_measure (population_real *mixed_pop, individual_real *a, individual_real *b, double sigma)
{
    int i, j;
    int flag;
    double temp_dist, dist_max, dist_min, distance1, distance2, distance_measure;
    double *max_obj, *min_obj;
    double *distance_array;

    max_obj = (double *) malloc (number_objective * sizeof(double));
    min_obj = (double *) malloc (number_objective * sizeof(double));
    distance_array = (double *) malloc (2 * popsize * sizeof(double));

    // Find the maximum and minimum values of each dimension.
    for (i = 0; i < number_objective; i++)
    {
        min_obj[i] = INF;
        max_obj[i] = -INF;
        for (j = 0; j < 2 * popsize; j++)
        {
            if (max_obj[i] < (mixed_pop->ind[j].obj[i]))
                max_obj[i] = mixed_pop->ind[j].obj[i];
            if (min_obj[i] > (mixed_pop->ind[j].obj[i]))
                min_obj[i] = mixed_pop->ind[j].obj[i];
        }
    }

    // Find the maximum and minimum distance with respect to the reference point(s)
    dist_min = INF;
    dist_max = -INF;
    for (i = 0; i < 2 * popsize; i++)
    {
        temp_dist = 0.0;
        for (j = 0; j < number_objective; j++)
            temp_dist += weights_obj[j] * pow ((mixed_pop->ind[i].obj[j] - reference_point[j])
                    / (max_obj[j] - min_obj[j]), 2);
        distance_array[i] = sqrt (temp_dist);
        if (distance_array[i] > dist_max)
            dist_max = distance_array[i];
        if (distance_array[i] < dist_min)
            dist_min = distance_array[i];
    }

    distance1 = 0;
    distance2 = 0;
    for (i = 0; i < number_objective; i++)
    {
        distance1 += weights_obj[i] * pow ((a->obj[i] - reference_point[i]) / (max_obj[i] - min_obj[i]), 2);
        distance2 += weights_obj[i] * pow ((b->obj[i] - reference_point[i]) / (max_obj[i] - min_obj[i]), 2);
    }
    distance1 = sqrt (distance1);
    distance2 = sqrt (distance2);

    // Discriminating dominance relationship.
    distance_measure = (distance1 - distance2) / (dist_max - dist_min);

    if (distance_measure < -sigma)
        flag = 1;
    else
    {
        if (distance_measure > sigma)
            flag = -1;
        else
            flag = 0;
    }

    free (max_obj);
    free (min_obj);
    free (distance_array);

    return flag;
}
