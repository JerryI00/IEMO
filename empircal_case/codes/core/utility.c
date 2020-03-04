/*
 * utility.c:
 *  This file contains the functions to facilitate some common usages.
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

# include "../header/utility.h"
# include "../header/rank_sort.h"

/* Calculate the L2-norm of a vector */
double norm_vector (double *a)
{
    int i;
    double sum;

    sum = 0;
    for (i = 0; i < number_objective; i++)
        sum += a[i] * a[i];

    return sqrt (sum);
}

/* Calculate the Euclidean distance between two points */
double euclidian_distance (double *a, double *b, int dimension)
{
    int i;
    double distance;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
        distance += (a[i] - b[i]) * (a[i] - b[i]);

    return sqrt(distance);
}

double normalised_euclidean_distance (const double *a, const double *b, const double* f_max, const double* f_min, int dimension)
{

    int i;
    double distance;

    distance = 0.0;
    for (i = 0; i < dimension; i++)
        distance += pow((a[i] - b[i]) / (f_max[i] - f_min[i]), 2);

    return sqrt(distance);
}


/* Build up multi-level directories */
void _mkdir (const char *dir)
{
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf (tmp, sizeof(tmp), "%s", dir);
    len = strlen (tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
    {
        if (*p == '/')
        {
            *p = 0;
            mkdir (tmp, S_IRWXU);
            *p = '/';
        }
    }
    mkdir (tmp, S_IRWXU);
}

/* Calculate the combinatorial number for n choose k */
int combination (int n, int k)
{
    int i;

    if (n < k)
        return -1;
    double ans = 1;
    for (i = k + 1; i <= n; i++)
    {
        ans = ans * i;
        ans = ans / (double) (i - k);
    }

    return (int) ans;
}

/* Shuffle the 'perm' array */
void random_permutation (int *perm, int size)
{
    int i, num, start;
    int *index, *flag;

    index = malloc (size * sizeof(int));
    flag  = malloc (size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        index[i] = i;
        flag[i]  = 1;
    }

    num = 0;
    while (num < size)
    {
        start = rnd (0, size - 1);
        while (1)
        {
            if (flag[start])
            {
                perm[num] = index[start];
                flag[start] = 0;
                num++;
                break;
            }
            if (start == (size - 1))
                start = 0;
            else
                start++;
        }
    }

    free (index);
    free (flag);

    return;
}

/* Update the current ideal point */
void update_ideal_point (individual_real *individual)
{
    int i;

    for (i = 0; i < number_objective; i++)
        if (individual->obj[i] < ideal_point[i])
            ideal_point[i] = individual->obj[i];

    return;
}


/* Update the current nadir point */
void update_nadir_point (individual_real *individual)
{
    int i;

    for (i = 0; i < number_objective; i++)
        if (individual->obj[i] > nadir_point[i])
            nadir_point[i] = individual->obj[i];

    return;
}

double weighted_euclidean_distance_ASF(const double *x, const double* y, const double* weights,
                                       const double* max, const double* min, int dimension)
{
    double sum, check_sum, frac;

    sum       = 0;
    check_sum = 0;
    for (int i = 0; i < dimension; i++)
    {
        if (weights[i] > 1 || weights[i] < 0)
            print_error (1, 3, "Weight  ", i, " exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

        check_sum += weights[i];
    }

    if (check_sum > 1)
        print_error (1, 1, "Weights exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

    for (int i = 0; i < dimension; i++)
    {
        frac = (x[i] - y[i]) / (max[i] - min[i]);
        sum += weights[i] * frac * frac;
    }

    return sqrt(sum);
}

double tchebycheff_ASF (const double *x, const double* y, const double* weights, int dimensions)
{
    double check_sum = 0;

    for (int i = 0; i < dimensions; i++)
    {
        if(weights[i] > 1 || weights[i] < 0)
            print_error(1, 3, "Weight  ", i, " exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

        check_sum += weights[i];
    }

    if (check_sum > 1)
        print_error (1, 1, "Weights exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

    double curr_diff = 0, max_diff = -INF, diff_sum = 0;
    for (int i = 0; i < dimensions; i++)
    {
        curr_diff = weights[i] * (x[i] - y[i]);

        if (curr_diff > max_diff)
            max_diff = curr_diff;

        diff_sum += curr_diff;
    }

    return max_diff + (0.00001 * diff_sum);
}



void normalise_vector (double* vector, int length)
{
    double z = 0;

    for (int i = 0; i < length; i++)
        z += vector[i];

    for (int i = 0; i < length; i++)
        vector[i] = vector[i] / z;

    return;
}