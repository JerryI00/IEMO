/*
 * igd.c:
 *  This file contains the functions to calculate inverted generation distance (IGD).
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

# include "../header/indicators.h"

static struct double_vector *record = NULL;

/* Calculate the IGD value of a population */
void record_igd (void *ptr, int id)
{
    double value;

    if (record == NULL)
    {
        record = (struct double_vector *) malloc (sizeof(struct double_vector));
        record->value = nan("1");
        record->next  = NULL;
    }

    // calculate IGD
    value = calculate_igd (ptr);
    double_vector_pushback (record, value);

    return;
}

void print_igd (char *file_name)
{
    int i;
    double value;
    FILE *fpt;

    fpt = fopen (file_name, "w");

    i = 0;
    while (1)
    {
        value = double_vector_get (record->next, i++);
        if (!isnan(value))
            fprintf (fpt, "%lf\n", value);
        else
            break;
    }

    fclose (fpt);
    double_vector_free (record);
    record = NULL;

    return;
}


/* Calculate IGD value */
double calculate_igd (void *ptr)
{
    int i, j;
    double igd_value;
    double min_distance, cur_distance;

    population_real *pop = (population_real *)ptr;

    igd_value = 0.0;
    for (i = 0; i < PF_size; i++)
    {
        min_distance = euclidian_distance (PF_data[i], pop->ind[0].obj, number_objective);
        for (j = 1; j < popsize; j++)
        {
            cur_distance = euclidian_distance (PF_data[i], pop->ind[j].obj, number_objective);
            if(cur_distance < min_distance)
                min_distance = cur_distance;
        }
        igd_value += min_distance;
    }
    igd_value = igd_value / PF_size;

    return igd_value;
}
