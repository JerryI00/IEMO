/*
 * gd.c:
 *  This file contains the functions to calculate generation distance (GD).
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
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

/* Calculate the GD value of a population */
void record_gd (void *ptr, int id)
{
    double value;

    if (record == NULL)
    {
        record = (struct double_vector *) malloc (sizeof(struct double_vector));
        record->value = nan("1");
        record->next  = NULL;
    }

    // calculate GD
    value = calculate_gd (ptr);
    double_vector_pushback (record, value);

    return;
}

void print_gd (char *file_name)
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


/* Calculate GD value */
double calculate_gd (void *ptr)
{
    int i, j;
    double gd_value;
    double min_distance, cur_distance;

    population_real *pop = (population_real *)ptr;

    gd_value = 0.0;
    for (i = 0; i < popsize; i++)
    {
        min_distance = euclidian_distance (pop->ind[i].obj, PF_data[0], number_objective);
        for (j = 1; j < PF_size; j++)
        {
            cur_distance = euclidian_distance (pop->ind[i].obj, PF_data[j], number_objective);
            if(cur_distance < min_distance)
                min_distance = cur_distance;
        }
        gd_value += min_distance;
    }
    gd_value = gd_value / popsize;

    return gd_value;
}
