/*
 * hv.c:
 *  This file contains the functions to calculate Hypervolume (HV). Here we employ the codes developed by the WFG.
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
# include "../externals/WFG_1.15/wfg.h"

static struct double_vector *record = NULL;

/* Calculate the hv value of a population */
void record_hv (void *ptr, int id)
{
    double value;
    if (record == NULL)
    {
        record = (struct double_vector *) malloc (sizeof(struct double_vector));
        record->value = nan ("1");
        record->next  = NULL;
    }

    // calculate hv
    value = calculate_hv (ptr);
    double_vector_pushback (record, value);

    return;
}

/* Calculate HV value */
double calculate_hv (void *ptr)
{
    double hv_value;

    // call wfg's hv calculation method
    hv_value = hv_wfg (ptr);

    return hv_value;
}

/* Print the HV value into a specified file */
void print_hv (char *file_name)
{
    int i;
    double value;
    FILE *fpt;

    i   = 0;
    fpt = fopen (file_name, "w");
    while (1)
    {
        value = double_vector_get (record->next, i++);
        if (!isnan (value))
            fprintf (fpt, "%lf\n", value);
        else break;
    }
    fclose (fpt);

    // free memory
    double_vector_free (record);
    record = NULL;

    return;
}
