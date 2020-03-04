/*
 * population.c:
 *  This file contains the functions to perform operations upon populations and individuals.
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

# include "../header/population.h"

/* Function to initialize an individual randomly */
void initialize_individual_real (individual_real *ind)
{
    int i;

    if (number_variable != 0)
        for (i = 0; i < number_variable; i++)
            ind->xreal[i] = rndreal (variable_lowerbound[i], variable_upperbound[i]);
    ind->cv = 0;    // initialze the CV of each solution to be 0 (assume all solutions are feasible)

    return;
}

/* Function to initialize a population of individuals */
void initialize_population_real (population_real *pop)
{
    int i;

    for (i = 0; i < popsize; i++)
        initialize_individual_real (&(pop->ind[i]));

    return;
}

/* Function to read static population from external file (Debug use) */
void read_population_real (population_real * pop, char * fileName)
{
    int i, j;
    FILE *popdata = NULL;
    popdata       = fopen (fileName, "r");

    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < number_variable; j++)
        {
            fscanf (popdata, "%lf", &(pop->ind[i].xreal[j]));
        }
    }
    fclose (popdata);

    return;
}

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void copy_ind (individual_real *ind1, individual_real *ind2)
{
    int i;

    ind2->rank    = ind1->rank;
    ind2->fitness = ind1->fitness;
    ind2->cv      = ind1->cv ;
    if (number_variable != 0)
        for (i = 0; i < number_variable; i++)
            ind2->xreal[i] = ind1->xreal[i];

    for (i = 0; i < number_objective; i++)
        ind2->obj[i] = ind1->obj[i];

    return;
}

/* Merge 'pop1' and 'pop2' into a new population 'pop3' */
void merge (population_real *pop1, population_real *pop2, population_real *pop3)
{
    int i, j;

    for (i = 0; i < popsize; i++)
        copy_ind (&(pop1->ind[i]), &(pop3->ind[i]));

    for (i = 0, j = popsize; i < popsize; i++, j++)
        copy_ind (&(pop2->ind[i]), &(pop3->ind[j]));

    return;
}
