/*
 * fill.c:
 *  This file implements the main function related to the environmental selection of HypE. Note that this code is coded
 *  according to the PISA (http://www.tik.ee.ethz.ch/sop/pisa/?page=pisa.php)
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

void fill_hype_sort (population_real *new_pop, population_real *mixed_pop, int size)
{
    int i, j;
    int flag;
    int end;
    int front_size = 0;
    int archieve_size = 0;
    int rank = 1;

    list *pool;
    list *elite;
    list *temp1, *temp2;

    pool  = (list *) malloc (sizeof(list));
    elite = (list *) malloc (sizeof(list));
    pool->index   = -1;
    pool->parent  = NULL;
    pool->child   = NULL;
    elite->index  = -1;
    elite->parent = NULL;
    elite->child  = NULL;

    temp1 = pool;
    for (i = 0; i < size; i++)
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    i = 0;
    do
    {
        temp1 = pool->child;
        insert (elite, temp1->index);
        front_size = 1;
        temp2 = elite->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        do
        {
            temp2 = elite->child;
            if (temp1 == NULL)
                break;

            do
            {
                end  = 0;
                flag = check_dominance (&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
                if (flag == 1)
                {
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    end = 1;
                }
            }
            while (end != 1 && temp2 != NULL);
            if (flag == 0 || flag == 1)
            {
                insert (elite, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);
        temp2 = elite->child;
        j = i;

        if ( (archieve_size + front_size) <= popsize)
        {
            //printf("\naddF%d:",rank);
            do
            {
                //printf("[%d]->[%d]",temp2->index,i);
                copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size += 1;
                temp2 = temp2->child;
                i += 1;

            }
            while (temp2 != NULL);
            rank += 1;
        }
        else
        {
            hype_fill (mixed_pop, new_pop, i, front_size, elite);
            archieve_size = popsize;
            for (j = i; j < popsize; j++)
            {
                new_pop->ind[j].rank = rank;
            }
        }
        temp2 = elite->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (archieve_size < popsize);

    // free memory
    while (pool!=NULL)
    {
        temp1 = pool;
        pool = pool->child;
        free (temp1);
    }
    while (elite!=NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }

    return;
}

/* check weakly dominate */
int weaklyDominates (double *point1, double *point2, int no_objectives)
{
    int i;
    int better;

    i      = 0;
    better = 1;
    while (i < no_objectives && better)
    {
        better = point1[i] <= point2[i];
        i++;
    }

    return better;
}

/* HypE sampling procedure for Hypervolume approximation */
double hypeSampling (struct double_with_index *val, int size, int nrOfSamples, int param_k, double *rho, void *ptr, list *elite)
{
    int i, s, k;
    int idx, domCount, counter;
    int *hitstat;
    double *sample;
    list * temp;
    population_real * pop = ptr;

    hitstat = malloc (sizeof(int) * size);
    sample  = malloc (sizeof(double) * number_objective);

    for (i = 0; i < size; i++)
        val[i].x = 0.0;
    for (s = 0; s < nrOfSamples; s++)
    {
        for(k = 0; k < number_objective; k++)
            sample[k] = rndreal (ideal_point[k], nadir_point[k]);

        domCount = 0;
        temp     = elite->child;
        counter  = 0;
        do{
            idx = temp->index;

            if (weaklyDominates (pop->ind[idx].obj, sample, number_objective))
            {
                domCount++;
                if(domCount > param_k)
                    break;
                hitstat[counter] = 1;
            }
            else
                hitstat[counter] = 0;
            counter++;
            temp = temp->child;
        } while (temp != NULL);

        if (domCount > 0 && domCount <= param_k)
        {
            temp = elite ->child;
            counter = 0;
            do {
                int idx = temp->index;
                if (hitstat[counter] == 1)
                    val[counter].x += rho[domCount];
                temp = temp->child;
                counter++;
            } while (temp != NULL);
        }
    }
    temp    = elite->child;
    counter = 0;
    do {
        idx = temp->index;

        val[counter].idx = idx;
        val[counter].x = val[counter].x  / (double) nrOfSamples;
        for (k = 0; k < number_objective; k++)
            val[counter].x *= (nadir_point[k]-ideal_point[k]);
        temp = temp->child;
        counter++;
    } while (temp != NULL);

    free (hitstat);
    free (sample);
}

/* Fill the population according to the non-domination levels and remove the individual with the least Hypervolume contribution */
void hype_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite)
{
    int i, j;
    int min_idx, param_k, nrOfSamples;
    double min_value;

    list *temp;
    double *rho;
    struct double_with_index *hv;

    hv  = malloc (front_size * sizeof(struct double_with_index));
    rho = malloc (sizeof(double) * popsize);

    nrOfSamples = 10000;
    while (front_size + count > popsize)
    {
        param_k = front_size + count - popsize;
        // set rho
        rho[0] = 0;
        for (i = 1; i <= param_k; i++)
        {
            rho[i] = 1.0 / (double) i;
            for (j = 1; j <= i - 1; j++)
                rho[i] *= (double) (param_k - j) / (double) (front_size - j);
        }

        hypeSampling (hv, front_size, nrOfSamples, param_k, rho, mixed_pop, elite);

        min_idx   = hv[0].idx;
        min_value = hv[0].x;
        for (i = 1; i < front_size; i++)
        {
            if (hv[i].x < min_value)
            {
                min_idx  = hv[i].idx;
                min_value = hv[i].x;
            }
        }
        // remove min_idx
        temp = elite->child;
        do {
            int idx = temp->index;
            if (min_idx == idx)break;
            temp = temp->child;
        } while (temp != NULL);
        temp->parent->child = temp->child;
        if (temp->child != NULL)
            temp->child->parent = temp->parent;
        free (temp);
        front_size--;

        fflush (stdout);
    }
    i = 0;
    while (count < popsize)
    {
        copy_ind (&(mixed_pop->ind[hv[i].idx]), &(new_pop->ind[count]));
        count++;
        i++;
    }

    // garbage collection
    free (hv);
    free (rho);

    return;
}


