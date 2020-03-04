/*
 * fillnds_hv.c:
 *  This file contains the functions to perform non-dominated sorting and environmental selection in SMS-EMOA.
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
#include "../../header/global.h"

void fill_constraint_hv_sort (FILECONTENTS *f, population_real* new_pop, population_real* mixed_pop, int size)
{
    int i, j;
    int flag;
    int end;
    int front_size = 0;
    int archieve_size = 0;
    int rank = 1;
    int cv_count = 0;

    list *pool;
    list *elite;
    list *temp1, *temp2;

    pool  = (list *)malloc(sizeof(list));
    elite = (list *)malloc(sizeof(list));
    pool->index   = -1;
    pool->parent  = NULL;
    pool->child   = NULL;
    elite->index  = -1;
    elite->parent = NULL;
    elite->child  = NULL;

    temp1 = pool;
    for (i = 0; i < size; i++)
    {
        if ( mixed_pop->ind[i].cv > -EPS) // cv = 0
        {
            cv_count ++;
            insert (temp1,i);
            temp1 = temp1->child;
            //  printf("%d ",i);
        }
    }
    i = 0;
    do
    {
        if(cv_count > popsize)
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
                //printf("flag%d\n",i);
                hv_fill (f,mixed_pop, new_pop, i, front_size, elite);
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
        else
        {
            temp1 = pool->child;
            archieve_size = 0;
            while(temp1!=NULL)
            {
                //printf("cp(a):%d->%d,%d",temp1->index,archieve_size,mixed_pop->ind[temp1->index].rank);
                copy_ind(&mixed_pop->ind[temp1->index], &new_pop->ind[archieve_size]);
                archieve_size ++;
                temp1 = temp1->child;
            }
            struct double_with_index* temp;
            int c = 0;
            temp  = malloc (sizeof(struct double_with_index) * (2 * popsize - archieve_size));
            for (i = 0; i < 2 * popsize; i++)
            {
                if(mixed_pop->ind[i].cv < -EPS)
                {
                    temp[c].idx = i;
                    temp[c].x   = mixed_pop->ind[i].cv;
                    c++;
                }
            }
            c = 0;
            qsort (temp, 2 * popsize - archieve_size, sizeof(struct double_with_index), double_with_index_smaller_cmp);
            while (archieve_size < popsize)
            {
                // printf("cp(a):%d->%d,%lf",temp[c].idx,archieve_size,temp[c].x);
                copy_ind (&mixed_pop->ind[temp[c].idx], &new_pop->ind[archieve_size]);
                c++;
                archieve_size++;
            }
        }

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

void fill_hv_sort (FILECONTENTS *f, population_real* new_pop, population_real* mixed_pop, int size)
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

    pool  = (list *)malloc(sizeof(list));
    elite = (list *)malloc(sizeof(list));
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
            //printf("flag%d\n",i);
            hv_fill (f,mixed_pop, new_pop, i, front_size, elite);
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

/* Fill the population according to the non-domination levels and remove the individual with the least Hypervolume contribution */
void hv_fill (FILECONTENTS *f,population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite)
{
    int i;
    int id, num_same;
    int *dist;
    list *temp;

    dist = (int *) malloc (front_size * sizeof(int));

    i_read_data (f, mixed_pop, elite->child, front_size);

    i_n = f->fronts[0].n;

    //printf("i_n:%d\n",i_n);

    double eh[i_n + 2];
    if (i_n == 2)
        i_ihv2 (f->fronts[0], eh);
    else
        i_ihv (f->fronts[0], eh);

    temp = elite->child;
    do
    {
        id = temp->index;
        num_same = 0;
        for (i = 0; i < number_objective; i++)
        {
            if (fabs ((nadir_point[i] - eh[i]) - mixed_pop->ind[id].obj[i]) < 1e-4)
                num_same++;
        }
        if (num_same == number_objective)
            break;
        temp = temp->child;
    } while (temp != NULL);

    list *temp2 = elite->child;
    do
    {
        if (temp2->index != id)
        {
            copy_ind (&(mixed_pop->ind[temp2->index]), &(new_pop->ind[count]));
            count++;
        }
        temp2 = temp2->child;
    } while (temp2 != NULL);

    // free memory
    free_file_content (f);
    free (dist);

    return;
}

