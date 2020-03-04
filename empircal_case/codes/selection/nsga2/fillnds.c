/*
 * fillnds.c:
 *  This file contains the functions to perform non-dominated sorting in NSGA-II and its variants.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization for Learning and Adaptive System (COLA) Laboratory @ University of Exeter
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

/* Environmental selection in NSGA-II */
void fill_nondominated_sort (population_real *new_pop, population_real *mixed_pop)
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
    for (i = 0; i < 2 * popsize; i++)
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
                    temp2 = temp2->child;
                if (flag == -1)
                    end = 1;
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
        if ((archieve_size + front_size) <= popsize)
        {
            do
            {
                copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size++;
                temp2 = temp2->child;
                i++;
            }
            while (temp2 != NULL);
            assign_crowding_distance_indices (new_pop, j, i - 1);
            rank++;
        }
        else
        {
            crowding_fill (mixed_pop, new_pop, i, front_size, elite);
            archieve_size = popsize;
            for (j = i; j < popsize; j++)
                new_pop->ind[j].rank = rank;
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

    // garbage collection
    while (pool != NULL)
    {
        temp1 = pool;
        pool  = pool->child;
        free (temp1);
    }
    while (elite != NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }

    return;
}

/* Environmental selection of r-NSGA-II */
void fill_r_nondominated_sort (population_real *new_pop, population_real *mixed_pop, double sigma)
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
    for (i = 0; i < 2 * popsize; i++)
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
                flag = check_r_dominance (&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]),mixed_pop, sigma);
                if (flag == 1)
                {
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                    temp2 = temp2->child;
                if (flag == -1)
                    end = 1;
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
        if ((archieve_size + front_size) <= popsize)
        {
            do
            {
                copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size++;
                temp2 = temp2->child;
                i++;
            }
            while (temp2 != NULL);
            assign_crowding_distance_indices (new_pop, j, i - 1);
            rank++;
        }
        else
        {
            crowding_fill (mixed_pop, new_pop, i, front_size, elite);
            archieve_size = popsize;
            for (j = i; j < popsize; j++)
                new_pop->ind[j].rank = rank;
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

    // garbage collection
    while (pool != NULL)
    {
        temp1 = pool;
        pool  = pool->child;
        free (temp1);
    }
    while (elite != NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }

    return;
}

/* Environmental selection in R-NSGA-II */
void fill_R_nondominated_sort (population_real *new_pop, population_real *mixed_pop, double* reference_point, double* weights, double epsilon)
{
    int i, j;
    int flag, end;
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
    for (i = 0; i < 2 * popsize; i++)
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
                    temp2 = temp2->child;
                if (flag == -1)
                    end = 1;
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
        if ((archieve_size + front_size) <= popsize)
        {
            do
            {
                copy_ind (&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size++;
                temp2 = temp2->child;
                i++;
            }
            while (temp2 != NULL);
            rank++;
        }
        else
        {
            r_crowding_fill (mixed_pop, new_pop, i, front_size, elite, reference_point, weights, epsilon);
            archieve_size = popsize;
            for (j = i; j < popsize; j++)
                new_pop->ind[j].rank = rank;
        }
        temp2 = elite->child;

        while(elite->child != NULL) {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
    }
    while (archieve_size < popsize);

    // garbage collection
    while (pool != NULL)
    {
        temp1 = pool;
        pool  = pool->child;
        free (temp1);
    }
    while (elite != NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }
}

/* Routine to fill a population with individuals in the decreasing order of crowding distance */
void crowding_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite)
{
    int i, j;
    int *dist;
    list *temp;

    dist = (int *) malloc (front_size * sizeof(int));

    assign_crowding_distance_list (mixed_pop, elite->child, front_size);

    temp = elite->child;
    for (j = 0; j < front_size; j++)
    {
        dist[j] = temp->index;
        temp    = temp->child;
    }
    quicksort_dist (mixed_pop, dist, front_size);
    for (i = count, j = front_size - 1; i < popsize; i++, j--)
        copy_ind(&mixed_pop->ind[dist[j]], &new_pop->ind[i]);

    // garbage collection
    free (dist);

    return;
}

/* Similar to the normal crowding distance assignment procedure in NSGA-II, but considering the preference information
 * by taking the distance to the reference point in calculation */
void r_crowding_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite, double* reference_point, double* weights, double epsilon)
{
    int i, j;
    int temp_idx,remaining_cnt;
    int *front_index;
    double temp, temp_dist;
    double *weighted_distance, *max_obj, *min_obj;
    list *head, *temp_item, *curr_item, *next_item, *r_item;
    individual_real ind, *r_ind, *c_ind;

    front_index       = (int *) malloc (front_size * sizeof(int));
    weighted_distance = (double *) malloc (front_size * sizeof(double));
    max_obj           = (double *) malloc (number_objective * sizeof(double));
    min_obj           = (double *) malloc (number_objective * sizeof(double));

    temp_item = elite->child;
    for (i = 0; i < front_size; i++)
    {
        front_index[i]   = temp_item->index;
        temp_item        = temp_item->child;
    }


    /*for (i = 0; i < number_objective; i++)
    {
        max_obj[i] = -INF;
        min_obj[i] = INF;
        for (j = 0; j < 2 * popsize; j++)
        {
            if (max_obj[i] < (mixed_pop->ind[j].obj[i]))
                max_obj[i] = mixed_pop->ind[j].obj[i];
            if (min_obj[i] > (mixed_pop->ind[j].obj[i]))
                min_obj[i] = mixed_pop->ind[j].obj[i];
        }
    }
*/// Find the maximum and minimum values of each dimension.

    elite = elite->child;
    for (i = 0; i < number_objective; i++)
    {
        max_obj[i] = -INF;
        min_obj[i] = INF;
        for (j = 0; j < front_size; j++)
        {
            ind = mixed_pop->ind[get_item (elite, j)->index];
            if (max_obj[i] < (ind.obj[i]))
                max_obj[i] = ind.obj[i];
            if (min_obj[i] > ind.obj[i])
                min_obj[i] = ind.obj[i];
        }
    }
    // Calculate the weighted Euclidean distance.

    for (i = 0; i < front_size; i++)
    {
        temp_dist = 0.0;
        ind = mixed_pop->ind[get_item (elite, i)->index];
        for (j = 0; j < number_objective; j++)
            temp_dist += weights[j] * pow ((ind.obj[j] - reference_point[j]) / (max_obj[j] - min_obj[j]), 2);
        weighted_distance[i] = sqrt (temp_dist);
        mixed_pop->ind[get_item (elite, i)->index].fitness = INF - weighted_distance[i];
    }
    // Step 3 in the R-NSGA-II paper
    if(epsilon>0)
    {
        head = elite;
        temp_item = head;
        temp_idx = 0;
        while (temp_item != NULL) {
            temp_item->index2 = temp_idx;
            temp_item = temp_item->child;
            temp_idx++;
        }

        remaining_cnt = front_size;
        while (remaining_cnt > 0) {
            int r_idx = rand() % remaining_cnt;    //choose one solution at random.
            r_item = get_item(head, r_idx);
            r_ind = &(mixed_pop->ind[r_item->index]);

            curr_item = head;
            while (curr_item != NULL) {
                c_ind = &(mixed_pop->ind[curr_item->index]);

                temp = 0.0;
                for (i = 0; i < number_objective; i++)
                    temp += fabs(r_ind->obj[i] - c_ind->obj[i]) / (nadir_point[i] - ideal_point[i]);

                next_item = curr_item->child;
                // Elimination of neighboring solutions within 'epsilon' radius
                if (temp < epsilon)
                {
                    if (curr_item->index != r_item->index)
                    {
                        weighted_distance[curr_item->index2] = INF;
                        mixed_pop->ind[curr_item->index].fitness = -INF;
                    }

                   if(curr_item->child != NULL)
                        curr_item->child->parent = curr_item->parent;

                    if(curr_item->parent != NULL)
                        curr_item->parent->child = curr_item->child;

                    if(curr_item->index == head->index)
                        head = next_item;

                    curr_item = next_item;
                    remaining_cnt--;
                }
                else
                    curr_item = next_item;
            }
        }
    }
    // Sort out the solution and eliminate the large value of crowding dist.

   quicksort_dist (mixed_pop, front_index, front_size);
   for (i = count, j = front_size - 1; i < popsize; i++, j--)
        copy_ind(&mixed_pop->ind[front_index[j]], &new_pop->ind[i]);

   free (weighted_distance);
   free (max_obj);
   free (min_obj);
   free (front_index);
   return;
}
