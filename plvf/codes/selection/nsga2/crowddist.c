/*
 * crowding_distance.c:
 *  This file contains the functions to perform crowding distance calculation in NSGA-II.
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

/* Routine to compute crowding distances */
void assign_crowding_distance (population_real *pop, int *dist, int **obj_array, int front_size)
{
    int i, j;

    for (i = 0; i < number_objective; i++)
    {
        for (j = 0; j < front_size; j++)
            obj_array[i][j] = dist[j];
        //在每一个目标i上进行排序，使用obj_array[i]返回排序结果
        quicksort_front_obj (pop, i, obj_array[i], front_size);
    }

    for (j = 0; j < front_size; j++)
        pop->ind[dist[j]].fitness = 0.0;

    for (i=0; i<number_objective; i++)
    {
        //对每一个目的边界点的拥挤距离给予一个无穷大数，使在多样性保持时能保留边界点
        pop->ind[obj_array[i][0]].fitness = INF;
        pop->ind[obj_array[i][front_size - 1]].fitness = INF;
    }
    for (i = 0; i < number_objective; i++)
    {
        for (j = 1; j < front_size-1; j++)
        {
            if (pop->ind[obj_array[i][j]].fitness != INF)
            {
                if (pop->ind[obj_array[i][front_size - 1]].obj[i] == pop->ind[obj_array[i][0]].obj[i])
                {
                    pop->ind[obj_array[i][j]].fitness += 0.0;
                }
                else
                {
                    pop->ind[obj_array[i][j]].fitness += (pop->ind[obj_array[i][j + 1]].obj[i] - pop->ind[obj_array[i][j - 1]].obj[i]) / (pop->ind[obj_array[i][front_size - 1]].obj[i] - pop->ind[obj_array[i][0]].obj[i]);
                }
            }
        }
    }
    for (j = 0; j < front_size; j++)
    {
        if (pop->ind[dist[j]].fitness != INF)
        {
            pop->ind[dist[j]].fitness = (pop->ind[dist[j]].fitness) / number_objective;
        }
    }

    return;
}

/* Routine to compute crowding distance based on ojbective function values when the population_real in in the form of a list */
void assign_crowding_distance_list (population_real *pop, list *lst, int front_size)
//lst为指向非支配集序号数组的指针
{
    int i;
    int *dist;
    int **obj_array;
    list *temp;

    temp = lst;
    if (front_size == 1)
    {
        pop->ind[lst->index].fitness = INF;
        return;
    }
    if (front_size == 2)
    {
        pop->ind[lst->index].fitness = INF;
        pop->ind[lst->child->index].fitness = INF;
        return;
    }
    obj_array = (int **)malloc(number_objective*sizeof(int *));
	//dist用于存储个体的原始序号，以供后面排序传参之用
    dist = (int *)malloc(front_size*sizeof(int));
    for (i = 0; i < number_objective; i++)
    {
		//obj_array[i]用存储个体在目标i上的个体排序序号
        obj_array[i] = (int *)malloc(front_size*sizeof(int));
    }
    for (i = 0; i < front_size; i++)
    {
        dist[i] = temp->index;
        temp    = temp->child;
    }
    assign_crowding_distance (pop, dist, obj_array, front_size);

    // free memory
    free (dist);
    for (i = 0; i < number_objective; i++)
        free (obj_array[i]);
    free (obj_array);

    return;
}

/* Routine to compute crowding distance based on objective function values when the population_real in in the form of an array */
void assign_crowding_distance_indices (population_real *pop, int c1, int c2)
{
    int i;
    int front_size;

    int *dist;
    int **obj_array;

    front_size = c2 - c1 + 1;
    if (front_size == 1)
    {
        pop->ind[c1].fitness = INF;
        return;
    }
    if (front_size == 2)
    {
        pop->ind[c1].fitness = INF;
        pop->ind[c2].fitness = INF;
        return;
    }
    obj_array = (int **)malloc(number_objective * sizeof(int *));
    dist      = (int *)malloc(front_size * sizeof(int));
    for (i = 0; i < number_objective; i++)
    {
        obj_array[i] = (int *)malloc(front_size * sizeof(int));
    }
    for (i = 0; i < front_size; i++)
    {
        dist[i] = c1++;
    }
    assign_crowding_distance (pop, dist, obj_array, front_size);

    // free memory
    free (dist);
    for (i = 0; i < number_objective; i++)
        free (obj_array[i]);
    free (obj_array);

    return;
}
