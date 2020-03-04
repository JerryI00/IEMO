/*
 * qsort.c:
 *  This file contains the functions to perform quick sort for NSGA-II's crowding distance.
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

# include "../../header/global.h"
# include "../../header/rand.h"
# include "../../header/selection.h"

/* Randomized quick sort routine to sort a population based on a particular objective chosen */
void quicksort_front_obj (population_real *pop, int objcount, int obj_array[], int obj_array_size)
{
    q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size - 1);

    return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
void q_sort_front_obj (population_real *pop, int objcount, int obj_array[], int left, int right)
{
    int i, j;
    int index;
    int temp;
    double pivot;

    if (left < right)
    {
        index = rnd (left, right);
		// obj_array[right]������ʱ�洢�������������ŵĸ���
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        pivot = pop->ind[obj_array[right]].obj[objcount];
        i = left - 1;		// i��¼��Ŀ��ֵ��pivotС�ĸ�����Ŀ
        for (j = left; j < right; j++)
        {
            if (pop->ind[obj_array[j]].obj[objcount] <= pivot)
            {
                i += 1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }
        index = i + 1;
		// �ҵ���pivotС�ĸ�����i������˿��Խ� �������������ŵĸ��� ����i+1��λ�ã�
		// ���ø��洢��obj_array[right]�У���˿��Խ�obj_array[i+1]��obj_array[right]�������ﵽĿ��
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        q_sort_front_obj (pop, objcount, obj_array, left, index-1);
        q_sort_front_obj (pop, objcount, obj_array, index+1, right);
    }

    return;
}

/* Randomized quick sort routine to sort a population based on crowding distance */
void quicksort_dist(population_real *pop, int *dist, int front_size)
{
    q_sort_dist (pop, dist, 0, front_size - 1);

    return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
void q_sort_dist(population_real *pop, int *dist, int left, int right)
{
    int i, j;
    int index;
    int temp;
    double pivot;

    if (left < right)
    {
        index = rnd (left, right);
        temp  = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = pop->ind[dist[right]].fitness;

        i = left - 1;
        for (j = left; j < right; j++)
        {
            if (pop->ind[dist[j]].fitness <= pivot)
            {
                i += 1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            }
        }
        index = i + 1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        q_sort_dist (pop, dist, left, index-1);
        q_sort_dist (pop, dist, index+1, right);
    }
    return;
}
