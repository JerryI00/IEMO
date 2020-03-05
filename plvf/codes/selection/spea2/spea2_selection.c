/*
 * spea2.c:
 *  This file contains the main procedures of the standard SPEA2.
 *
 * Authors:
 *  Qi Xu <qixu.student@gmail.com>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Qi Xu, Renzhi Chen, Ke Li
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

double search_k_minimum (double *a, int k, int size)
{
    int i, j;
    int temp_int;
    double k_min, temp;

    print_error(k > size, 1, "Error in search_k_minimum: k must be smaller or equal to size!");

    for (i = 0; i < k; i++)
    {
        k_min    = a[i];
        temp_int = i;

        for (j = i + 1; j < size; j++)
        {
            if (k_min > a[j])
            {
                k_min    = a[j];
                temp_int = j;
            }
        }
        if (i != temp_int)
        {
            temp        = a[i];
            a[i]        = a[temp_int];
            a[temp_int] = temp;
        }
    }
    return (k_min);
}

void fitness_spea2 (population_real *pop, int total_size, int k_min, int *dominated_Num, int *bedominated_Num, int **dominated_Matrix, int *R_i, double **distance_Matrix, double* D_i, double *kth_distance)
{
    int i, j, flag;

    for (i = 0; i < total_size; i++)
    {
        dominated_Num[i]   = 0;
        bedominated_Num[i] = 0;
        R_i[i]             = 0;
        D_i[i]             = 0.0;
    }

    // Calculate R(i)
    for (i = 0; i < total_size - 1; i++)
    {
        for (j = i + 1; j < total_size; j++)
        {
            flag = check_dominance(&(pop->ind[i]), &(pop->ind[j]));
            dominated_Matrix[i][j] = flag;
            dominated_Matrix[j][i] = -(flag);

            if (flag == 1)
            {
                dominated_Num[i]   += 1;
                bedominated_Num[j] += 1;
            }
            if (flag == -1)
            {
                bedominated_Num[i] += 1;
                dominated_Num[j]   += 1;
            }
        }
    }
    for (i = 0; i < total_size; i++)
    {
        for (j = 0; j < total_size; j++)
        {
            if (dominated_Matrix[i][j] == -1)
            {
                R_i[i] += dominated_Num[j];
            }
        }
    }

    // Calculate D(i)
    for (i = 0; i < total_size; i++)
    {
        distance_Matrix[i][i] = INF;
        for (j = i + 1; j < total_size; j++)
        {
            distance_Matrix[i][j] = distance_Matrix[j][i] = euclidian_distance(pop->ind[i].xreal, pop->ind[j].xreal, number_variable);
        }
    }
    for (i = 0; i < total_size; i++)
    {
        kth_distance[i] = search_k_minimum(distance_Matrix[i], k_min, total_size);
    }
    for (i = 0; i < total_size; i++)
    {
        D_i[i] = 1.0 / (2.0 + kth_distance[i]);
    }
    for (i = 0; i < total_size; i++)
    {
        pop->ind[i].fitness = R_i[i] + D_i[i];
    }
}

void truncate_pop (population_real *tmppop, int tmppop_size, population_real *pop2, int archive_size, double **distance_Matrix)
{
    int i, j, flag;
    list *elist, *temp_list, *plist;
    double temp1, temp2;

    elist         = (list *) malloc(sizeof(list));
    elist->index  = -1;
    elist->parent = NULL;
    elist->child  = NULL;
    temp_list     = elist;

    for (i = 0; i < tmppop_size; i++)
    {
        insert(temp_list, i);
        temp_list = temp_list->child;
    }

    for (i = 0; i < tmppop_size; i++)
    {
        distance_Matrix[i][i] = INF;
        for (j = i + 1; j < tmppop_size; j++)
        {
            distance_Matrix[i][j] = distance_Matrix[j][i] = euclidian_distance(tmppop->ind[i].xreal, tmppop->ind[j].xreal, number_variable);
        }
    }

    for (i = 0; i < tmppop_size - archive_size; i++)
    {
        temp_list = elist->child;

        while (temp_list != NULL)
        {
            flag = 0;
            temp1 = distance_Matrix[temp_list->index][0];

            for (j = 1; j < tmppop_size; j++)
            {
                if (distance_Matrix[temp_list->index][j] < temp1)
                {
                    temp1 = distance_Matrix[temp_list->index][j];
                    flag  = j;
                }
            }
            temp_list->index2 = flag;
            temp_list = temp_list->child;
        }

        temp_list = elist->child;
        temp1 = distance_Matrix[temp_list->index][temp_list->index2];
        temp_list = elist->child;
        plist = temp_list;
        while (temp_list != NULL)
        {
            if (distance_Matrix[temp_list->index][temp_list->index2] < temp1)
            {
                temp1 = distance_Matrix[temp_list->index][temp_list->index2];
                plist = temp_list;
                temp_list = temp_list->child;
            }
            else
            {
                temp_list = temp_list->child;
            }
        }

        if (plist != NULL)
        {
            if (plist->index == 0)
            {
                temp1 = distance_Matrix[plist->index2][plist->index + 1];

                for (j = 2; j < tmppop_size; j++) {
                    if (distance_Matrix[plist->index2][j] < temp1) {
                        temp1 = distance_Matrix[plist->index2][j];
                    }
                }
            }
            else
            {
                temp1 = distance_Matrix[plist->index2][0];
                for (j = 1; j < tmppop_size; j++)
                {
                    if ((j != plist->index) && (distance_Matrix[plist->index2][j] < temp1))
                    {
                        temp1 = distance_Matrix[plist->index2][j];
                    }
                }
            }
            if (plist->index2 == 0)
            {
                temp2 = distance_Matrix[plist->index][plist->index2 + 1];
                for (j = 2; j < tmppop_size; j++)
                {
                    if (distance_Matrix[plist->index][j] < temp2)
                    {
                        temp2 = distance_Matrix[plist->index][j];
                    }
                }
            } else {
                temp2 = distance_Matrix[plist->index2][0];
                for (j = 1; j < tmppop_size; j++)
                {
                    if ((j != plist->index2) && (distance_Matrix[plist->index][j] < temp2))
                    {
                        temp2 = distance_Matrix[plist->index][j];
                    }
                }
            }
            temp_list = elist->child;
            if (temp1 < temp2) {
                while (temp_list != NULL) {
                    if (temp_list->index != plist->index2)
                    {
                        temp_list = temp_list->child;
                    }
                    else
                    {
                        for (j = 0; j < tmppop_size; j++)
                        {
                            distance_Matrix[j][plist->index2] = INF;
                        }
                        del (temp_list);
                        break;
                    }
                }
            }
            else
            {
                while (temp_list != NULL)
                {
                    if (temp_list->index != plist->index)
                    {
                        temp_list = temp_list->child;
                    }
                    else
                    {
                        for (j = 0; j < tmppop_size; j++)
                        {
                            distance_Matrix[j][plist->index] = INF;
                        }
                        del (temp_list);
                        break;
                    }
                }
            }
        }
    }

    temp_list = elist->child;
    i = 0;
    while (temp_list != NULL)
    {
        copy_ind(&(tmppop->ind[temp_list->index]), &(pop2->ind[i]));
        i++;
        temp_list = temp_list->child;
    }
    while (elist != NULL)
    {
        temp_list = elist;
        elist = elist->child;
        free(temp_list);
    }
}

void selection_spea2 (population_real *mixed_pop,int total_size,population_real *archive,int archive_size, individual_real *temp_ind, population_real *temp_pop, double **distance_Matrix)
{
    int i, j, flag;
    int num_nondominated;
    individual_real *ind, *ind1, *ind2;
    double min, temp;

    num_nondominated = 0;
    for (i = 0; i < total_size; i++) // Count num_nondominated
    {
        ind = &(mixed_pop->ind[i]);
        if (ind->fitness < 1.0)
        {
            num_nondominated += 1;
        }
    }
    if (num_nondominated <= archive_size)
    {
        for (i = 0; i < archive_size; i++)
        {
            ind1 = &(mixed_pop->ind[i]);
            min = ind1->fitness;
            flag = i;

            for (j = i + 1; j < total_size; j++)
            {
                ind2 = &(mixed_pop->ind[j]);
                temp = ind2->fitness;
                if (min > temp)
                {
                    flag = j;
                    min  = temp;
                }
            }
            copy_ind (&(mixed_pop->ind[i]), temp_ind);
            copy_ind (&(mixed_pop->ind[flag]), &(mixed_pop->ind[i]));
            copy_ind (temp_ind, &(mixed_pop->ind[flag]));
            copy_ind (&(mixed_pop->ind[i]), &(archive->ind[i]));
        }
    }
    else // Requires truncation: num_nondominated > archive_size
    {
        j = 0;
        for (i = 0; i < total_size; i++)
        {
            if (mixed_pop->ind[i].fitness < 1.0)
            {
                copy_ind (&(mixed_pop->ind[i]), &(temp_pop->ind[j]));
                j += 1;
            }
        }
        truncate_pop (temp_pop, j, archive, archive_size,distance_Matrix);
    }

}









