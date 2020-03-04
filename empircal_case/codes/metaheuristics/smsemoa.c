/*
 * smsemoa.c:
 *  This file implements the main procedures of SMS-EMOA. It is based on the following reference:
 *
 * N. Beume, B. Naujoks, M. Emmerich, "SMS-EMOA: Multiobjective selection based on dominated hypervolume",
 * European Journal of Operational Research. 181(3): 1653-1669, 2007.
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

# include "../header/metaheuristics.h"

void SMSEMOA (population_real *parent_pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i, j;
    int generation;
    int maxdepth, maxStackSize;
    FILECONTENTS *f = malloc (sizeof(FILECONTENTS));

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    // initialize process
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);
    initialize_nadirpoint (parent_pop);

    // preparation for IWFG algorithm, which is used for calculating the individual Hypervolume contribution
    i_maxn   = number_objective;
    i_maxm   = popsize + 1;
    maxdepth = i_maxn - 2;

    i_fs         = malloc (sizeof(FRONT) * maxdepth);
    partial      = malloc (sizeof(double) * i_maxm);
    heap         = malloc (sizeof(int) * i_maxm);
    stacksize    = malloc (sizeof(int) * i_maxm);
    stacks       = malloc (sizeof(SLICE*) * i_maxm);
    fsorted      = malloc (sizeof(FRONT) * i_maxn);
    torder       = malloc (sizeof(int *) * MAX(i_maxm, i_maxn));
    tcompare     = malloc (sizeof(int *) * i_maxm);
    maxStackSize = MIN(i_maxn - 2, i_slicingDepth (i_maxn)) + 1;
    for (i = 0; i < maxdepth; i++) {
        i_fs[i].points = malloc(sizeof(POINT) * i_maxm);
        for (j = 0; j < i_maxm; j++) {
            i_fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (i_maxn - i - 1));
        }
    }
    for (i = 0; i < i_maxm; i++)
    {
        stacks[i] = malloc (sizeof(SLICE) * maxStackSize);
        for (j = 1; j < maxStackSize; j++)
            stacks[i][j].front.points = malloc (sizeof(POINT) * i_maxm);
    }
    for (i = 0; i < i_maxn; i++)
        fsorted[i].points = malloc(sizeof(POINT) * i_maxm);
    for (i = 0; i < MAX(i_maxn, i_maxm); i++)
        torder[i] = malloc (sizeof(int) * i_maxn);
    for (i = 0; i < i_maxm; i++)
        tcompare[i] = malloc (sizeof(int) * i_maxn);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while (evaluation_count < max_evaluation)
    {
        generation++;
        print_progress ();
        for (i = 0; i < popsize; i++)
        {
            fflush (stdout);
            // reproduction (crossover and mutation)
            crossover_real_steadystate (parent_pop, &(offspring_pop->ind[0]), &(offspring_pop->ind[1]));
            mutation_ind (&(offspring_pop->ind[0]));
            evaluate_individual (&(offspring_pop->ind[0]));

            // update nadir point
            update_nadir_point (&(offspring_pop->ind[0]));

            // environmental selection
            merge (parent_pop, offspring_pop, mixed_pop);
            fill_hv_sort (f, parent_pop, mixed_pop, popsize + 1);
        }
            // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    // garbage collection
    for (i = 0; i < maxdepth; i++)
    {
        for (j = 0; j < i_maxm; j++)
            free (i_fs[i].points[j].objectives);
        free (i_fs[i].points);
    }
    free (i_fs);

    for (i = 0; i < i_maxm; i++)
        free (stacks[i]);

    free (partial);
    free (heap);
    free (stacksize);
    free (stacks);

    for (i = 0; i < i_maxn; i++)
        free (fsorted[i].points);
    free (fsorted);
    for (i = 0; i < MAX(i_maxn, i_maxm); i++)
        free (torder[i]);
    for (i = 0; i < i_maxm; i++)
        free (tcompare[i]);

    free (torder);
    free (tcompare);

    return;
}
