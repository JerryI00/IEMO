/*
 * nsga3.c:
 *  This file implements the main procedures of NSGA-III. It is based on the following reference:
 *
 * K. Deb, H. Jain, "An Evolutionary Many-Objective Optimization Algorithm Using Reference-Point-Based Nondominated
 * Sorting Approach, Part I: Solving Problems With Box Constraints",
 * IEEE Trans. Evol. Comput. 18(4): 577-601, 2014.
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

int RemainingCount;
int max_rank;
int num_candidates;
int selected_count;
int cur_size;

int *fronts_size;
int *candidate_flag;
int *association_flag;
int *associationCount;
int *lastAssociationCount;
int *ReferenceDirection;
double *maxObjValues;
double *nadirPoint;
double *u;
double *previousMaxValues;
double *maxArr;
double *PerpendicularDistance;
double *intercepts;
double **Z;
double **unitDirections;
double **distanceMatrix;

list **fronts;

population_real *candidate_pop;
population_real *previousExtreme_pop;
population_real *extreme_pop;

/* Association procedure */
void association (population_real *pop, int size)
{
    int i, j, k, min_idx;
    double d1, d2, lam, min_distance;

    // calculate perpendicular distances towards each reference point
    for (i = 0; i < number_weight; i++)
    {
        for (j = 0; j < size; j++)
        {
            d1  = 0.0;
            lam = 0.0;
            for (k = 0; k < number_objective; k++)
            {
                d1 += (pop->ind[j].obj[k] - ideal_point[k]) * lambda[i][k] / intercepts[k];
                lam += lambda[i][k] * lambda[i][k];
            }
            lam = sqrt(lam);
            d1  = d1 / lam;
            d2  = 0.0;
            for (k = 0; k < number_objective; k++)
                d2 += pow(((pop->ind[j].obj[k] - ideal_point[k]) / intercepts[k] - d1 * lambda[i][k] / lam), 2.0);

            // Store the distance in the matrix and in the individual object
            distanceMatrix[j][i] = sqrt(d2);
        }
    }

    // find the closest reference direction to each individual
    for (i = 0; i < size; i++)
    {
        min_idx      = 0;
        min_distance = distanceMatrix[i][0];
        for (j = 1; j < number_weight; j++)
        {
            if (distanceMatrix[i][j] < min_distance)
            {
                min_idx            = j;
                min_distance       = distanceMatrix[i][j];
            }
        }
        ReferenceDirection[i]    = min_idx;
        PerpendicularDistance[i] = min_distance;
    }
}

/* Non-dominated sorting procedure */
void assign_rank (population_real *pop, int size)
{
    int i, j;
    int flag, end;
    int front_size = 0;
    int rank = 0;

    cur_size = 0;

    list *pool;
    list *elite;
    list *temp1, *temp2, *temp3;

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
                flag = check_dominance (&(pop->ind[temp1->index]), &(pop->ind[temp2->index]));
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

        // copy each level into candidate
        if (cur_size < size)
        {
            fronts_size[rank] = front_size;

            if (fronts[rank] != NULL)    // need to clean up
            {
                temp3 = fronts[rank]->child;
                do
                {
                    temp3 = del (temp3);
                    temp3 = temp3->child;
                }
                while (temp3 !=NULL);
            }
            else
            {
                fronts[rank] = malloc (sizeof(list));
            }

            fronts[rank]->index  = -1;
            fronts[rank]->parent = NULL;
            fronts[rank]->child  = NULL;

            temp3    = fronts[rank];
            max_rank = rank;
            temp2    = elite->child;

            do
            {
                pop->ind[temp2->index].rank = rank;
                cur_size++;
                insert (temp3, temp2->index);
                temp2 = temp2->child;
            }
            while (temp2 != NULL);
            rank++;
        }

        temp2 = elite->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (cur_size < size);

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

/* Fill the next parent population one front by the other. Note that if |F_1| + ... + |F_l| is larger than the popsize,
 * we only store F_1 ... F_{l-1} into 'new_pop'; while F_1 ... F_l are stored in 'candidate_pop'. */
void fill_nd_pop (population_real *mixed_pop, population_real *new_pop)
{
    int rank;
    list *temp;

    rank           = 0;
    num_candidates = 0;
    while (num_candidates < popsize && rank < max_rank)
    {
        temp           = fronts[rank]->child;
        selected_count = num_candidates;
        do
        {
            copy_ind (&mixed_pop->ind[temp->index], &candidate_pop->ind[num_candidates]);
            if (num_candidates < popsize)
                copy_ind (&mixed_pop->ind[temp->index], &new_pop->ind[num_candidates]);
            num_candidates++;
            temp = temp->child;

        } while (temp != NULL);
        rank++;
    }
}

void getExtremePoints (population_real *pop, int size)
{
    // calculate unit directions
    int i, j, k;
    int min_idx;
    double max, nextValue;
    double *wDirection;

    for (i = 0; i < number_objective; i++)
    {
        for (j = 0; j < number_objective; j++)
        {
            if (i == j)
                unitDirections[i][j] = 1;
            else
                unitDirections[i][j] = EPS;     // 1e-6
        }
    }
    // re-calculate the previous MAX values of the previous extreme points
    for (i = 0; i < number_objective; i++)
    {
        // set the unit direction (unit direction j)
        wDirection           = unitDirections[i];
        previousMaxValues[i] = (previousExtreme_pop->ind[i].obj[0] - ideal_point[0]) / wDirection[0];
        for (k = 1; k < number_objective; k++)
        {
            nextValue = (previousExtreme_pop->ind[i].obj[k] - ideal_point[k]) / wDirection[k];
            if (nextValue > previousMaxValues[i])
                previousMaxValues[i] = nextValue;
        }
    }

    for (i = 0; i < number_objective; i++)
    {
        wDirection = unitDirections[i];

        // iterate over all the members of the populations
        for ( j = 0; j < size; j++)
        {
            max = (pop->ind[j].obj[0] - ideal_point[0]) / wDirection[0];
            for (k = 1; k < number_objective; k++)
            {
                nextValue = (pop->ind[j].obj[k] - ideal_point[k]) / wDirection[k];
                if (nextValue > max)
                    max = nextValue;
            }
            maxArr[j] = max;
        }
        // select the minimum value out of maxArr
        min_idx = 0;
        for (j = 1; j < size; j++)
        {
            if (maxArr[j] < maxArr[min_idx])
                min_idx = j;
        }

        if (previousMaxValues[i] < maxArr[min_idx])
        {
            /* This means that the previous extreme point was better than the current extreme point and we should retain
             * the previous extreme point instead of replacing it with a new weaker one. extremePoints[i] = previousExtreme_pop[i]; */
            copy_ind (&previousExtreme_pop->ind[i], &extreme_pop->ind[i]);
        }
        else
        {
            /* Now the individual whose index in minIndex in the population is the one representing the extreme factor in the current directions.
             * extremePoints[i] = new Individual(optimizationProblem, individuals[minIndex], individualEvaluator); */
            copy_ind (&pop->ind[min_idx], &extreme_pop->ind[i]);
        }
    }
}

/* Solve the linear system Ax = b */
double* gaussianElimination (double **A, double *b, double *x)
{
    int i, j, p;
    int N, max;
    double alpha, sum, t;
    double *temp;

    N = number_objective;
    for (p = 0; p < N; p++)
    {
        // find pivot row and swap
        max = p;
        for (i = p + 1; i < N; i++)
            if (fabs(A[i][p]) > fabs(A[max][p]))
                max = i;
        temp   = A[p];
        A[p]   = A[max];
        A[max] = temp;
        t      = b[p];
        b[p]   = b[max];
        b[max] = t;

        // singular or nearly singular
        if (fabs(A[p][p]) <= EPS)
            return NULL;

        // pivot within A and b
        for (i = p + 1; i < N; i++)
        {
            alpha = A[i][p] / A[p][p];
            b[i] -= alpha * b[p];
            for ( j = p; j < N; j++)
                A[i][j] -= alpha * A[p][j];
        }
    }

    // back substitution
    for (i = N - 1; i >= 0; i--)
    {
        sum = 0.0;
        for (j = i + 1; j < N; j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

/* Fill the next parent population by adding one member from the last front */
void niching (population_real *mix_pop, population_real *new_pop)
{
    int i, j, rand;
    int min_idx, min_size, direction_idx;
    int selected_count2, newMemberIndex, remainingIndvsCount, minDirClusterSize;
    list *temp1;
    double min_dist;
    struct int_vector *minlist;

    for (i = 0; i < number_weight; i++)
    {
        associationCount[i] = 0;
        association_flag[i] = 1;
        lastAssociationCount[i] = 0;
    }

    for (i = 0; i < selected_count; i++)
        associationCount[ReferenceDirection[i]] = associationCount[ReferenceDirection[i]] + 1;

    for (i = selected_count; i < num_candidates; i++)
    {
        lastAssociationCount[ReferenceDirection[i]] = lastAssociationCount[ReferenceDirection[i]] + 1;
        candidate_flag[i] = 1;
    }

    remainingIndvsCount = popsize - selected_count;
    minlist = malloc (sizeof(struct int_vector));
    minlist->next  = NULL;
    minlist->value = -1;

    selected_count2 = 0;
    while (remainingIndvsCount > 0)
    {
        minDirClusterSize = -1;
        for (i = 0; i < number_weight; i++)
        {
            if (association_flag[i] == 0)
                continue;
            if (associationCount[i] <= minDirClusterSize || minDirClusterSize == -1)
                minDirClusterSize = associationCount[i];
        }
        min_size = 0;
        minlist  = malloc (sizeof(struct int_vector));
        minlist->next  = NULL;
        minlist->value = -1;

        for (i = 0; i < number_weight; i++)
        {
            if (association_flag[i] == 0)
                continue;
            if (associationCount[i] == minDirClusterSize)
            {
                int_vector_pushback (minlist, i) ;
                min_size++;
            }
        }

        rand          = rnd (1, min_size);
        direction_idx = int_vector_get(minlist, rand);

        int_vector_free (minlist);
        if (lastAssociationCount[direction_idx] == 0)
            association_flag[direction_idx] = 0;
        else
        {
            newMemberIndex = -1;
            if (associationCount[direction_idx] == 0)
            {
                // select the one with the smallest perpendicular distance to the current reference vector
                min_idx  = -1;
                min_dist = -1;
                for (i = selected_count; i < num_candidates; i++)
                {
                    if (direction_idx == ReferenceDirection[i] && candidate_flag[i] && (min_idx == -1 || PerpendicularDistance[i] < min_dist))
                    {
                        min_dist = PerpendicularDistance[i];
                        min_idx  = i;
                    }
                }
                newMemberIndex = min_idx;
            }
            else
            {
                for (i = selected_count; i < num_candidates; i++)
                    if (direction_idx == ReferenceDirection[i] && candidate_flag[i])
                        break;
                newMemberIndex = i;
            }
            copy_ind (&mix_pop->ind[newMemberIndex], &new_pop->ind[selected_count + selected_count2]);

            associationCount[direction_idx]     = associationCount[direction_idx] + 1;
            lastAssociationCount[direction_idx] = lastAssociationCount[direction_idx] -1;
            candidate_flag[newMemberIndex]      = 0;
            selected_count2++;
            remainingIndvsCount--;
        }
    }
}

void getIntercepts (population_real *pop, int size)
{
    /* calculate the vector of maximum objective values & the nadir point Initialize the structures & set all their
     * initial values to negative infinity */
    int i, j;
    int flag;

    for (i = 0; i < number_objective;i++)
        Z[i] = malloc (number_objective * sizeof(double));

    for (i = 0; i < number_objective; i++)
    {
        maxObjValues[i] = -EPS ;
        nadirPoint[i] = -EPS ;  //??
        //nadirPoint[i] =  1e4;
    }
    /* traverse all the individuals of the population and get their maximum value of objective (The simplest way of
     * calculating the nadir point is to get these maximum values among the first front individuals) */
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < number_objective; j++)
        {
            if (maxObjValues[j] < pop->ind[i].obj[j] - ideal_point[j])
                maxObjValues[j] = pop->ind[i].obj[j] - ideal_point[j];
            if (pop->ind[i].rank == 0)
            {
                if (nadirPoint[j] < pop->ind[i].obj[j] - ideal_point[j])
                    nadirPoint[j] = pop->ind[i].obj[j] - ideal_point[j];
            }
        }
    }

    // caculate the intercepts, create the hyperplane, prepare your arrays for gaussian elimination
    for (i = 0; i < number_objective; i++)
        for ( j = 0; j < number_objective; j++)
            Z[i][j] = extreme_pop->ind[i].obj[j] - ideal_point[j];

    for (i = 0; i < number_objective; i++)
        u[i] = 1;

    flag = 0;   // false: do not use nadir point
    // Solve the system of equations using gaussian elimination
    if (gaussianElimination (Z, u,intercepts) == NULL)
        flag = 1;

    if (!flag)
    {
        for (i = 0; i < number_objective; i++)
            intercepts[i] = 1 / intercepts[i];
    }
    else // If the follwing condition is true this means that you have to resort to the nadir point
    {
        for (i = 0; i < number_objective; i++)
            intercepts[i] = nadirPoint[i];
    }

    /* If any of the intercepts is still Zero (which means that one of the nadir values is Zero), then use the maximum
     * value of each objective instead (remember that these values were calculated among all the individuals, not just
     * the first-front individuals) */
    for (i = 0; i < number_objective; i++)
    {
        if (intercepts[i] < EPS)
        {
            for (j = 0; j < number_objective; j++)
                intercepts[j] = maxObjValues[j];
            break;
        }
    }
}

void NSGA3 (population_real *parent_pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i;
    int generation;

    nums_weight();
    //initialize_uniform_weight ();
    //read_uniform_weight("1.dat");
    // for assign_rank
    fronts      = (list **) malloc (2 * popsize * sizeof(list *));
    fronts_size = (int *) malloc (2 * popsize * sizeof(int));;
    for (i = 0; i < 2 * popsize; i++)
        fronts[i] = NULL;

    // for gaussian
    maxObjValues = malloc (number_objective * sizeof(double));
    nadirPoint   = malloc (number_objective * sizeof(double));
    u = malloc (number_objective * sizeof(double));
    Z = malloc (number_objective * sizeof(double *));

    // for niching
    ReferenceDirection    = (int *) malloc (2 * popsize * sizeof(int));
    candidate_flag        = (int *) malloc (2 * popsize * sizeof(int));
    association_flag      = (int *) malloc (number_weight * sizeof(int));
    associationCount      = (int *) malloc (number_weight * sizeof(int));
    lastAssociationCount  = (int *) malloc (number_weight * sizeof(int));
    intercepts            = (double *) malloc (number_objective * sizeof(double));
    previousMaxValues     = (double *) malloc (number_objective * sizeof(double));
    maxArr                = (double *) malloc (2 * popsize * sizeof(double));
    PerpendicularDistance = (double *) malloc (2 * popsize * sizeof(double));
    distanceMatrix        = (double **) malloc (sizeof(double *) * 2 * popsize);
    for (i = 0; i < 2 * popsize; i++)
        distanceMatrix[i] = (double *) malloc (sizeof(double) * popsize);

    unitDirections = (double **) malloc (number_objective * sizeof(double *));
    for (i = 0; i < number_objective; i++)
        unitDirections[i] = (double *) malloc (sizeof(double) * popsize);

    previousExtreme_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (previousExtreme_pop, number_objective);
    extreme_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (extreme_pop, number_objective);
    candidate_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (candidate_pop, 2 * popsize);

    generation       = 1;
    evaluation_count = 0;
    printf ("Progress: 1%%");

    initialize_neighborhood ();
    initialize_population_real (parent_pop);
    evaluate_population (parent_pop);
    initialize_idealpoint (parent_pop);

    assign_rank(parent_pop, popsize);

    getExtremePoints (parent_pop, popsize);
    getIntercepts (parent_pop, popsize);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, generation, 0);
    while (evaluation_count < max_evaluation)
    {
        generation++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_real (parent_pop, offspring_pop);

        mutation_real (offspring_pop);

        // population evaluations
        evaluate_population (offspring_pop);

        // environmental selection
        merge (parent_pop, offspring_pop, mixed_pop);
        assign_rank (mixed_pop, 2 * popsize);

        for (i = 0; i < 2 * popsize; i++)
            update_ideal_point (&mixed_pop->ind[i]);

        fill_nd_pop (mixed_pop, parent_pop);

        getExtremePoints (candidate_pop, num_candidates);
        getIntercepts (candidate_pop, num_candidates);

        association (candidate_pop, num_candidates);
        niching (candidate_pop, parent_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    // garbage collection
    free (fronts_size);
    deallocate_memory_pop (candidate_pop, 2 * popsize);
    free (candidate_pop);
    free (associationCount);
    free (lastAssociationCount);
    free (association_flag);
    for(i = 0; i < number_objective; i++)
        free (unitDirections[i]);
    for (i = 0; i < 2 * popsize; i++)
        free (distanceMatrix[i]);
    free (candidate_flag);
    free (unitDirections);
    free (distanceMatrix);
    free (fronts);
    free ( maxObjValues);
    free (nadirPoint);
    free (u);
    free (Z);

    return;
}
