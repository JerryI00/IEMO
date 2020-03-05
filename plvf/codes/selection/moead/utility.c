/*
 * utility.c:
 *  This file contains some utility functions used in MOEA/D for initialization and memory control related.
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



double** my_pinv(double** n, int row, int col);
int C;
#define NS 10
/* Memory release for MOEA/D related pointers */
void moead_free ()
{
    int i;

    free (ideal_point);

    for (i = 0; i < number_weight; i++)
        free (lambda[i]);
    free (lambda);

    for (i = 0 ; i < popsize ; i++)
        free(neighborhood[i]);
    free(neighborhood);

    lambda       = NULL;
    ideal_point  = NULL;
    neighborhood = NULL;

    return;
}

void uniform_random_weights_const(double **lambda, int num_weights, int dimensions, double lower_bound, double upper_bound) {
    for(int i = 0; i < num_weights; i++) {

        for(int j = 0; j < dimensions; j++) {
            double scale = rand() / (double) RAND_MAX;
            lambda[i][j] = scale * (upper_bound - lower_bound) + lower_bound;
        }

        normalise_vector(lambda[i], dimensions);
    }
}

void uniform_random_weights(double **lambda, int num_weights, int dimensions, const double* lower_bounds, const double* upper_bounds) {
    for(int i = 0; i < num_weights; i++) {

        for(int j = 0; j < dimensions; j++) {
            double scale = rand() / (double) RAND_MAX;
            lambda[i][j] = scale * (upper_bounds[j] - lower_bounds[j]) + lower_bounds[j];
        }

        normalise_vector(lambda[i], dimensions);
    }
}

void initialize_uniform_weight ()
{
    if(weight_file != NULL)
    {
        read_uniform_weight(weight_file);
        free(weight_file);
        return;
    }
    int i, j, l;
    int ptr;
    int layer;

    int layer_size;

    double shrink;
    double *Vec;


    number_weight = 0;

    int gaps = 1;
    while(1)
    {
        layer_size  = combination (number_objective + gaps - 1, gaps);
        //printf("[%d]%d\n",gaps,layer_size);
        if(layer_size>popsize) break;
        gaps = gaps + 1;
        number_weight = layer_size;
    }
    gaps = gaps - 1;
    print_error(gaps == 0,1, "popsize is too small, cannot generate weight");
    lambda = (double **) malloc (number_weight * sizeof(double *));
    for (i = 0; i < number_weight; i++)
        lambda[i] = (double *) malloc(number_objective * sizeof(double));

    C   = 0;

    Vec = (double *) malloc (number_objective * sizeof(double));
    for (i = 0; i < number_objective; i++)
    Vec[i] = 0;
    set_weight (Vec, gaps, 0, number_objective);
    for (i = 0; i < number_weight; i++)
        for (j = 0; j < number_objective; j++) {
            lambda[i][j] = lambda[i][j] / gaps;
        }
    free (Vec);


    return;
}
/* Initialize the weight vectors according to the Das and Dennis' method */
void initialize_layers_weight ()
{
    int i, j, l;
    int ptr;
    int layer;

    int *layer_size;
    int *gaps_table;

    double shrink;
    double *Vec;

    gaps_table = weight_gaps_table[number_objective];
    for (layer = 0; layer < number_objective; layer++)
        if (gaps_table[layer] <= 0)
            break;

    number_weight = 0;
    layer_size = (int *) malloc (sizeof(int) * layer);
    for (i = 0; i < layer; i++)
    {
        layer_size[i]  = combination (number_objective + gaps_table[i] - 1, gaps_table[i]);
        number_weight = number_weight + layer_size[i];
    }

    lambda = (double **) malloc (number_weight * sizeof(double *));
    for (i = 0; i < number_weight; i++)
        lambda[i] = (double *) malloc(number_objective * sizeof(double));

    C   = 0;
    ptr = 0;

    shrink = 1;
    for (l = 0; l < layer; l++) {
        Vec = (double *) malloc (number_objective * sizeof(double));
        for (i = 0; i < number_objective; i++)
            Vec[i] = 0;
        set_weight (Vec, gaps_table[l], 0, number_objective);
        for (i = ptr; i < ptr + layer_size[l]; i++)
            for (j = 0; j < number_objective; j++) {
                lambda[i][j] = lambda[i][j] / gaps_table[l];
                lambda[i][j] = (1 - shrink) / number_objective + shrink * lambda[i][j];
            }
        ptr    = ptr + layer_size[l];
        shrink = shrink * 0.8;
        free (Vec);
    }
    free (layer_size);

    return;
}

void initialize_miu_weight (int miu)
{
    if(weight_file != NULL)
    {
        read_uniform_weight(weight_file);
        free(weight_file);
        return;
    }
    int i, j, l;
    int ptr;
    int layer;
    int number_miu_weight;

    int layer_size;

    double shrink;
    double *Vec;


    number_miu_weight = 0;

    int gaps = 1;
    while(1)
    {
        layer_size  = combination (number_objective + gaps - 1, gaps);
        //printf("[%d]%d\n",gaps,layer_size);
        if(layer_size>miu) break;
        gaps = gaps + 1;
        number_miu_weight = layer_size;
    }
    gaps = gaps - 1;
    print_error(gaps == 0,1, "miu is too small, cannot generate weight");
    rbf_lambda = (double **) malloc (number_miu_weight * sizeof(double *));
    for (i = 0; i < number_miu_weight; i++)
        rbf_lambda[i] = (double *) malloc(number_objective * sizeof(double));

    C   = 0;

    Vec = (double *) malloc (number_objective * sizeof(double));
    for (i = 0; i < number_objective; i++)
        Vec[i] = 0;
    rbf_set_weight (Vec, gaps, 0, number_objective);
    for (i = 0; i < number_miu_weight; i++)
    {
        for (j = 0; j < number_objective; j++)
        {
            rbf_lambda[i][j] = rbf_lambda[i][j] / gaps;
            printf("%lf\t", rbf_lambda[i][j]);
        }
        printf("\n");
    }
    free (Vec);


    return;
}
/* Read weight vectors from an external file */
void read_uniform_weight (char *file)
{
    int i, j;

    FILE *weight  = NULL;
    number_weight = popsize;
    weight = fopen (file, "r");
    printf("weight file:%s",file);

    lambda = (double **) malloc (number_weight * sizeof(double *));
    for (i = 0; i < number_weight; i++)
        lambda[i] = (double *) malloc (number_objective * sizeof(double));

    for (i = 0; i < number_weight; i++)
        for (j = 0; j < number_objective; j++)
            fscanf (weight, "%lf", &(lambda[i][j]));
    return;
}

/* Sample uniformly distributed weight vectors according to Das and Denis' method */
void set_weight (double *weight, double unit, double sum, int dim)
{
    int i;

    if (dim == number_objective)
    {
        for ( i = 0; i < number_objective; i++)
            weight[i] = 0;
    }

    if (dim == 1)
    {
        weight[0] = unit - sum;
        for ( i = 0; i < number_objective; i++)
            lambda[C][i] = weight[i];
        C = C + 1;
        return;
    }
    for (i = 0; i <= unit - sum; i++)
    {
        weight[dim - 1] = i;
        set_weight (weight, unit, sum + i, dim - 1);
    }

    return;
}

void rbf_set_weight (double *weight, double unit, double sum, int dim)
{
    int i;

    if (dim == number_objective)
    {
        for ( i = 0; i < number_objective; i++)
            weight[i] = 0;
    }

    if (dim == 1)
    {
        weight[0] = unit - sum;
        for ( i = 0; i < number_objective; i++)
            rbf_lambda[C][i] = weight[i];
        C = C + 1;
        return;
    }
    for (i = 0; i <= unit - sum; i++)
    {
        weight[dim - 1] = i;
        rbf_set_weight (weight, unit, sum + i, dim - 1);
    }

    return;
}

/* Initialize the neighborhood structure of weight vectors */
void initialize_neighborhood ()
{
    int i, j, k;
    struct double_with_index *dist;

    dist         = (struct double_with_index*) malloc (sizeof(struct double_with_index) * number_weight);
    neighborhood = (int **) malloc (popsize * sizeof(int *));   // free in moead_free

    for (i = 0; i < popsize; i++)
        neighborhood[i] = (int *) malloc (neighbor_size * sizeof(int));
    for (i = 0; i < popsize; i++)
    {
        int id = i % number_weight;
        // calculate the distances based on weight vectors
        for (j = 0; j < number_weight; j++)
        {
            dist[j].x   = euclidian_distance (lambda[id], lambda[j], number_objective);
            dist[j].idx = j;
        }
        for (k = 0; k < neighbor_size; k++)
        {
            for (j = k + 1; j < number_weight; j++)
            {
                if (dist[k].x - dist[j].x >EPS)
                {
                    double temp = dist[k].x;
                    dist[k].x   = dist[j].x;
                    dist[j].x   = temp;
                    int id = dist[k].idx;
                    dist[k].idx = dist[j].idx;
                    dist[j].idx  = id;
                }
            }
        }
        for(j = 0; j < neighbor_size; j ++)
            neighborhood[i][j] = dist[j].idx;
    }

    free (dist);

    return;
}

/* Tournament selection to pick the most active subproblems to evolve (based on the utility) */
void tour_selection_subproblem (int depth)
{
    int i;
    int selected_size, candidate_size;
    int best_candidate_id, best_subproblem_id;
    int random, temp_candidate_id;

    for (i = 0; i < number_objective; i++)
        int_vector_pushback (selected, i);
    selected_size = int_vector_size (selected);

    for (i = number_objective; i < number_weight; i++)
        int_vector_pushback (candidate, i);
    candidate_size = int_vector_size (candidate);

    while (selected_size < (int) (number_weight / 5.0))
    {
        best_candidate_id  = (int) (rndreal(0,1) * candidate_size);
        best_subproblem_id = int_vector_get (candidate, best_candidate_id + 1);
        for (i = 1; i < depth; i++)
        {
            random            = (int) (rndreal (0, 1) * candidate_size);
            temp_candidate_id = int_vector_get (candidate, random + 1);
            if (utility[temp_candidate_id] > utility[best_subproblem_id])
            {
                best_candidate_id  = random;
                best_subproblem_id = temp_candidate_id;
            }
        }
        int_vector_pushback (selected, best_subproblem_id);
        selected_size++;

        int_vector_remove (candidate, best_candidate_id + 1);
        candidate_size--;
    }
}

/* Update the utility of each subproblem */
void comp_utility (population_real *pop, population_real *saved_values)
{
    int i;
    double f1, f2, cur_utility, delta;

    for (i = 0; i < popsize; i++)
    {
        f1 = fitnessFunction (&(pop->ind[i]), lambda[i]);
        f2 = fitnessFunction (&(saved_values->ind[i]), lambda[i]);

        delta = f2 - f1;
        if (delta > 0.001)
            utility[i] = 1.0;
        else
        {
            cur_utility = (0.95 + (0.05 * delta / 0.001)) * utility[i];
            utility[i] = cur_utility < 1.0 ? cur_utility : 1.0;
        }
        copy_ind (&(pop->ind[i]), &(saved_values->ind[i]));
    }

    return;
}
void nums_weight()
{
    int i,j;
    number_weight = 275;
    lambda = (double **) malloc (number_weight * sizeof(double *));
    for (i = 0; i < number_weight; i++)
        lambda[i] = (double *) malloc(number_objective * sizeof(double));
    FILE *fp;
    if ((fp = fopen ("m10.txt","r"))==NULL)
    {
        printf("error\n");
        exit(1);
    }
    for (i = 0; i < number_weight; i++)
        for (j = 0; j < number_objective; j++) {
            fscanf (fp,"%lf",&(lambda[i][j]));

        }
}

void update_neighborhood (double *weight,int miu,int number,double eta)
{
    int i, j, k;
    struct double_with_index *dist;


    if((number_weight - 10) % 10 == 0)
        neighbor_size = (number_weight - 10) / 10;
    else
        neighbor_size = (number_weight - 10) / 10 + 1;
    if ( (number==9) && (number_weight - 10) % 10 != 0)
        neighbor_size = number_weight -10- 9 * neighbor_size;

    int *neighbor_index1;
    neighbor_index1 = (int *) malloc (neighbor_size * sizeof(int));
    dist         = (struct double_with_index*) malloc (sizeof(struct double_with_index) * number_weight);

        for (j = 0; j < number_weight; j++)
        {
            dist[j].x   = euclidian_distance (weight, lambda[j], number_objective);
            if (dist[j].x == 0)
            {
                dist[j].x = INF;
            }
            dist[j].idx = j;

        }

        for (k = 0; k < neighbor_size; k++)
        {
            double temp = INF;
            int n;
            for (j = 0; j < number_weight; j++)
            {

                if (dist[j].x < temp)
                {
                    temp = dist[j].x;
                    neighbor_index1[k] = dist[j].idx;
                    n =j;
                }

            }
            dist[n].x = 10000;
        }
        for(i = 0;i < number_objective;i++)
        {
            for (j = 0; j < neighbor_size; j++)
            {
                lambda[neighbor_index1[j]][i] +=  eta * (weight[i] - lambda[neighbor_index1[j]][i]) +10000 ;
            }

           weight[i] +=10000;
        }


    free (dist);
    free(neighbor_index1);

}
   void rbf_train(double **a,double *rbf,double *weights_rbf,double sigma)
{
    int i,j,k;
    double **active_function;
    active_function = (double **)malloc(sizeof(double *) * 10);
    for(i = 0; i < 10; i++)
    {
        active_function[i] = (double *)malloc(sizeof(double) * 10);
    }


    double *center;
    center = (double *)malloc(number_objective * sizeof(double));


     for(i = 0; i < 10; i++)
    {
         for(j = 0;j < 10;j++)
        active_function[i][j] = 0;
    }

      for (i = 0;i < 10;i++)
    {
        for (j = 0; j < number_objective; j++)
            center[j] = a[i][j];

        for (j = 0; j < 10; j++)
        {
            for (k = 0; k < number_objective; k++)
            {
                active_function[j][i] += pow(center[k] - a[j][k], 2);
            }
            active_function[j][i] = exp(sqrt(active_function[j][i])*(-sigma));
           // printf("%lf\t",active_function[j][i]);
        }
      //printf("\n");
    }
    active_function = my_pinv(active_function,10,10);

       for(i = 0; i < 10; i++)
       {
           weights_rbf[i]=0.0;
           for(j = 0;j < 10;j++)
              weights_rbf[i] += rbf[j] * active_function[i][j];

       }

}
double rbf_test(individual_real *individual,double *weights_rbf,double **a,double sigma)
{
    int i,j,k;
    double *active_function;
    double rbf_fitness;
    rbf_fitness = 0.0;
    active_function = (double *)malloc(10 * sizeof(double));
    double *center;
    center = (double *)malloc(number_objective * sizeof(double));


    for(i = 0; i < 10; i++)
    {

            active_function[i] = 0;
    }

    for (i = 0;i < 10;i++)
    {
        for (j = 0; j < number_objective; j++)
        {
            center[j] = a[i][j];
        }

            for (k = 0; k < number_objective; k++)
            {
                active_function[i] += pow(center[k] - individual->obj[k], 2);
            }
     active_function[i] = exp(sqrt(active_function[i])*(-sigma));

    }
    for (i = 0;i < 10;i++)
        rbf_fitness += weights_rbf[i] * active_function[i];
    return rbf_fitness;
}

void rbf_first_train(double **a,double *rbf,double *weights_rbf,double sigma,int miu)
{
    int i,j,k;
    double **active_function;
    active_function = (double **)malloc(sizeof(double *) * miu);
    for(i = 0; i < miu; i++)
    {
        active_function[i] = (double *)malloc(sizeof(double) * miu);
    }

    double *center;
    center = (double *)malloc(number_objective * sizeof(double));


    for(i = 0; i < miu; i++)
    {
        for(j = 0;j < miu;j++)
            active_function[i][j] = 0;
    }

    for (i = 0;i < miu;i++)
    {
        for (j = 0; j < number_objective; j++)
            center[j] = a[i][j];

        for (j = 0; j < miu; j++)
        {
            for (k = 0; k < number_objective; k++)
            {
                active_function[j][i] += pow(center[k] - a[j][k], 2);
            }
            active_function[j][i] = exp(sqrt(active_function[j][i])*(-sigma));


        }

    }

    active_function = my_pinv(active_function,miu,miu);
    for(i = 0; i < miu; i++)
    {
        weights_rbf[i]=0.0;
        for(j = 0;j < miu;j++)
            weights_rbf[i] += rbf[j] * active_function[i][j];
    }

}