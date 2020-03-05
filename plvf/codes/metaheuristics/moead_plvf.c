/*
 * moead.c:
 *  This file contains the main procedures of the standard MOEAD.
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
 * MERCHANTABILITY or FITNESS FOR central_vector PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received central_vector copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

# include "../header/metaheuristics.h"

double *weights_obj,*reference_point;
//double *reference_point;
void PLVF (population_real *pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i, j ,k = 0;
    int generation;
    int subproblem_id, neighbor_type;
    double *rbf;
    int *rbf_index;
    int best_index=0;

/* common paramters */
    int miu;                                //number of incumbent candidates.
    int tau;                                // number of generations between two consecutive consultation.
    int first_tau;                          //  number of generations in first consecutive consultation.
    double eta;                             //step size.
    double sigma;                            // the width of the Gaussian function.
    generation       = 1;
    evaluation_count = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", run_index);

   nums_weight();
    //initialize_uniform_weight ();


    print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_neighborhood ();
    initialize_population_real (pop);
    evaluate_population (pop);
    initialize_idealpoint (pop);

    track_evolution (pop, generation, 0);
    double *weights_rbf;
    weights_rbf=(double *)malloc(40 * sizeof(double));




    // initialize parameter settings
    permutation = malloc (number_weight * sizeof(int));
    miu = 2 * number_objective + 1;
    tau = 25,first_tau=25;
    eta = 0.05;
    sigma = 0.8;
    individual_real* offspring = &(offspring_pop->ind[0]);
    rbf = (double *) malloc (popsize * sizeof(double));
    double **central_vector;
    central_vector = (double **) malloc (40 * sizeof(double *));
    for (i = 0; i < 40; i++)
        central_vector[i] = (double *) malloc(number_objective * sizeof(double));
    double temp;

    rbf_index = (int *) malloc (10 * sizeof(int));
    double train_fitness[10];
    double **seed;
    seed = (double **) malloc (miu * sizeof(double *));
    for (i = 0; i < miu; i++)
        seed[i] = (double *) malloc(number_objective * sizeof(double));
    double angle[popsize];
    double d_temp;
    while (evaluation_count < max_evaluation)
    {
        print_progress ();
        if (generation == first_tau)
        {

           // set miu "seed" reference points;
            for (j = 0;j < number_objective;j++)
                seed[miu-1][j] = 1.0 / (double) number_objective ;
            for (i = 0;i < number_objective;i++)
            {
                for (j = 0;j < number_objective;j++)
                    seed[i][j]=0;
            }
            for (i = 0;i < number_objective;i++)
            {
                seed[i][i] = 1;
            }

            for (i = number_objective; i <  2 * number_objective; i++)
            {
                for (j = 0;j < number_objective;j++)
                    seed[i][j] = (seed[(i-number_objective)][j]+seed[miu-1][j]) / (2.0);
            }
            //Find  the nearest neighbor from the reference points.
            for(i = 0; i < miu; i++)
            {
                d_temp = -INF;
                for (j = 0;j < popsize; j++)
                {
                    angle[j] = cos_angle(seed[i], pop->ind[j].obj, number_objective);
                    if (angle[j] > d_temp)
                    {
                        d_temp = angle[j];
                        k = j;
                    }
                }
               // rbf[i] = Rbf_fitnessFunction((&pop->ind[k]), weights_obj);
                rbf[i] = cos_angle_inverse(pop->ind[k].obj,reference_point,number_objective);
                for (j = 0; j < number_objective; j++)
                {
                    central_vector[i][j] = pop->ind[k].obj[j];

                }

            }

            rbf_first_train(central_vector, rbf, weights_rbf,sigma,miu);

        }
        //Consultation
        if (generation > first_tau)
        {
            if ((generation) % tau == 0)
            {



                double temp_rbf = INF;
                int temp_index=0 ;
                for (i = 0; i < popsize; i++)
                {
                    rbf[i] = rbf_test((&pop->ind[i]),weights_rbf,central_vector,sigma);
                }

                // Find the best 10 incumbent candidates.
                for (i = 0; i < 10; i++)
                {
                    temp = INF;
                    for (j = 0; j < popsize; j++)
                    {
                        if (rbf[j] < temp)
                        {
                            temp = rbf[j];
                            rbf_index[i] = j;
                        }
                    }
                    rbf[rbf_index[i]] = INF;
                }
                    for (i = 0; i < 10; i++)
                    {
                        best_index = rbf_index[i];
                        lambda[best_index][0] += 10000;
                    }

                //Tune the other weights to the best 10 weights.

                for (i = 0; i < 10; i++)
                {
                    j = rbf_index[i];
                    lambda[j][0] -= 10000;
                    update_neighborhood(lambda[j], 10, i,eta);
                }
                neighbor_size = 20;
                for (i = 0; i < number_weight; i++)
                {
                    for (j = 0; j < number_objective; j++)
                    {
                        lambda[i][j] -= 10000;
                    }

                }
                //RBF neural network training.
                for (i = 0; i < popsize; i++)
                {
                    rbf[i] = rbf_test((&pop->ind[i]),weights_rbf,central_vector,sigma);
                }

                for (i = 0;i < 10; i++)
                {
                    temp_rbf = INF;
                    for(j = 0;j < popsize;j++)
                    {
                        if( rbf[j] < temp_rbf)
                        {
                            temp_rbf = rbf[j];
                            temp_index = j;
                        }
                    }
                    for (k = 0; k < number_objective; k++)
                    {
                        central_vector[i][k] = pop->ind[temp_index].obj[k];

                    }
                    train_fitness[i] = cos_angle_inverse(pop->ind[temp_index].obj,reference_point,number_objective);
                    rbf[temp_index] = INF;
                }

                rbf_train(central_vector, train_fitness, weights_rbf,sigma);
            }
        }
        random_permutation (permutation, number_weight);
        for (i = 0; i < number_weight; i++)
        {
            subproblem_id = permutation[i];

            // crossover and mutation
            crossover_moead_real (pop, offspring, subproblem_id, &neighbor_type);
            mutation_ind (offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point (offspring);

            // update subproblem
            update_subproblem (pop, offspring, subproblem_id, neighbor_type);
        }
        generation++;

        track_evolution (pop, generation, evaluation_count >= max_evaluation);

    }

    free (permutation);
    moead_free ();
    free (rbf_index);
    free (rbf);
    free (weights_rbf);
    for (i = 0; i < 30; i++)
        free (central_vector[i]);
    free (central_vector);
    for (i = 0; i < miu; i++)
        free (seed[i]);
    free (seed);
    return;
}
