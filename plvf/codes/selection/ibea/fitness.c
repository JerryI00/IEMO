/*
 * fitness.c:
 *  This file contains the functions for calculating the fitness values for IBEA.
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

# define rho 1.1
# define kappa 0.05

/* Calculate the epsilon indicator values */
double cal_eps_indicator (individual_real *ind1, individual_real *ind2,population_real *mixed_pop)
{
    int i, j;
    double r, max_eps, temp;

    r = nadir_point[0] - ideal_point[0];
    max_eps = (ind1->obj[0] - ideal_point[0]) / r - (ind2->obj[0] - ideal_point[0]) / r;
    for (i = 1; i < number_objective; i++)
    {
        r = nadir_point[i] - ideal_point[i];
        temp = (ind1->obj[i] - ideal_point[i]) / r - (ind2->obj[i] - ideal_point[i]) / r;
        if (temp > max_eps)
            max_eps = temp;
    }


    return max_eps; //Find the maximum value of eps which satisfies the indicator-formula:I_value.
}

/* Fitness calculation in PBEA */
void pbea_calcFitness (population_real *mixed_pop, double **Indicator_value, int size, double* weights, double* reference_point, double specificity)
{
    int i, j;
    double sum, asf_min, maxAbs_epsIndicatorValue;
    double asf_values[size];
    population_real* pop;

    pop = mixed_pop;

    // determine indicator values and their maximum
    maxAbs_epsIndicatorValue = 0;

    asf_min = INF;
    for (i = 0; i < size; i++)
    {
        asf_values[i] = asf_pbea (pop->ind[i].obj, reference_point, weights);

        if (asf_values[i] < asf_min)
            asf_min = asf_values[i];
    }

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            Indicator_value[i][j] = cal_eps_indicator (&pop->ind[i], &pop->ind[j],mixed_pop);
            Indicator_value[i][j] = Indicator_value[i][j] / (asf_values[j] + specificity - asf_min);

            if (maxAbs_epsIndicatorValue < fabs(Indicator_value[i][j]))
                maxAbs_epsIndicatorValue = fabs(Indicator_value[i][j]);
        }
    }

    // calculate for each pair of individuals the corresponding fitness component
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            Indicator_value[i][j] = -exp((-Indicator_value[i][j] / maxAbs_epsIndicatorValue) / kappa);

    // calculate the solutions' fitness.
    for(i = 0; i < size; i++)
    {
        sum = 0;
        for (j = 0; j < size; j++)
            if (i != j)
                sum += Indicator_value[j][i];
        pop->ind[i].fitness = sum;
    }

    return;
}

/* ASF used in the PBEA */
double asf_pbea (double *a, double *reference_point, double *weights)
{
    int i;
    double sum, ASF_value;
    double temp, weights_s, ASF_max;

    sum     = 0;
    ASF_max = -INF;
    for (i = 0; i < number_objective; i++)
    {
        temp = a[i] - reference_point[i];
        sum += temp;
        weights_s = weights[i] * temp;
        if (ASF_max < weights_s)
            ASF_max = weights_s;
    }
    ASF_value = ASF_max + 0.00001 * sum;

    return ASF_value;
}

void calcFitness (population_real *mixed_pop, double **Indicator_value, int size)
{
    int i, j;
    double sum, maxAbs_epsIndicatorValue;
    population_real* pop;

    pop = mixed_pop;
    // determine indicator values and their maximum
    maxAbs_epsIndicatorValue = 0;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            Indicator_value[i][j] = cal_eps_indicator (&pop->ind[i], &pop->ind[j],mixed_pop);
            if (maxAbs_epsIndicatorValue < fabs (Indicator_value[i][j]))
                maxAbs_epsIndicatorValue = fabs (Indicator_value[i][j]);
        }
    }

    // calculate for each pair of individuals the corresponding fitness component
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            Indicator_value[i][j] = -exp ((-Indicator_value[i][j] / maxAbs_epsIndicatorValue) / kappa);
    // calculate the solutions' fitness.
    for (i = 0; i < size; i++)
    {
        sum = 0;
        for (j = 0; j < size; j++)
            if (i != j)
                sum += Indicator_value[j][i];
        pop->ind[i].fitness = sum;
    }

    return;
}

/* Assign the fitness values to the solutions within a population */
void cal_fitnesses (void *ptr, double **fitcomp, int size)
{
    int i, j;
    double sum;
    population_real *pop = (population_real*) ptr;

    for (i = 0; i < size; i++)
    {
        sum = 0;
        for (j = 0; j < size; j++)
            if (i != j)
                sum += fitcomp[j][i];
        pop->ind[i].fitness = sum;
    }

    return;
}
