/*
 * fitness.c:
 *  This file contains the aggregation function (subproblem) used in MOEA/D.
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

/* Calculate the aggregation function of a solution with respect to a particular subproblem */
double fitnessFunction (individual_real *individual, double *lambda)
{
    int i;
    double fitness, diff, feval, sum, maxFun;

    if (function_type == WS)
    {
        sum = 0.0;
        for (i = 0; i < number_objective; i++)
            sum += lambda[i] * individual->obj[i];

        fitness = sum;
    }
    else if (function_type == N_WS)
    {
        sum = 0.0;
        for (i = 0; i < number_objective; i++)
           sum += lambda[i] * (individual->obj[i] - ideal_point[i]) / (nadir_point[i] - ideal_point[i]);

        fitness = sum;
    }
    else if (function_type == TCH)
    {
        maxFun = -1.0e+30;

        for (i = 0; i < number_objective; i++)
        {
            diff = fabs (individual->obj[i] - ideal_point[i]);
            if (lambda[i] < EPS)
                feval = 0.00001 * diff;
            else
                feval = diff * lambda[i];

            if (feval > maxFun)
                maxFun = feval;
        }

        fitness = maxFun;
    }
    else if (function_type == N_TCH)
    {
        maxFun = -1.0e+30;

        for (i = 0; i < number_objective; i++)
        {
            diff = fabs ((individual->obj[i] - ideal_point[i]) / (nadir_point[i] - ideal_point[i]));
            if (lambda[i] < EPS)
                feval = 0.00001 * diff;
            else
                feval = diff * lambda[i];

            if (feval > maxFun)
                maxFun = feval;
        }

        fitness = maxFun;
    }
    else if (function_type == ITCH)
    {
        maxFun = -1.0e+30;

        for (i = 0; i < number_objective; i++)
        {
            diff = fabs (individual->obj[i] - ideal_point[i]);
            if (lambda[i] < EPS)
                feval = diff / 0.000001;
            else
                feval = diff / lambda[i];

            if (feval > maxFun)
                maxFun = feval;
        }

        fitness = maxFun;
    }
    else if (function_type == N_ITCH)
    {
        maxFun = -1.0e+30;

        for (i = 0; i < number_objective; i++)
        {
            diff = fabs ((individual->obj[i] - ideal_point[i]) / (nadir_point[i] - ideal_point[i]));
            if (lambda[i] < EPS)
                feval = diff / 0.000001;
            else
                feval = diff / lambda[i];

            if (feval > maxFun)
                maxFun = feval;
        }

        fitness = maxFun;
    }
    else if (function_type == PBI)
    {
        double theta;
        double d1, d2, nl;

        theta = 5.0;
        d1 = d2 = nl = 0.0;

        for (i = 0; i < number_objective; i++)
        {
            d1 += (individual->obj[i] - ideal_point[i]) * lambda[i];
            nl += pow (lambda[i], 2.0);
        }
        nl = sqrt (nl);
        d1 = fabs (d1) / nl;

        for (i = 0; i < number_objective; i++)
            d2 += pow ((individual->obj[i] - ideal_point[i]) - d1 * (lambda[i] / nl), 2.0);
        d2 = sqrt (d2);

        fitness = (d1 + theta * d2);
    }
    else if (function_type == N_PBI)
    {
        double theta;
        double d1, d2, nl;

        theta = 5.0;
        d1 = d2 = nl = 0.0;

        for (i = 0; i < number_objective; i++)
        {
            d1 += (individual->obj[i] - ideal_point[i]) / (nadir_point[i] - ideal_point[i]) * lambda[i];
            nl += pow (lambda[i], 2.0);
        }
        nl = sqrt (nl);
        d1 = fabs (d1) / nl;

        for (i = 0; i < number_objective; i++)
            d2 += pow ((individual->obj[i] - ideal_point[i]) / (nadir_point[i] - ideal_point[i]) - d1 * (lambda[i] / nl), 2.0);
        d2 = sqrt (d2);

        fitness = (d1 + theta * d2);
    }
    else
        print_error (1, 1, "Unknown aggregation function type!");

    return fitness;
}
