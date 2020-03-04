/*
 * sbx.c:
 *  This file contains the functions to perform Simulated Binary Crossover (SBX) operation.
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

# include "../../header/reproduction.h"

/* Routine for real variable SBX crossover */
void sbx_crossover (individual_real *parent1, individual_real *parent2, individual_real *child1, individual_real *child2)
{
    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;

    if (randomperc() <= pcross_real)
    {

        for (i = 0; i < number_variable; i++)
        {
            if (randomperc() <= 0.5)
            {
                if (fabs (parent1->xreal[i]-parent2->xreal[i]) > 1e-9)
                {
                    if (parent1->xreal[i] < parent2->xreal[i])
                    {
                        y1 = parent1->xreal[i];
                        y2 = parent2->xreal[i];
                    }
                    else
                    {
                        y1 = parent2->xreal[i];
                        y2 = parent1->xreal[i];
                    }
                    yl    = variable_lowerbound[i];
                    yu    = variable_upperbound[i];
                    rand  = randomperc ();
                    beta  = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
                    alpha = 2.0 - pow (beta, -(eta_c + 1.0));
                    if (rand <= (1.0 / alpha))
                    {
                        betaq = pow ((rand * alpha), (1.0 / (eta_c + 1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                    }
                    c1    = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    beta  = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                    alpha = 2.0 - pow (beta, -(eta_c + 1.0));
                    if (rand <= (1.0 / alpha))
                    {
                        betaq = pow ((rand * alpha), (1.0 / (eta_c + 1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                    }
                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
                    if (c1 < yl)
                        c1 = yl;
                    if (c2 < yl)
                        c2 = yl;
                    if (c1 > yu)
                        c1 = yu;
                    if (c2 > yu)
                        c2 = yu;
                    if (randomperc () <= 0.5)
                    {
                        child1->xreal[i] = c2;
                        child2->xreal[i] = c1;
                    }
                    else
                    {
                        child1->xreal[i] = c1;
                        child2->xreal[i] = c2;
                    }
                }
                else
                {
                    child1->xreal[i] = parent1->xreal[i];
                    child2->xreal[i] = parent2->xreal[i];
                }
            }
            else
            {
                child1->xreal[i] = parent1->xreal[i];
                child2->xreal[i] = parent2->xreal[i];
            }
        }
    }
    else
    {
        for (i = 0; i < number_variable; i++)
        {
            child1->xreal[i] = parent1->xreal[i];
            child2->xreal[i] = parent2->xreal[i];
        }
    }

    return;
}
