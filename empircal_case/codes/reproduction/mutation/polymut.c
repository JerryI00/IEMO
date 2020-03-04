/*
 * polymut.c:
 *  This file contains the functions to perform polynomial mutation operation.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Ke Li, Renzhi Chen
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
# include "../../header/reproduction.h"

/* Routine for real polynomial mutation of an individual */
void polymut_ind (individual_real *ind)
{
    int i;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    for (i = 0; i < number_variable; i++)
    {
        if (randomperc() <= pmut_real)
        {
            y       = ind->xreal[i];
            yl      = variable_lowerbound[i];
            yu      = variable_upperbound[i];
            delta1  = (y - yl) / (yu - yl);
            delta2  = (yu - y) / (yu - yl);
            rnd     = randomperc();
            mut_pow = 1.0 / (eta_m + 1.0);
            if (rnd <= 0.5)
            {
                xy     = 1.0 - delta1;
                val    = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (eta_m + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            }
            else
            {
                xy     = 1.0 - delta2;
                val    = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (eta_m + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }
            y = y + deltaq * (yu - yl);
            if (y < yl)
                y = yl;
            if (y > yu)
                y = yu;
            ind->xreal[i] = y;

        }
    }

    return;
}
