/*
 * UF1.c
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
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

# include "../../header/problems.h"

void uf1 (individual_real *ind)
{
    int i, count1, count2;
    double sum1, sum2, yj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->xreal;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    for (i = 2; i <= number_variable; i++)
    {
        yj = xreal[i - 1] - sin (6.0 * PI * xreal[0] + i * PI / number_variable);
        yj = yj * yj;
        if (i % 2 == 0)
        {
            sum2 += yj;
            count2++;
        }
        else
        {
            sum1 += yj;
            count1++;
        }
    }

    obj[0] = xreal[0] + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - sqrt (xreal[0]) + 2.0 * sum2 / (double)count2;
}
