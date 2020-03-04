/*
 * UF2.c
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

void uf2 (individual_real *ind)
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
        if (i % 2 == 0)
        {
            yj   = xreal[i - 1] - 0.3 * xreal[0] * (xreal[0] * cos (24.0 * PI * xreal[0] + 4.0 * i * PI / number_variable) + 2.0) * sin (6.0 * PI * xreal[0] + i * PI / number_variable);
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            yj   = xreal[i - 1] - 0.3 * xreal[0] * (xreal[0] * cos (24.0 * PI * xreal[0] + 4.0 * i * PI / number_variable) + 2.0) * cos (6.0 * PI * xreal[0] + i * PI / number_variable);
            sum1 += yj * yj;
            count1++;
        }
    }

    obj[0] = xreal[0] + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - sqrt (xreal[0]) + 2.0 * sum2 / (double)count2;
}
