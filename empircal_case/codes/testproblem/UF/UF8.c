/*
 * UF8.c
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

void uf8 (individual_real *ind)
{
    int i, count1, count2, count3;
    double sum1, sum2, sum3, yj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->xreal;

    sum1   = sum2   = sum3   = 0.0;
    count1 = count2 = count3 = 0;
    for (i = 3; i <= number_variable; i++)
    {
        yj = xreal[i - 1] - 2.0 * xreal[1] * sin (2.0 * PI * xreal[0] + i * PI / number_variable);
        if (i % 3 == 1)
        {
            sum1  += yj * yj;
            count1++;
        }
        else if (i % 3 == 2)
        {
            sum2  += yj * yj;
            count2++;
        }
        else
        {
            sum3  += yj * yj;
            count3++;
        }
    }

    obj[0] = cos (0.5 * PI * xreal[0]) * cos (0.5 * PI * xreal[1]) + 2.0 * sum1 / (double)count1;
    obj[1] = cos (0.5 * PI * xreal[0]) * sin (0.5 * PI * xreal[1]) + 2.0 * sum2 / (double)count2;
    obj[2] = sin (0.5 * PI * xreal[0]) + 2.0 * sum3 / (double)count3;
}
