/*
 * C2-DTLZ2.c
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
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

#include "../../header/problems.h"

void c2dtlz2 (individual_real *ind)
{
    int i, j, k, aux;
    double re, gx, fsum, sum1, sum2, sqr;
    double r;  // r determines the size, i.e., radius, of each feasible region
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->xreal;

    if (obj <= 3)   // recommended by C-NSGA-III paper
        r = 0.4;
    else
        r = 0.5;

    gx = 0.0;
    k  = number_variable - number_objective + 1;
    for(i = number_variable - k; i < number_variable; i++)
        gx += pow ((xreal[i] - 0.5), 2.0);

    for (i = 0; i < number_objective; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < number_objective; i++)
    {
        for (j = 0; j < number_objective - (i + 1); j++)
            obj[i] *= cos (PI * 0.5 * xreal[j]);
        if (i != 0)
        {
            aux     = number_objective - (i + 1);
            obj[i] *= sin (PI * 0.5 * xreal[aux]);
        }
    }
    fsum = 0.0;
    for(i = 0; i < number_objective; i++)
        fsum = fsum + obj[i] * obj[i];

    re   = (obj[0] - 1) * (obj[0] - 1) + fsum - obj[0] * obj[0] - r * r;
    sum2 = 0;
    sqr  = sqrt (number_objective);
    for (i = 0; i < number_objective; i++)
    {
        sum1 = (obj[i] - 1) * (obj[i] - 1) + fsum - obj[i] * obj[i] - r * r;
        if (sum1 < re)
            re = sum1;
        sum2 = sum2 + (obj[i] - 1.0 / sqr) * (obj[i] - 1.0 / sqr);
    }
    sum2 = sum2 - r * r;

    if (sum2 < re) re = sum2;
    if (re > 0) re = -re;
    else re = 0;
    ind->cv = re;
}
