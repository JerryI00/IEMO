/*
 * DTLZ7.c
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

void dtlz7 (individual_real *ind)
{
    int i, j, k;
    int aux;
    double h, gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->xreal;

    gx = 0.0;
    k  = number_variable - number_objective + 1;
    for(i = number_variable - k; i < number_variable; i++)
        gx += xreal[i];
    gx = 1.0 + (9.0 * gx) / k;

    for (i = 0; i < number_objective; i++)
        obj[i] = xreal[i];

    h = 0.0;
    for (i = 0; i < number_objective - 1; i++)
        h += (obj[i] / (1.0 + gx)) * (1.0 + sin (3.0 * PI * obj[i]));
    h = number_objective - h;

    obj[number_objective - 1] = (1 + gx) * h;
}
