/*
 * utility.h:
 *  This is the header file for the utility functions.
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

#ifndef SAMARITAN_UTILITY_H
#define SAMARITAN_UTILITY_H

# include "../header/global.h"
# include "../header/rand.h"

void _mkdir (const char *dir);
double euclidian_distance (double *a, double *b, int dimension);
double normalised_euclidean_distance(const double *a, const double *b, const double* f_max, const double* f_min, int dimension);
int combination(int n, int k);
void random_permutation(int* perm, int size);
void update_ideal_point (individual_real *individual);
void update_nadir_point (individual_real *individual);
double weighted_euclidean_distance_ASF(const double *x, const double* y, const double* weights,
                                       const double* max, const double* min, int dimension);
double tchebycheff_ASF(const double *x, const double* y, const double* weights, int dimensions);



void normalise_vector(double* vector, int length);
#endif //SAMARITAN_UTILITY_H
