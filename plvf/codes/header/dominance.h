/*
 * dominance.h:
 *  This is the header file for dominance check.
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

# ifndef Samaritan_DOMINANCE_H
# define Samaritan_DOMINANCE_H

# include "global.h"

int check_dominance (individual_real *a, individual_real *b);
int check_g_dominance (individual_real *a);
int check_r_dominance (individual_real *a, individual_real *b, population_real *mixed_pop, double sigma);
int dist_measure (population_real *mixed_pop, individual_real *a, individual_real *b, double sigma);

list** nondominated_sort_idxs (population_real* pop, int popsize);

# endif // Samaritan_DOMINANCE_H
