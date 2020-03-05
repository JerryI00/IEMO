/*
 * initialization.h:
 *  This is the header file for initialization operations.
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

# ifndef Samaritan_INITIALIZATION_H
# define Samaritan_INITIALIZATION_H

# include "global.h"
# include "rand.h"

void calc_gvalue (individual_real *ind);
void initialize_individual_real (individual_real *ind);
void initialize_population_real (population_real *pop);
void read_population_real (population_real * pop, char * fileName);
void copy_ind (individual_real *ind1, individual_real *ind2);
void merge (population_real *pop1, population_real *pop2, population_real *pop3);

# endif // Samaritan_INITIALIZATION_H
