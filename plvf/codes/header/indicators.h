/*
 * indicator.h:
 *  This is the header file for performance metrics.
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

# ifndef SAMARITAN_INDICATOR_H
# define SAMARITAN_INDICATOR_H

# include "global.h"
# include "analyse.h"
# include "print.h"
# include "vector.h"
# include "utility.h"

void record_gd (void *ptr, int id);
void print_gd (char *file_name);
double calculate_gd (void *ptr);

void record_igd (void *ptr, int id);
double calculate_igd (void *ptr);
void print_igd (char * file_name);

void record_hv (void *ptr, int id);
double calculate_hv (void *ptr);
void print_hv (char * file_name);

# endif // SAMARITAN_INDICATOR_H
