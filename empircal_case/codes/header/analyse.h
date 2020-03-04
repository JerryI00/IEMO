/*
 * dominance.h:
 *  This is the header file for dominance check.
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

# ifndef SAMARITAN_ANALYSE_H
# define SAMARITAN_ANALYSE_H

# include "global.h"
# include "population.h"
# include "print.h"
# include "indicators.h"

void track_evolution (void *ptr, int id,int end);
void gnu_plot (char *file, char * title);
void* py_plot(void *ptr,int gen);

# endif // SAMARITAN_ANALYSE_H
