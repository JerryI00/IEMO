/*
 * memory.h:
 *  This is the header file for memory operations.
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

# ifndef Samritan_MEMORY_H
# define Samritan_MEMORY_H

# include "../header/global.h"

void allocate_memory_ind (individual_real *ind);
void allocate_memory_pop (population_real *pop, int size);
void deallocate_memory_ind (individual_real *ind);
void deallocate_memory_pop (population_real *pop, int size);

# endif // Samritan_MEMORY_H
