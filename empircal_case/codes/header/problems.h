/*
 * problems.h:
 *  This is the header file for benchmark problems.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
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

# ifndef Samaritan_PROBLEMS_H
# define Samaritan_PROBLEMS_H

# include "../header/global.h"

void evaluate_population (population_real* pop);
void evaluate_individual (individual_real* ind);
void zdt1 (individual_real* ind);
void zdt2 (individual_real* ind);
void zdt3 (individual_real* ind);
void zdt4 (individual_real* ind);
void zdt6 (individual_real* ind);
void dtlz1 (individual_real* ind);
void dtlz2 (individual_real* ind);
void dtlz3 (individual_real* ind);
void dtlz4 (individual_real* ind);
void dtlz5 (individual_real* ind);
void dtlz6 (individual_real* ind);
void dtlz7 (individual_real* ind);
void uf1 (individual_real* ind);
void uf2 (individual_real* ind);
void uf3 (individual_real* ind);
void uf4 (individual_real* ind);
void uf5 (individual_real* ind);
void uf6 (individual_real* ind);
void uf7 (individual_real* ind);
void uf8 (individual_real* ind);
void uf9 (individual_real* ind);
void uf10 (individual_real* ind);
void wfg1 (individual_real* ind);
void wfg2 (individual_real* ind);
void wfg3 (individual_real* ind);
void wfg4 (individual_real* ind);
void wfg41 (individual_real* ind);
void wfg42 (individual_real* ind);
void wfg43 (individual_real* ind);
void wfg44 (individual_real* ind);
void wfg45 (individual_real* ind);
void wfg46 (individual_real* ind);
void wfg47 (individual_real* ind);
void wfg48 (individual_real* ind);
void wfg5 (individual_real* ind);
void wfg6 (individual_real* ind);
void wfg7 (individual_real* ind);
void wfg8 (individual_real* ind);
void wfg9 (individual_real* ind);

void c1dtlz1 (individual_real* ind);
void c1dtlz3 (individual_real* ind);
void c2dtlz2 (individual_real* ind);
void c3dtlz1 (individual_real* ind);
void c3dtlz4 (individual_real* ind);

# endif // Samaritan_PROBLEMS_H
