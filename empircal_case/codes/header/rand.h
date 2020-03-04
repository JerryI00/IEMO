/*
 * rand.h:
 *  This is the header file for random number generation functions.
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

# ifndef Samaritan_RANDOM_H
# define Samaritan_RANDOM_H
/* Definition of random number generation routines */

# include "global.h"

extern double seed;
extern double oldrand[55];
extern int jrand;

void advance_random ();
void warmup_random (double seed);
void randomize ();
double randomperc ();
double read_randomperc();
int rnd (int low, int high);
double rndreal (double low, double high);

# endif // Samaritan_RANDOM_H
