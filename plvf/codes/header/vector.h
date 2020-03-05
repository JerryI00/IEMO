/*
 * vector.h:
 *  This is the header file for the functions related to the data structure 'vector'.
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

#ifndef SAMARITAN_VECTOR_H
#define SAMARITAN_VECTOR_H

#include "global.h"
#include "print.h"

typedef struct int_vector
{
    int value;
    struct int_vector* next;
}i_vector;


typedef struct double_vector
{
    double value;
    struct double_vector* next;
}d_vector;

int int_vector_size (struct int_vector* head);
void int_vector_pushback (struct int_vector* head, int value);
int int_vector_pop (struct int_vector* head);
int int_vector_get (struct int_vector* head, int index);
void int_vector_set (struct int_vector* head, int index, int value);
void int_vector_free (struct int_vector* head);
void int_vector_print (struct int_vector* head);
void int_vector_remove(struct int_vector* head, int index);

int double_vector_size (struct double_vector* head);
void double_vector_pushback (struct double_vector* head, double value);
double double_vector_pop (struct double_vector* head);
double double_vector_get (struct double_vector* head, int index);
void double_vector_free (struct double_vector* head);
void double_vector_print (struct double_vector* head);

#endif //SAMARITAN_VECTOR_H
