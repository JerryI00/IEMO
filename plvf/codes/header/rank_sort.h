/*
 * rank_sort.h:
 *  This file contains the sorting cmp and struct for sort with rank.
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

#ifndef SAMARITAN_RANK_SORT_H
#define SAMARITAN_RANK_SORT_H

typedef struct double_with_index
{
    int idx;
    double x;
} mysort;


int double_with_index_smaller_cmp (const void * a, const void * b);

int double_with_index_greater_cmp (const void * a, const void * b);


#endif //SAMARITAN_RANK_SORT_H
