/*
 * tournament_selection.c:
 *  This file contains the functions to perform the mating selection according to a binary tournament procedure.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Ke Li, Renzhi Chen
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

# include "../../header/reproduction.h"

/* Routine for binary tournament */
individual_real* tournament (individual_real *ind1, individual_real *ind2)
{
    int flag;

    flag = check_dominance (ind1, ind2);
    if (flag == 1)
        return (ind1);
    else if (flag == -1)
        return (ind2);
    else
    {
        if (ind1->fitness > ind2->fitness)
            return(ind1);
        else if (ind2->fitness > ind1->fitness)
            return(ind2);
        else
        {
            if ((randomperc()) <= 0.5)
                return (ind1);
            else
                return (ind2);
        }
    }
}

individual_real* tournament_min (individual_real *ind1, individual_real *ind2) // For SPEA2.
{
    if (ind1->fitness < ind2->fitness)
        return(ind1);
    else if (ind2->fitness < ind1->fitness)
        return(ind2);
    else
    {
        if ((randomperc()) <= 0.5)
            return (ind1);
        else
            return (ind2);
    }
}
