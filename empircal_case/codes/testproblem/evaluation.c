/*
 * evaluation.c:
 *  This file contains the functions to perform function evaluations.
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

# include "../header/problems.h"
# include "../header/print.h"

void evaluate_population (population_real *pop)
{
    int i;

    for (i = 0; i < popsize; i++)
        evaluate_individual (&(pop->ind[i]));

    return;
}

void evaluate_individual (individual_real *ind)
{
    int flag;

    flag    = 0;
    (strcmp (problem_name, "ZDT1")  != 0)? :(zdt1 (ind), flag = 1);
    (strcmp (problem_name, "ZDT2")  != 0)? :(zdt2 (ind), flag = 1);
    (strcmp (problem_name, "ZDT3")  != 0)? :(zdt3 (ind), flag = 1);
    (strcmp (problem_name, "ZDT4")  != 0)? :(zdt4 (ind), flag = 1);
    (strcmp (problem_name, "ZDT6")  != 0)? :(zdt6 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ1") != 0)? :(dtlz1 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ2") != 0)? :(dtlz2 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ3") != 0)? :(dtlz3 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ4") != 0)? :(dtlz4 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ5") != 0)? :(dtlz5 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ6") != 0)? :(dtlz6 (ind), flag = 1);
    (strcmp (problem_name, "DTLZ7") != 0)? :(dtlz7 (ind), flag = 1);
    (strcmp (problem_name, "UF1") != 0)? :(uf1 (ind), flag = 1);
    (strcmp (problem_name, "UF2") != 0)? :(uf2 (ind), flag = 1);
    (strcmp (problem_name, "UF3") != 0)? :(uf3 (ind), flag = 1);
    (strcmp (problem_name, "UF4") != 0)? :(uf4 (ind), flag = 1);
    (strcmp (problem_name, "UF5") != 0)? :(uf5 (ind), flag = 1);
    (strcmp (problem_name, "UF6") != 0)? :(uf6 (ind), flag = 1);
    (strcmp (problem_name, "UF7") != 0)? :(uf7 (ind), flag = 1);
    (strcmp (problem_name, "UF8") != 0)? :(uf8 (ind), flag = 1);
    (strcmp (problem_name, "UF9") != 0)? :(uf9 (ind), flag = 1);
    (strcmp (problem_name, "UF10") != 0)? :(uf10 (ind), flag = 1);
    (strcmp (problem_name, "WFG1") != 0)? :(wfg1 (ind), flag = 1);
    (strcmp (problem_name, "WFG2") != 0)? :(wfg2 (ind), flag = 1);
    (strcmp (problem_name, "WFG3") != 0)? :(wfg3 (ind), flag = 1);
    (strcmp (problem_name, "WFG4") != 0)? :(wfg4 (ind), flag = 1);
    (strcmp (problem_name, "WFG41") != 0)? :(wfg41 (ind), flag = 1);
    (strcmp (problem_name, "WFG42") != 0)? :(wfg42 (ind), flag = 1);
    (strcmp (problem_name, "WFG43") != 0)? :(wfg43 (ind), flag = 1);
    (strcmp (problem_name, "WFG44") != 0)? :(wfg44 (ind), flag = 1);
    (strcmp (problem_name, "WFG45") != 0)? :(wfg45 (ind), flag = 1);
    (strcmp (problem_name, "WFG46") != 0)? :(wfg46 (ind), flag = 1);
    (strcmp (problem_name, "WFG47") != 0)? :(wfg47 (ind), flag = 1);
    (strcmp (problem_name, "WFG48") != 0)? :(wfg48 (ind), flag = 1);
    (strcmp (problem_name, "WFG5") != 0)? :(wfg5 (ind), flag = 1);
    (strcmp (problem_name, "WFG6") != 0)? :(wfg6 (ind), flag = 1);
    (strcmp (problem_name, "WFG7") != 0)? :(wfg7 (ind), flag = 1);
    (strcmp (problem_name, "WFG8") != 0)? :(wfg8 (ind), flag = 1);
    (strcmp (problem_name, "WFG9") != 0)? :(wfg9 (ind), flag = 1);
    (strcmp (problem_name, "C1DTLZ1")  != 0)? :(c1dtlz1 (ind), flag = 1);
    (strcmp (problem_name, "C1DTLZ3")  != 0)? :(c1dtlz3 (ind), flag = 1);
    (strcmp (problem_name, "C2DTLZ2")  != 0)? :(c2dtlz2 (ind), flag = 1);
    (strcmp (problem_name, "C3DTLZ1")  != 0)? :(c3dtlz1 (ind), flag = 1);
    (strcmp (problem_name, "C3DTLZ4")  != 0)? :(c3dtlz4 (ind), flag = 1);

    print_error (flag == 0, 2, "UNKNOWN test problem: ", problem_name);

    evaluation_count++;

    return;
}
