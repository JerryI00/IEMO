/*
 * plot.c:
 *  This file contains the functions to plot the population distribution (Python or GNU Plot).
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

# include "../header/analyse.h"

void *py_plot (void *pop,int gen)
{
    int i, j;

    if (pop == NULL && pythonplot != NULL)
    {
        fprintf (pythonplot, "%d \n", -1);
        return NULL;
    }
    population_real * ptr = pop;

    if (number_objective == 3)
    {
        if(pythonplot == NULL)
            pythonplot = popen ("python -W ignore ./script/plot3D.py", "w");
        print_error (pythonplot == NULL, 1, "cannot open python script");

        fprintf (pythonplot, "%d \n", popsize);
        fprintf (pythonplot, "%d \n", gen);
        for(i = 0; i < popsize; i++)
        {
            for(j = 0; j < number_objective; j++)
            {
                fprintf(pythonplot,"%lf \n", ptr->ind[i].obj[j]);
            }
        }
    }
    if (number_objective == 2)
    {
        if (pythonplot == NULL)
            pythonplot = popen ("python -W ignore ./script/plot2D.py", "w");

        print_error (pythonplot == NULL, 1, "cannot open python script");

        fprintf (pythonplot, "%d \n", popsize);
        fprintf (pythonplot, "%d \n", gen);

        for (i = 0; i < popsize; i++)
        {
            for (j = 0; j < number_objective; j++)
            {
                fprintf (pythonplot, "%lf \n", ptr->ind[i].obj[j]);
            }
        }
    }
    fflush (pythonplot);

    return NULL;
}

/* plot the population by GNU Plot */
void gnu_plot (char *filename, char *title)
{
    char command[BUFSIZE_L];
    FILE * gnuplot = NULL;

    gnuplot = popen ("gnuplot -persistent", "w");
    print_error (gnuplot == NULL, 1, "Error: gnuplot not found.");

    if (number_objective == 2)
    {
        sprintf (command, "set title \"%s\"\n", title);
        fprintf (gnuplot, "%s \n", command);
        sprintf (command, "plot \'%s\'\n", filename);
        fprintf (gnuplot, "%s \n", command);
        fflush (gnuplot);
    }
    else if (number_objective == 3)
    {
        sprintf (command,"set title \"%s\"\n",title);
        fprintf (gnuplot, "%s \n", command);
        fprintf (gnuplot, "set view 60, 130 \n");
        fprintf (gnuplot, "set ticslevel 0\n");
        fprintf (gnuplot, "set xrange [0:]\n");
        fprintf (gnuplot, "set yrange [0:]\n");
        fprintf (gnuplot, "set zrange [0:]\n");
        sprintf (command, "splot \'%s\'\n", filename);
        fprintf (gnuplot, "%s \n", command);
        fflush (gnuplot);
    }
    else
        print_error (1, 1, "Error: Unsupported number of objective in plot\n");
}