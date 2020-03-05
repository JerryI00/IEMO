/*
 * analyse.c:
 *  This file contains the functions to record the results, including the population and metric values for analysis.
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

static double t[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

void track_evolution (void *ptr, int id, int end)
{
    int i, j;
    int read_ptr;

    char name[20];
    char id_char[10];

    char output_dir_level1[BUFSIZE_L];    // upper level directory
    char output_dir_level2[BUFSIZE_L];    // lower level directory
    char output_file[BUFSIZE_L];

    sprintf (id_char, "%d", id);
    // set the output directory
    sprintf (output_dir_level1, "./out/%s_M%d_D%d/%s/",
             problem_name,
             number_objective,
             number_variable,
             algorithm_name
    );
    sprintf (output_dir_level2, "./out/%s_M%d_D%d/%s/%d/",
             problem_name,
             number_objective,
             number_variable,
             algorithm_name,
             run_index
    );
    _mkdir (output_dir_level2);

    // first time output, init the parameter and output var and fun
    if (id == 1)
    {
        // set analyse list
        for (i = 0; i < BUFSIZE_S; i++)
            analyse_list[i] = 0;

        read_ptr = 0;
        while (1)
        {
            int name_c = 0;
            while (analyse_stream[read_ptr] != ' '
                   && analyse_stream[read_ptr] != '\t'
                   && analyse_stream[read_ptr] != '\n'
                   && analyse_stream[read_ptr] != 0)
            {
                name[name_c] = analyse_stream[read_ptr];
                name_c++;
                read_ptr++;
            }
            if (analyse_stream[read_ptr] == 0)
                name[name_c] = 0;
            name[name_c] = 0;

            if (!strcmp(name, "VAR"))
                analyse_list[VAR] = 1;
            else if (!strcmp(name, "FUN"))
                analyse_list[FUN] = 1;
            else if (!strcmp(name, "GD"))
                analyse_list[GD] = 1;
            else if (!strcmp(name, "IGD"))
                analyse_list[IGD] = 1;
            else if (!strcmp(name, "HV"))
                analyse_list[HV] = 1;
            else if (!strcmp(name, "PLOT"))
                analyse_list[PLOT] = 1;
            else if (!strcmp(name, "analyse:"))
                ;
            else
                print_error(1,2,"unknown setting for analyse ",name);

            if (analyse_stream[read_ptr] == 0)
                break;
            read_ptr++;
        }
    }

    if (runtime_output == 1 && (id % output_interval == 0 || id == 1 || end == 1))
    {

        if (analyse_list[VAR])
        {
            sprintf (output_file, "%smedium_VAR_%s.out", output_dir_level2, id_char);
            print_variable (output_file, ptr);
        }
        if (analyse_list[FUN])
        {
            sprintf (output_file, "%smedium_FUN_%s.out", output_dir_level2, id_char);
            print_objective (output_file, ptr);
        }
        if (analyse_list[GD])
        {
            record_gd (ptr, id);
        }
        if (analyse_list[IGD])
        {
            record_igd (ptr, id);
        }
        if (analyse_list[HV])
        {
            record_hv (ptr,id);
        }
        if (analyse_list[PLOT])
        {
            py_plot(ptr,id);
        }
    }

    if (end == 1)
    {
        if (analyse_list[VAR])
        {
            sprintf (output_file, "%sVAR%d.out", output_dir_level1, run_index);
            print_variable (output_file, ptr);
        }
        if (analyse_list[FUN])
        {
            sprintf (output_file, "%sFUN%d.out", output_dir_level1, run_index);
            print_objective (output_file, ptr);
        }
        if (analyse_list[GD])
        {
            if (runtime_output == 1)
            {
                sprintf (output_file, "%sGD_%d.txt", output_dir_level2, run_index);
                print_gd (output_file);
            }
        }
        if (analyse_list[IGD])
        {
            if (runtime_output == 1)
            {
                sprintf (output_file, "%sIGD_%d.txt", output_dir_level2, run_index);
                print_igd (output_file);
            }
        }
        if (analyse_list[HV])
        {
            if (runtime_output == 1)
            {
                sprintf (output_file, "%sHV_%d.txt", output_dir_level2, run_index);
                print_hv (output_file);
            }
        }
        if (analyse_list[PLOT])
        {
            py_plot(NULL,0);
            // for gnuplot
            sprintf (output_file, "%sFUN%d.out", output_dir_level1, run_index);
            gnu_plot(output_file, "FUN");
        }
    }
}
