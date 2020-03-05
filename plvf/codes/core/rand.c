/*
 * rand.c:
 *  This file contains the functions to generate random number (create next batch of 55 random numbers).
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

# include "../header/rand.h"

double seed;
double oldrand[55];
int jrand;
static FILE * rnd_input = NULL;

/* Get seed number for random and start it up */
void randomize ()
{
    int j1;
    for(j1 = 0; j1 <= 54; j1++)
    {
        oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random (seed);

    return;
}

/* Get randomize off and running */
void warmup_random (double seed)
{
    int j1, ii;
    double new_random, prev_random;

    oldrand[54] = seed;
    new_random = 0.000000001;
    prev_random = seed;
    for(j1 = 1; j1 <= 54; j1++)
    {
        ii = (21 * j1) % 54;
        oldrand[ii] = new_random;
        new_random  = prev_random-new_random;
        if(new_random < 0.0)
        {
            new_random += 1.0;
        }
        prev_random = oldrand[ii];
    }
    advance_random ();
    advance_random ();
    advance_random ();
    jrand = 0;

    return;
}

/* Create next batch of 55 random numbers */
void advance_random ()
{
    int j1;
    double new_random;
    for(j1 = 0; j1 < 24; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1 + 31];
        if(new_random < 0.0)
        {
            new_random = new_random+1.0;
        }
        oldrand[j1] = new_random;
    }
    for(j1 = 24; j1 < 55; j1++)
    {
        new_random = oldrand[j1] - oldrand[j1 - 24];
        if(new_random < 0.0)
        {
            new_random = new_random + 1.0;
        }
        oldrand[j1] = new_random;
    }
}

/* Fetch a single random number between 0.0 and 1.0 */
double randomperc ()
{
    jrand++;
    if(jrand >= 55)
    {
        jrand = 1;
        advance_random();
    }
    return((double)oldrand[jrand]);
}

double read_randomperc ()
{
    if(rnd_input==NULL)
        rnd_input = fopen("rnd.out","r");
    double temp;
    if(fscanf(rnd_input,"%lf",&temp)==EOF){
        fseek(rnd_input,0,SEEK_SET);
        fscanf(rnd_input,"%lf",&temp);
    };

   // printf("r: %lf\n",temp);
    return temp;
}

/* Fetch a single random integer between low and high including the bounds */
int rnd (int low, int high)
{
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        //res = low + (int)(randomperc() * (high - low));
        res = low + (unsigned int)(randomperc() * (high - low + 1));
        if (res > high)
        {
            res = high;
        }
        else if (res < low)
        {
            res = low;
        }
    }
    return (res);
}

/* Fetch a single random real number between low and high including the bounds */
double rndreal (double low, double high)
{
    return (low + (high - low) * randomperc());
}
