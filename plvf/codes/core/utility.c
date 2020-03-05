/*
 * utility.c:
 *  This file contains the functions to facilitate some common usages.
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

# include "../header/utility.h"
# include "../header/rank_sort.h"
#define NS 10
double getA(double arcs[NS][NS],int n)//按第一行展开计算|A|
{
   if(n==1)
   {
        return arcs[0][0];
    }
    double ans = 0.0;
    double temp[NS][NS];
    int i,j,k;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n-1;j++)
        {
            for(k=0;k<n-1;k++)
            {
                temp[j][k] = arcs[j+1][(k>=i)?k+1:k];

            }
        }
        double t = getA(temp,n-1);
        if(i%2==0)
        {
            ans += arcs[0][i]*t;
        }
        else
        {
            ans -=  arcs[0][i]*t;
        }
    }
    return ans;
}
void getAStart(double arcs[NS][NS],int n,double ans[NS][NS])//计算每一行每一列的每个元素所对应的余子式，组成A*
{
    if(n==1)
    {
        ans[0][0] = 1;
        return;
    }
    int i,j,k,t;
    double temp[NS][NS];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n-1;k++)
            {
                for(t=0;t<n-1;t++)
                {
                    temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
                }
            }


            ans[j][i]  =  getA(temp,n-1);
            if((i+j)%2 == 1)
            {
                ans[j][i] = - ans[j][i];
            }
        }
    }
}

/* Calculate the L2-norm of a vector */
double norm_vector (double *a)
{
    int i;
    double sum;

    sum = 0;
    for (i = 0; i < number_objective; i++)
        sum += a[i] * a[i];

    return sqrt (sum);
}

/* Calculate the Euclidean distance between two points */
double euclidian_distance (double *a, double *b, int dimension)
{
    int i;
    double distance;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
        distance += (a[i] - b[i]) * (a[i] - b[i]);

    return sqrt(distance);
}

double cos_angle(double *a,double *b,int dimension)
{
    int i;
    double dot_product=0.0;
    double length_a = 0.0,length_b = 0.0;
    double cos;

    for(i = 0; i < dimension; i++)
    {

        dot_product += a[i] * b[i];
    }
    for (i = 0; i < dimension; i++)
    {
        length_a +=  pow(a[i] , 2);
        length_b +=  pow(b[i] , 2);
    }
    length_a = sqrt(length_a);
    length_b = sqrt(length_b);
    cos = dot_product / (length_a * length_b);
    return cos;

}

double cos_angle_inverse(double *a,double *b,int dimension)
{
    int i;
    double dot_product=0.0;
    double length_a = 0.0,length_b = 0.0;
    double cos;

    for(i = 0; i < dimension; i++)
    {

        dot_product += a[i] * b[i];
    }
    for (i = 0; i < dimension; i++)
    {
        length_a +=  pow(a[i] , 2);
        length_b +=  pow(b[i] , 2);
    }
    length_a = sqrt(length_a);
    length_b = sqrt(length_b);
    cos = 1-dot_product / (length_a * length_b);
    return cos;

}
double normalised_euclidean_distance(const double *a, const double *b, const double* f_max, const double* f_min, int dimension) {

    int i;
    double distance;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
        distance += pow((a[i] - b[i]) / (f_max[i] - f_min[i]), 2);

    return sqrt(distance);

}


/* Build up multi-level directories */
void _mkdir (const char *dir)
{
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf (tmp, sizeof(tmp), "%s", dir);
    len = strlen (tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
    {
        if (*p == '/')
        {
            *p = 0;
            mkdir (tmp, S_IRWXU);
            *p = '/';
        }
    }
    mkdir (tmp, S_IRWXU);
}

/* Calculate the combinatorial number for n choose k */
int combination (int n, int k)
{
    int i;

    if (n < k)
        return -1;
    double ans = 1;
    for (i = k + 1; i <= n; i++)
    {
        ans = ans * i;
        ans = ans / (double) (i - k);
    }

    return (int) ans;
}

/* Shuffle the 'perm' array */
void random_permutation (int *perm, int size)
{
    int i, num, start;
    int *index, *flag;

    index = malloc (size * sizeof(int));
    flag  = malloc (size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        index[i] = i;
        flag[i]  = 1;
    }

    num = 0;
    while (num < size)
    {
        start = rnd (0, size - 1);
        while (1)
        {
            if (flag[start])
            {
                perm[num] = index[start];
                flag[start] = 0;
                num++;
                break;
            }
            if (start == (size - 1))
                start = 0;
            else
                start++;
        }
    }

    free (index);
    free (flag);
    return;
}

/* Update the current ideal point */
void update_ideal_point (individual_real *individual)
{
    int i;

    for (i = 0; i < number_objective; i++)
        if (individual->obj[i] < ideal_point[i])
            ideal_point[i] = individual->obj[i];

    return;
}


/* Update the current nadir point */
void update_nadir_point (individual_real *individual)
{
    int i;

    for (i = 0; i < number_objective; i++)
        if (individual->obj[i] > nadir_point[i])
            nadir_point[i] = individual->obj[i];

    return;
}

double weighted_euclidean_distance_ASF(const double *x, const double* y, const double* weights,
                                       const double* max, const double* min, int dimension) {

    double check_sum = 0;
    for(int i = 0; i < dimension; i++) {
        if(weights[i] > 1 || weights[i] < 0)
            print_error(1, 3, "Weight  ", i, " exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

        check_sum += weights[i];
    }

    if(check_sum > 1)
        print_error(1, 1, "Weights exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

    double sum = 0;
    for(int i = 0; i < dimension; i++) {
        double frac = (x[i] - y[i]) / (max[i] - min[i]);
        sum += weights[i] * frac * frac;
    }

    return sqrt(sum);

}

double tchebycheff_ASF(const double *x, const double* y, const double* weights, int dimensions) {

    double check_sum = 0;
    for(int i = 0; i < dimensions; i++) {
        if(weights[i] > 1 || weights[i] < 0)
            print_error(1, 3, "Weight  ", i, " exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

        check_sum += weights[i];
    }

    if(check_sum > 1)
        print_error(1, 1, "Weights exceeded bounds 0 to -1 at weighted_euclidean_distance_ASF");

    double curr_diff = 0, max_diff = -INF, diff_sum = 0;
    for(int i = 0; i < dimensions; i++) {

        curr_diff = weights[i] * (x[i] - y[i]);

        if(curr_diff > max_diff)
            max_diff = curr_diff;

        diff_sum += curr_diff;
    }

    return max_diff + (0.00001 * diff_sum);

}



void normalise_vector(double* vector, int length) {
    double z = 0;

    for(int i = 0; i < length; i++)
        z += vector[i];

    for (int i = 0; i < length; i++)
        vector[i] = vector[i] / z;
}
 void get_inverse_matrix(double arcs[NS][NS],int number)
 {

     double astar[NS][NS];
     int i,j;
     int n;
     n = number;
//     printf("n=\n");
//     while(scanf("%d",&n)!=EOF && n)
//     {
//         for(i=0;i<n;i++)
//         {
//             for(j=0;j<n;j++)
//             {
//                 scanf("%d",&arcs[i][j]);
//             }
//         }

         double a = getA(arcs,n);
         if(a==0)
         {
             printf("can not transform!\n");
         }
         else
         {
             getAStart(arcs,n,astar);
             for(i=0;i<n;i++)
             {
                 for(j=0;j<n;j++)
                 {
//                     printf("%.3lf ",astar[i][j]/a);
               arcs[i][j] = astar[i][j]/a;
                 }
//                 printf("\n");
             }
         }
//         printf("\n");


       return;

 }