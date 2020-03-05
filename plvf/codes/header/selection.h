/*
 * selection.h:
 *  This is the header file for the environmental selection operations.
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

# ifndef Samaritan_SELECTION_H
# define Samaritan_SELECTION_H

# include "../header/global.h"
# include "../header/utility.h"
# include "../header/print.h"
# include "../header/rank_sort.h"
# include "../header/population.h"
# include "../header/dominance.h"
# include "../externals/IWFG/iwfg.h"

/* IBEA */
typedef enum {HYPERVOLUME, EPSILON} indicator_type ;
extern indicator_type indicator;
double asf_pbea (double *a, double *reference_point, double *weights);
double cal_eps_indicator (individual_real *ind1, individual_real *ind2, population_real *mixed_pop);
void calcFitness (population_real *mixed_pop, double **Indicator_value, int size);
void cal_fitnesses(void *ptr, double **fitcomp, int size);
void environmental_selection (population_real *mixed_pop, population_real *new_pop, int *flag, double **Indicator_value, int size);
void ibea_selection (population_real *mixed_pop, population_real *new_pop, int *flag, double **Indicator_value);
void pbea_calcFitness(population_real *mixed_pop, double **Indicator_value, int size, double* weights, double* reference_point, double specificity);
void pbea_selection(population_real *mixed_pop, population_real *new_pop, int *flag, double **Indicator_value, double* weights, double* reference_point, double specificity);

/* NSGA-II */
void fill_nondominated_sort (population_real *new_pop, population_real *mixed_pop);
void fill_r_nondominated_sort (population_real *new_pop, population_real *mixed_pop, double sigma);
void fill_R_nondominated_sort (population_real *new_pop, population_real *mixed_pop, double *reference_point, double *weights, double epsilon);
void crowding_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite);
void r_crowding_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite, double* reference_point, double* weights, double epsilon);

void assign_crowding_distance (population_real *pop, int *dist, int **obj_array, int front_size);
void assign_crowding_distance_list (population_real *pop, list *lst, int front_size);
void assign_crowding_distance_indices (population_real *pop, int c1, int c2);

void quicksort_front_obj(population_real *pop, int objcount, int obj_array[], int obj_array_size);
void q_sort_front_obj(population_real *pop, int objcount, int obj_array[], int left, int right);
void quicksort_dist(population_real *pop, int *dist, int front_size);
void q_sort_dist(population_real *pop, int *dist, int left, int right);

/* MOEA/D */
double fitnessFunction (individual_real *individual, double *lambda);
double Rbf_fitnessFunction (individual_real *individual, double *lambda);
void comp_utility (population_real *pop, population_real *saved_values);
void moead_free ();
void nums_weight();
void uniform_random_weights(double **lambda, int num_weights, int dimensions, const double* lower_bounds, const double* upper_bounds);
void uniform_random_weights_const(double **lambda, int num_weights, int dimensions, double lower_bound, double upper_bound);
void initialize_uniform_weight ();
void initialize_miu_weight (int miu);
void initialize_layers_weight ();
void read_uniform_weight(char *file);
void initialize_neighborhood ();
void set_weight (double *weight, double unit, double sum, int dim);
void tour_selection_subproblem (int depth);
void update_subproblem (population_real *pop, individual_real *individual, int subProblem_id, int neighborType);
void update_subproblem_constraint (population_real *pop, individual_real *individual, int subProblem_id, int neighborType);
void update_neighborhood (double *weight,int miu,int number,double eta);
void rbf_train(double **a,double *rbf,double *weights_rbf,double digma);
double rbf_test(individual_real *individual,double *weights_rbf,double **a,double sigma);
void rbf_first_train(double **a,double *rbf,double *weights_rbf,double sigma,int miu);
/* MOEA/D-STM */
double calculateDistance2 (individual_real* individual, double* lambda);
double norm_vector (double* z);
int prefers (int x, int y, struct double_with_index* womanPref, int size);
void stableMatching (int *statusMan,int * statusWoman , int * next, struct double_with_index** man_pref, struct double_with_index** woman_pref, int menSize, int womenSize);
void stm_selection (population_real *parent_pop, population_real *mixed_pop);
void stm_dra_selection (population_real *parent_pop, population_real *mixed_pop, int size);
void rbf_set_weight (double *weight, double unit, double sum, int dim);
/* SMS-EMOA */
void fill_constraint_hv_sort (FILECONTENTS *f,population_real* new_pop, population_real* mixed_pop,int size);
void fill_hv_sort (FILECONTENTS *f,population_real* new_pop, population_real* mixed_pop,int size);
void hv_fill (FILECONTENTS *f,population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite);

/* hype */
void fill_hype_sort (population_real* new_pop, population_real* mixed_pop, int size);
void hype_fill (population_real *mixed_pop, population_real *new_pop, int count, int front_size, list *elite);

/* SPEA2 */
void fitness_spea2 (population_real *pop, int total_size, int k_min, int *dominated_Num, int *bedominated_Num, int **dominated_Matrix, int *R_i, double **distance_Matrix, double* D_i, double *kth_distance);
void selection_spea2 (population_real *mixed_pop,int total_size,population_real *archive,int archive_size, individual_real *temp_ind, population_real *temp_pop, double **distance_Matrix);

# endif // Samaritan_SELECTION_H
