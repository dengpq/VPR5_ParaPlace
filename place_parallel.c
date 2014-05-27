/*#include <stdlib.h> */
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h> //gettimeofday

#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "mst.h"
#include "place_parallel.h"
#include "read_place.h"
#include "draw.h"
#include "place_and_route.h"
#include "net_delay.h"
#include "path_delay_parallel.h"
#include "path_delay2_parallel.h"
#include "timing_place_lookup.h"
#include "timing_place.h"
#include "place_stats.h"

#define _XOPEN_SOURCE 600
#define T_CONSTANT_GENERATOR -1000  /* Essentially -ve infinity */
#define PROB 10
#define TIMING_UPDATE 5
#define PARITION_UPDATE 80  /* why? */


/* Partitioning parameter.
 * The grid will be partitioned into x_partition rows *
 * and y_partition columns.                *
 * Feel free to change the numbers below.  *
 * For 8-threads, I'd like to set 4x2 regions. *
 * FIXME, when I set regions, the number of horizontal regions
 * should >= the number of vertical ones. That is *
 * x_partition[NUM_OF_THREADS] >= y_partition[NUM_OF_THREADS]. */
int x_partition[64] = {
    1, 1, 1, 2, 1, /*  1-5   */
    3, 0, 4, 3, 0, /*  6-10  */
    0, 0, 0, 0, 3, /* 11-15  */
    4, 0, 3, 0, 4, /* 16-20  */
    0, 2, 0, 4, 5, /* 21-25  */
    0, 0, 0, 0, 0, /* 26-30  */
    0, 0, 0, 0, 0, /* 31-35  */
    6, 0, 0, 0, 0, /* 36-40  */
    0, 0, 0, 0, 0, /* 41-45  */
    0, 0, 0, 7, 0, /* 46-50  */
    0, 0, 0, 0, 0, /* 51-55  */
    0, 0, 0, 0, 6, /* 56-60  */
    0, 0, 7, 8   /* 61-64 */
};

int y_partition[64] = {
    1, 2, 3, 2, 5, /*  1-5  */
    2, 0, 2, 3, 0, /*  6-10 */
    0, 0, 0, 0, 5, /* 11-15 */
    4, 0, 6, 0, 5, /* 16-20 */
    0, 11, 0, 6, 5,/* 21-25 */
    0, 0, 0, 0, 0, /* 26-30 */
    0, 0, 0, 0, 0, /* 31-35 */
    6, 0, 0, 0, 0, /* 36-40 */
    0, 0, 0, 0, 0, /* 41-45 */
    0, 0, 0, 7, 0, /* 46-50 */
    0, 0, 0, 0, 0, /* 51-55 */
    0, 0, 0, 0, 10,/* 56-60 */
    0, 0, 9, 8     /* 61-64 */
};

/*****************************************************************************
 *       The following function were designed for Parallel Placememt         *
 ****************************************************************************/
static void* try_place_parallel(void* args);

static void  try_place_a_subregion(const int  kthread_id,
                                   const int  krow,
                                   const int  kcol,
                                   const int  prob,
                                   int*  move_counter,
                                   thread_local_common_paras_t*  common_paras_ptr,
                                   thread_local_data_for_swap_t* swap_data_ptr);

/* try_place_parallel directed-called functions */
static void tp_initialize(pthread_data_t*  input_args,
                          int*  max_pins_per_fb,
                          thread_local_common_paras_t*  common_paras_ptr);

static void tp_alloc_mem(const int max_pins_per_fb,
                         const placer_opts_t placer_opts,
                         thread_local_data_for_swap_t* swap_data_ptr);


static void tp_init_local(const int* region_x_boundary,
                          thread_local_data_for_swap_t*  swap_data_ptr);

static void tp_init_localvert_grid();

static void tp_timing_update_full(int thread_id,
                                  double*** net_delay,
                                  pthread_data_t* input_args,
                                  double*  max_delay);

static void tp_compute_net_slack(int thread_id,
                                 double*** net_slack,
                                 double*   timing_cost,
                                 double*   delay_cost,
                                 double    crit_exponent,
                                 double    max_delay,
                                 pthread_data_t* input_args);

static void tp_timing_calc(int thread_id,
                           double* timing_cost,
                           double* delay_cost,
                           pthread_data_t* input_args);

static void tp_iter_data_update(const int kiter,
                                pthread_data_t* input_args,
                                thread_local_common_paras_t*  common_paras_ptr,
                                thread_local_data_for_swap_t* swap_data_ptr);

static void balance_two_consecutive_threads_edge(int thread_id);

static void balance_two_consecutive_threads_sinks(int thread_id);

static void comp_delta_td_cost_parallel(int from_block,
                                        int to_block,
                                        int num_of_pins,
                                        double* delta_timing,
                                        double* delta_delay,
                                        local_block_t* local_block,
                                        double** local_temp_point_to_point_timing_cost,
                                        double** local_temp_point_to_point_delay_cost);

static void tp_local_data_update(const int* region_x_boundary,
                                 const int* region_y_boundary,
                                 grid_tile_t** local_grid, local_block_t* local_block,
                                 const int krow, const int kcol);

/* try swap a pair of blocks in each Extend-SubRegion by each thread parallely */
static int  try_swap_parallel(const double t,
                              const int  kthread_id,
                              const int  x_from,
                              const int  y_from,
                              const int  z_from,
                              const int  from_block,
                              const int  place_cost_type,
                              const place_algorithm_t place_algorithm,
                              const double timing_tradeoff,
                              double*  total_cost,
                              double*  bb_cost,
                              double*  timing_cost,
                              double*  delay_cost,
                              double  inverse_prev_bb_cost,
                              double  inverse_prev_timing_cost,
                              double  range_limit,
                              int xMin, int xMax,
                              int yMin, int yMax,
                              thread_local_data_for_swap_t* swap_data_ptr);

static void get_bb_from_scratch_parallel(int inet,
                                         bbox_t* coords,
                                         bbox_t* num_on_edges,
                                         local_block_t* local_block);

static void get_non_updateable_bb_parallel(int inet,
                                           bbox_t* bb_coord_new,
                                           local_block_t* local_block);

static void update_bb_parallel(int inet,
                               bbox_t* local_bb_coord,
                               bbox_t* local_bb_edge,
                               bbox_t* bb_coord_new,
                               bbox_t* bb_edge_new,
                               int xold,
                               int yold,
                               int xnew,
                               int ynew,
                               local_block_t* local_block);

static boolean find_to_block_parallel(int x_from,
                                      int y_from,
                                      block_type_ptr type,
                                      int* x_to,
                                      int* y_to,
                                      int thread_id,
                                      int xmin, int xmax,
                                      int ymin, int ymax);

static int find_affected_nets_parallel(int* nets_to_update,
                                       int* net_block_moved,
                                       int from_block,
                                       int to_block,
                                       int num_of_pins,
                                       double* local_temp_net_cost);

static int assess_swap_parallel(double delta_c,
                                double t,
                                int local_seed);

/* 4 functions for data broadcast*/
static void update_from_local_to_global(local_block_t* local_block,
                                        grid_tile_t** local_grid,
                                        int x_start, int x_end,
                                        int y_start, int y_end);

static void update_from_global_to_local_hori(local_block_t* local_block,
                                             grid_tile_t** local_grid,
                                             int x_start, int x_end,
                                             int y_start, int y_end);

static void update_from_global_to_local_vert(local_block_t* local_block,
                                             grid_tile_t** local_grid,
                                             int x_start, int x_end,
                                             int y_start, int y_end);

static void update_from_global_to_local_grid_only(grid_tile_t** local_grid,
                                                  const int x_start,
                                                  const int x_end,
                                                  const int y_start,
                                                  const int y_end);
/* calculates the time difference*/
static double my_difftime2(struct timeval* start,
                           struct timeval* end);

static void print_grid(void);
/******************   Parallel Functions  Ending   ********************/


/************** Types and defines local to place.c ***************************/
#define SMALL_NET 4   /* Cut off for incremental bounding box updates. */
/* 4 is fastest -- I checked.                    */


/* For comp_cost.  NORMAL means use the method that generates updateable  *
 * bounding boxes for speed.  CHECK means compute all bounding boxes from *
 * scratch using a very simple routine to allow checks of the other       *
 * costs.                                                                 */
enum cost_methods {
    NORMAL,
    CHECK
};

#define FROM 0          /* What block connected to a net has moved? */
#define TO 1
#define FROM_AND_TO 2

#define ERROR_TOL .0025
#define MAX_MOVES_BEFORE_RECOMPUTE 50000

/********************** Variables local to place.c ***************************/
/* grid_tile_t** localvert_grid[num_grid_columns+2][num_grid_rows+2] */
grid_tile_t** localvert_grid = NULL;


/* [0..num_nets-1]  0 if net never connects to the same block more than  *
 *  once, otherwise it gives the number of duplicate connections.        */
static int* duplicate_pins;


/* [0..num_nets-1][0..num_unique_blocks-1]  Contains a list of blocks with *
 * no duplicated blocks for ONLY those nets that had duplicates.           */
static int** unique_pin_list;


/* Cost of a net, and a temporary cost of a net used during move assessment. */
static double* net_cost = NULL, *temp_net_cost = NULL;   /* [0..num_nets-1] */

/* [0..num_nets-1][1..num_pins-1]. What is the value of the timing   */
/* driven portion of the cost function. These arrays will be set to  */
/* (criticality * Tdel) for each point to point connection. */
static double** point_to_point_timing_cost = NULL;
static double** temp_point_to_point_timing_cost = NULL;



/* [0..num_nets-1][1..num_pins-1]. What is the value of the Tdel */
/* for each connection in the circuit */
static double** point_to_point_delay_cost = NULL;
static double** temp_point_to_point_delay_cost = NULL;


/* [0..num_blocks-1][0..pins_per_clb-1]. Indicates which pin on the net */
/* this block corresponds to, this is only required during timing-driven */
/* placement. It is used to allow us to update individual connections on */
/* each net */
static int** net_pin_index = NULL;


/* [0..num_nets-1].  Store the bounding box coordinates and the number of *
 * blocks on each of a net's bounding box (to allow efficient updates),   *
 * respectively.                                                          */
static bbox_t* bb_coords = NULL;
static bbox_t* bb_num_on_edges = NULL;

/* Stores the maximum and expected occupancies, plus the cost, of each   *
 * region in the placement.  Used only by the NONLINEAR_CONG cost        *
 * function.  [0..num_region-1][0..num_region-1].  Place_region_x and    *
 * y give the situation for the x and y directed channels, respectively. */
static place_region_t** place_region_x, ** place_region_y;

/* Used only with nonlinear congestion.  [0..num_regions].            */
static double* place_region_bounds_x, *place_region_bounds_y;

/* The arrays below are used to precompute the inverse of the average   *
 * number of tracks per channel between [subhigh] and [sublow].  Access *
 * them as chan?_place_cost_fac[subhigh][sublow].  They are used to     *
 * speed up the computation of the cost function that takes the length  *
 * of the net bounding box in each dimension, divided by the average    *
 * number of tracks in that direction; for other cost functions they    *
 * will never be used.                                                  */
static double** chanx_place_cost_fac;
static double** chany_place_cost_fac;


/* Expected crossing counts for nets with different #'s of pins.  From *
 * ICCAD 94 pp. 690 - 695 (with linear interpolation applied by me).   */
static const double cross_count[50] = {  /* [0..49] */
    1.0, 1.0, 1.0, 1.0828, 1.1536, 1.2206, 1.2823, 1.3385, 1.3991, 1.4493,
    1.4974, 1.5455, 1.5937, 1.6418, 1.6899, 1.7304, 1.7709, 1.8114, 1.8519,
    1.8924,
    1.9288, 1.9652, 2.0015, 2.0379, 2.0743, 2.1061, 2.1379, 2.1698, 2.2016,
    2.2334,
    2.2646, 2.2958, 2.3271, 2.3583, 2.3895, 2.4187, 2.4479, 2.4772, 2.5064,
    2.5356,
    2.5610, 2.5864, 2.6117, 2.6371, 2.6625, 2.6887, 2.7148, 2.7410, 2.7671,
    2.7933
};

/********************* Static subroutines local to place.c *******************/
static void alloc_and_load_unique_pin_list(void);

static void free_unique_pin_list(void);

static void alloc_place_regions(int num_regions);

static void load_place_regions(int num_regions);

static void free_place_regions(int num_regions);

static void alloc_and_load_placement_structs(int place_cost_type,
                                             int num_regions,
                                             double place_cost_exp,
                                             double** *old_region_occ_x,
                                             double** *old_region_occ_y,
                                             placer_opts_t
                                             placer_opts);

static void free_placement_structs(int place_cost_type,
                                   int num_regions,
                                   double** old_region_occ_x,
                                   double** old_region_occ_y,
                                   placer_opts_t placer_opts);

static void alloc_and_load_for_fast_cost_update(double place_cost_exp);

static void initial_placement(pad_loc_t pad_loc_type,
                              char* pad_loc_file);

static double comp_bb_cost_parallel(int start_net,
                                    int end_net);

static double comp_bb_cost(int method,
                           int place_cost_type,
                           int num_regions);

/* using double instead of double to alleviate floating point round-off error */
static int try_swap(double t,
                    double* cost,
                    double* bb_cost,
                    double* timing_cost,
                    double range_limit,
                    int place_cost_type,
                    double** old_region_occ_x,
                    double** old_region_occ_y,
                    int num_regions,
                    boolean fixed_pins,
                    place_algorithm_t place_algorithm,
                    double timing_tradeoff,
                    double inverse_prev_bb_cost,
                    double inverse_prev_timing_cost,
                    double* delay_cost,
                    int* x_lookup);

static double check_place(double bb_cost,
                         double timing_cost,
                         int place_cost_type,
                         int num_regions,
                         place_algorithm_t place_algorithm,
                         double delay_cost);

static double starting_t(double* cost_ptr,
                        double* bb_cost_ptr,
                        double* timing_cost_ptr,
                        int place_cost_type,
                        double** old_region_occ_x,
                        double** old_region_occ_y,
                        int num_regions,
                        boolean fixed_pins,
                        annealing_sched_t annealing_sched,
                        int max_moves,
                        double range_limit,
                        place_algorithm_t place_algorithm,
                        double timing_tradeoff,
                        double inverse_prev_bb_cost,
                        double inverse_prev_timing_cost,
                        double* delay_cost_ptr);

static void update_t_parallel(double* t,
                              int* inner_iter_num,
                              double range_limit,
                              double success_rat,
                              annealing_sched_t annealing_sched);

static void update_rlim(double* range_limit,
                        double success_rat);

static int exit_crit(double t,
                     double cost,
                     annealing_sched_t annealing_sched);

static int count_connections(void);

static void compute_net_pin_index_values(void);

static double get_std_dev(int n,
                          double sum_x_squared,
                          double av_x);

static void free_fast_cost_update_structs(void);

static double comp_td_point_to_point_delay(int inet,
                                           int ipin);

static double comp_td_point_to_point_delay_parallel(int inet,
                                                    int ipin,
                                                    local_block_t* local_block);

static void update_td_cost(int from_block,
                           int to_block,
                           int num_of_pins);

static void comp_delta_td_cost(int from_block,
                               int to_block,
                               int num_of_pins,
                               double* delta_timing,
                               double* delta_delay);

static unsigned long compute_td_costs_parallel_without_update_crit(double* timing_cost,
                                                                   double* connection_delay_sum,
                                                                   int start_nets,
                                                                   int finish_nets);

static unsigned long compute_td_costs_parallel_with_update_crit(double* timing_cost,
                                                                double* connection_delay_sum,
                                                                double** net_slack,
                                                                double max_delay,
                                                                double crit_exponent,
                                                                int start,
                                                                int finish);

static void comp_td_costs(double* timing_cost,
                          double* connection_delay_sum);

static int assess_swap(double delta_c,
                       double t);

static boolean find_to(int x_from,
                       int y_from,
                       block_type_ptr type,
                       double range_limit,
                       int* x_lookup,
                       int* x_to,
                       int* y_to);

static void get_non_updateable_bb(int inet,
                                  bbox_t* bb_coord_new);

static void update_bb(int inet,
                      bbox_t* bb_coord_new,
                      bbox_t* bb_edge_new,
                      int xold,
                      int yold,
                      int xnew,
                      int ynew);

static int find_affected_nets(int* nets_to_update,
                              int* net_block_moved,
                              int from_block,
                              int to_block,
                              int num_of_pins);

static double get_net_cost(int inet,
                          bbox_t* bb_ptr);

static double nonlinear_cong_cost(int num_regions);

static void update_region_occ(int inet,
                              bbox_t* coords,
                              int add_or_sub,
                              int num_regions);

static void save_region_occ(double** old_region_occ_x,
                            double** old_region_occ_y,
                            int num_regions);

static void restore_region_occ(double** old_region_occ_x,
                               double** old_region_occ_y,
                               int num_regions);

static void get_bb_from_scratch(int inet,
                                bbox_t* coords,
                                bbox_t* num_on_edges);

static double get_net_wirelength_estimate(int inet,
                                          bbox_t* bbptr);
/*===================   ORIGIN Functions in old vpr4.3    ===============*/

static void tp_data_print_to_screen(int thread_id,
                             double* timing_cost,
                             double* delay_cost,
                             double* crit_exponent,
                             double max_delay,
                             pthread_data_t* input_args,
                             placer_opts_t placer_opts,
                             double* cost,
                             double* bb_cost,
                             int* success_sum,
                             double* sum_of_squares,
                             int* move_counter,
                             int* tot_iter, double* success_rat,
                             double* av_cost, double* av_bb_cost,
                             double* av_timing_cost,
                             double* av_delay_cost, double* std_dev,
                             double* range_limit,
                             double* final_rlim,
                             double* inverse_delta_rlim,
                             double* t,
                             double* oldt, int* inner_iter_num,
                             double place_delay_value);

/***********************************************************************/
/* RESEARCH TODO: Bounding Box and range_limit need to be redone for    *
 * heterogeneous to prevent a QoR penalty */
/* Does almost all the work of placing a circuit. Width_fac gives the   *
 * width of the widest channel. Place_cost_exp says what exponent the   *
 * width should be taken to when calculating costs. This allows a       *
 * greater bias for anisotropic architectures. Place_cost_type          *
 * determines which cost function is used. num_regions is used only     *
 * the place_cost_type is NONLINEAR_CONG.                               */
void try_place_use_multi_threads(placer_opts_t      placer_opts,
                                 annealing_sched_t  annealing_sched,
                                 chan_width_distr_t chan_width_dist,
                                 router_opts_t      router_opts,
                                 detail_routing_arch_t det_routing_arch,
                                 segment_info_t*    segment_inf,
                                 timing_info_t      timing_inf,
                                 subblock_data_t*   subblock_data_ptr,
                                 t_mst_edge***      mst,
                                 operation_types_t  operation)
{
    /* Allocated here because it goes into timing critical code where each memory allocation is expensive */
    /* Used to quickly determine valid swap columns */
    printf("Before placement, Let's print the grid location...\n");
    print_grid();

    int* x_lookup = (int*)my_malloc(num_grid_columns * sizeof(int));

    /*used to free net_delay if it is re-assigned */
    double** remember_net_delay_original_ptr = NULL;
    double** net_slack, **net_delay;
    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
        placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE ||
        placer_opts.enable_timing_computations) {
        /*do this before the initial placement to avoid messing up the initial placement */
        alloc_lookups_and_criticalities(chan_width_dist,
                                        router_opts,
                                        det_routing_arch,
                                        segment_inf,
                                        timing_inf,
                                        *subblock_data_ptr,
                                        &net_delay,
                                        &net_slack);
        remember_net_delay_original_ptr = net_delay;
    }

    int width_fac = placer_opts.place_chan_width;
    boolean fixed_pins; /* Can pads move or not? */
    if (placer_opts.pad_loc_type == FREE) {
        fixed_pins = FALSE;
    } else {
        fixed_pins = TRUE;
    }

    init_chan(width_fac, chan_width_dist);

    double** old_region_occ_x, **old_region_occ_y;
    alloc_and_load_placement_structs(placer_opts.place_cost_type,
                                     placer_opts.num_regions,
                                     placer_opts.place_cost_exp,
                                     &old_region_occ_x,
                                     &old_region_occ_y,
                                     placer_opts);

    initial_placement(placer_opts.pad_loc_type, placer_opts.pad_loc_file);
    init_draw_coords((double)width_fac);

    /* Storing the number of pins on each type of block makes the swap routine *
     * slightly more efficient.                                                */

    /* Gets initial cost and loads bounding boxes. */
    double bb_cost = comp_bb_cost(NORMAL,
                                 placer_opts.place_cost_type,
                                 placer_opts.num_regions);
    int    num_connections = 0;
    double crit_exponent = 0.0;
    double max_delay = 0.0;
    double place_delay_value = 0.0;
    double timing_cost, delay_cost;
    double  inverse_prev_bb_cost = 0.0;
    double inverse_prev_timing_cost = 0.0;
    double cost = 0.0;
    int inet = 0;
    int ipin = 0;
    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        crit_exponent = placer_opts.td_place_exp_first; /*this will be modified when range_limit starts to change */

        compute_net_pin_index_values();

        num_connections = count_connections();
        printf("\nThere are %d point to point connections in this circuit\n\n",
                 num_connections);

        if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE) {
            for (inet = 0; inet < num_nets; ++inet)
                for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                    timing_place_crit[inet][ipin] = 0;    /*dummy crit values */
                }
           /* first pass gets delay_cost, which is used in criticality *
            * computations in the next call to comp_td_costs. */
            comp_td_costs(&timing_cost,
                          &delay_cost);
            place_delay_value = delay_cost / num_connections; /*used for computing criticalities */
            load_constant_net_delay(net_delay,
                                    place_delay_value);
        } else {
            place_delay_value = 0;
        }

        if (placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        /* This keeps net_delay up to date with the same values that * 
         * the placer is using point_to_point_delay_cost is computed *
         * each that comp_td_costs is called, and is also updated    *
         * after any swap is accepted   */
            net_delay = point_to_point_delay_cost;  
        }

        load_timing_graph_net_delays(net_delay);
        max_delay = load_net_slack(net_slack, 0);
        load_criticalities(placer_opts, net_slack, max_delay, crit_exponent);

        /*now we can properly compute costs  */
        comp_td_costs(&timing_cost,
                      &delay_cost);

        inverse_prev_timing_cost = 1 / timing_cost;
        inverse_prev_bb_cost = 1 / bb_cost;
        cost = 1; /* our new cost function uses normalized values of  */
        /*bb_cost and timing_cost, the value of cost will be reset  */
        /*to 1 at each temperature when *_TIMING_DRIVEN_PLACE is true */
    } else {
        /*BOUNDING_BOX_PLACE */
        cost = bb_cost;
        timing_cost = 0;
        delay_cost = 0;
        num_connections = 0;
        max_delay = 0;
        place_delay_value = 0;
        crit_exponent = 0;

        inverse_prev_timing_cost = 0;   /*inverses not used */
        inverse_prev_bb_cost = 0;
    }

    /* Sometimes I want to run the router with a random placement.  Avoid *
     * using 0 moves to stop division by 0 and 0 length vector_t problems,  *
     * by setting move_lim to 1 (which is still too small to do any       *
     * significant optimization).                                         */
    int move_lim = (int)(10 * pow(num_blocks, 1.3333));
    if (move_lim <= 0) {
        move_lim = 1;
    }

    double range_limit = (double)max(num_grid_columns, num_grid_rows);
    double t = starting_t(&cost, &bb_cost, &timing_cost,
                          placer_opts.place_cost_type,
                          old_region_occ_x, old_region_occ_y,
                          placer_opts.num_regions, fixed_pins, annealing_sched,
                          move_lim, range_limit, placer_opts.place_algorithm,
                          placer_opts.timing_tradeoff, inverse_prev_bb_cost,
                          inverse_prev_timing_cost, &delay_cost);
    int tot_iter = 0;
    printf("Initial Placement Cost: %g bb_cost: %g td_cost: %g delay_cost: %g\n\n",
            cost, bb_cost, timing_cost, delay_cost);

#ifndef SPEC
    printf
    ("%11s  %10s %11s  %11s  %11s %11s  %11s %9s %8s  %7s  %7s  %10s  %7s\n",
     "T", "Cost", "Av. BB Cost", "Av. TD Cost", "Av Tot Del",
     "P to P Del", "max_delay", "Ac Rate", "Std Dev", "R limit", "Exp",
     "Tot. Moves", "Alpha");
    printf
    ("%11s  %10s %11s  %11s  %11s %11s  %11s %9s %8s  %7s  %7s  %10s  %7s\n",
     "--------", "----------", "-----------", "-----------",
     "---------", "----------", "-----", "-------", "-------",
     "-------", "-------", "----------", "-----");
#endif

    char msg[BUFSIZE];
    sprintf(msg,
            "Initial Placement.  Cost: %g  BB Cost: %g  TD Cost %g  Delay Cost: %g "
            "\t max_delay %g Channel Factor: %d", cost, bb_cost, timing_cost,
            delay_cost, max_delay, width_fac);
    update_screen(MAJOR, msg, PLACEMENT, FALSE);

    clock_t start_cpu = clock();
    struct  timeval start;
    gettimeofday(&start, NULL);
    my_srandom(0);

    int inner_iter_num = annealing_sched.inner_num;

    pthread_mutex_init(&global_data_access.mutex, NULL);
    barrier_polling_reset();

    /* horizon_regions is regions in horizontal direction, verti_regions is  *
     * regions in vertical direction. Horizon_regions and verti_regions also *
     * means that threads used in horizontal or vertical direction.          */
    const int horizon_regions = x_partition[NUM_OF_THREADS - 1];
    const int verti_regions = y_partition[NUM_OF_THREADS - 1];
    if (horizon_regions < 1 || verti_regions < 1 || NUM_OF_THREADS > 64) {
        printf("This program cannot be ran with %d threads, try another configuration.\n", NUM_OF_THREADS);
        assert (0);
    }
    assert(horizon_regions * verti_regions == NUM_OF_THREADS);

    /* iSize{X,Y} is the number of cols/rows assigned to each threads.
     * extra_cols is used to distribute the extra rows/cols when {x,y}_part is
     * not evenly divisible by {x,y}_part */
    const int cols_assign_to_thread = (num_grid_columns + 2) / horizon_regions;
    const int rows_assign_to_thread = (num_grid_rows + 2) / verti_regions;
    const int extra_cols =
                (num_grid_columns + 2) - cols_assign_to_thread * horizon_regions;
    const int extra_rows =
                (num_grid_rows + 2) - rows_assign_to_thread * verti_regions;

    int success_sum, exit;
    double success_rat;
    double av_cost, av_bb_cost, av_timing_cost, av_delay_cost, sum_of_squares;
    double std_dev;

    pthread_data_t thread_data_array[NUM_OF_THREADS];
    pthread_t place_threads[NUM_OF_THREADS];
    int count = -1;
    for (count = NUM_OF_THREADS - 1; count >= 0 ; --count) {
        /* for x_start */
        if (count % horizon_regions ==  0 || horizon_regions == 1) {
            thread_data_array[count].x_start = 0;
        } else if (count % horizon_regions == 1) {
            thread_data_array[count].x_start = cols_assign_to_thread;
        } else if (count % horizon_regions < extra_cols + 1) {
            thread_data_array[count].x_start =
                (count % horizon_regions - 1) * (cols_assign_to_thread + 1)
                    + cols_assign_to_thread;
        } else {
            thread_data_array[count].x_start =
                (count % horizon_regions) * cols_assign_to_thread + extra_cols;
        }
        /* for x_end */
        if (count % horizon_regions == horizon_regions - 1 || horizon_regions == 1) {
            thread_data_array[count].x_end = num_grid_columns + 2;
        } else if (count % horizon_regions == 0) {
            thread_data_array[count].x_end = cols_assign_to_thread;
        } else if (count % horizon_regions < extra_cols + 1) {
            thread_data_array[count].x_end =
                (count % horizon_regions) * (cols_assign_to_thread + 1)
                    + cols_assign_to_thread;
        } else {
            thread_data_array[count].x_end =
             (count % horizon_regions + 1) * cols_assign_to_thread + extra_cols;
        }
        /* for y_start */
        if (count / horizon_regions == 0 || verti_regions == 1) {
            thread_data_array[count].y_start = 0;
        } else if (count / horizon_regions == 1) {
            thread_data_array[count].y_start = rows_assign_to_thread;
        } else if (count / horizon_regions < extra_rows + 1) {
            thread_data_array[count].y_start =
                (count / horizon_regions - 1) * (rows_assign_to_thread + 1)
                    + rows_assign_to_thread;
        } else {
            thread_data_array[count].y_start =
                (count / horizon_regions) * rows_assign_to_thread + extra_rows;
        }
        /* for y_end */
        if (count / horizon_regions == verti_regions - 1 || verti_regions == 1) {
            thread_data_array[count].y_end = num_grid_rows + 2;
        } else if (count / horizon_regions == 0) {
            thread_data_array[count].y_end = rows_assign_to_thread;
        } else if (count / horizon_regions < extra_rows + 1) {
            thread_data_array[count].y_end =
                (count / horizon_regions) * (rows_assign_to_thread + 1)
                     + rows_assign_to_thread;
        } else {
            thread_data_array[count].y_end =
                ((count / horizon_regions) + 1) * rows_assign_to_thread + extra_rows;
        }
        /* Initial region boundary OK */
        /* each region must be at least 8 X 8 */
        assert(thread_data_array[count].x_end - thread_data_array[count].x_start > 8);
        assert(thread_data_array[count].y_end - thread_data_array[count].y_start > 8);

        thread_data_array[count].thread_id = count;
        /* For each thread, it had these following global varibbles local copy, *
         * so it must synchronous with pthread mutex */
        thread_data_array[count].placer_opts = placer_opts;
        thread_data_array[count].annealing_sched = annealing_sched;
        thread_data_array[count].fixed_pins = fixed_pins;
        thread_data_array[count].net_slack = net_slack;
        thread_data_array[count].net_delay = net_delay;

        thread_data_array[count].t = &t;

        thread_data_array[count].bb_cost = &bb_cost;
        thread_data_array[count].timing_cost = &timing_cost;
        thread_data_array[count].delay_cost = &delay_cost;
        thread_data_array[count].inverse_prev_bb_cost = &inverse_prev_bb_cost;
        thread_data_array[count].inverse_prev_timing_cost = &inverse_prev_timing_cost;
        thread_data_array[count].av_cost = &av_cost;
        thread_data_array[count].av_bb_cost = &av_bb_cost;
        thread_data_array[count].av_timing_cost = &av_timing_cost;
        thread_data_array[count].av_delay_cost = &av_delay_cost;
        thread_data_array[count].cost = &cost;

        thread_data_array[count].move_lim = &move_lim;
        thread_data_array[count].inner_iter_num = &inner_iter_num;
        thread_data_array[count].tot_iter = &tot_iter;
        thread_data_array[count].success_sum = &success_sum;
        thread_data_array[count].success_rat = &success_rat;
        thread_data_array[count].sum_of_squares = &sum_of_squares;
        thread_data_array[count].std_dev = &std_dev;

        thread_data_array[count].place_delay_value = &place_delay_value;
        thread_data_array[count].max_delay = &max_delay;
        thread_data_array[count].num_connections = &num_connections;

        thread_data_array[count].crit_exponent = &crit_exponent;
        thread_data_array[count].exit = &exit;
        thread_data_array[count].range_limit = &range_limit;

        assert(thread_data_array[count].y_end > thread_data_array[count].y_start);

        if (count != 0) {
             pthread_create(&place_threads[count],
                            NULL,
                            try_place_parallel,
                            (void*)&thread_data_array[count]);
        } else {
            /* Why count == 0, it needn't create pthread? */
            try_place_parallel(&thread_data_array[count]);
        }
    } /* end of for(count = NUM_OF_THREADS-1; count >= 0; --count) */

    for (count = NUM_OF_THREADS - 1; count > 0; --count) {
        pthread_join(place_threads[count], NULL);
    }

    pthread_mutex_destroy(&global_data_access.mutex);

    clock_t finish_cpu = clock();
    struct  timeval finish;
    gettimeofday(&finish, NULL);
    bb_cost = check_place(bb_cost,
                          timing_cost,
                          placer_opts.place_cost_type,
                          placer_opts.num_regions,
                          placer_opts.place_algorithm,
                          delay_cost);

    if (placer_opts.enable_timing_computations &&
            placer_opts.place_algorithm == BOUNDING_BOX_PLACE) {
        /*need this done since the timing data has not been kept up to date*
         *in bounding_box mode */
        for (inet = 0; inet < num_nets; ++inet)
            for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                timing_place_crit[inet][ipin] = 0;    /*dummy crit values */
            }

        comp_td_costs(&timing_cost, &delay_cost);   /*computes point_to_point_delay_cost */
    }

    double place_est_crit_delay = 0.0;
    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE ||
            placer_opts.enable_timing_computations) {
        net_delay = point_to_point_delay_cost;  /*this makes net_delay up to date with    *
                             *the same values that the placer is using*/
        load_timing_graph_net_delays(net_delay);
        place_est_crit_delay = load_net_slack(net_slack, 0);
#ifdef CREATE_ECHO_FILES
        /*      print_sink_delays("placement_sink_delays.echo"); */
        print_net_slack("placement_net_slacks.echo", net_slack);
        print_critical_path("placement_crit_path.echo",
                            *subblock_data_ptr);
#endif /* CREATE_ECHO_FILES */
        printf("Placement Estimated Crit Path Delay: %g\n\n",
               place_est_crit_delay);
    }

    sprintf(msg,
            "Placement. Cost: %g  bb_cost: %g td_cost: %g Channel Factor: %d max_delay: %g",
            cost, bb_cost, timing_cost, width_fac, max_delay);
    printf("Placement. Cost: %g  bb_cost: %g  td_cost: %g  delay_cost: %g.\n",
           cost, bb_cost, timing_cost, delay_cost);

    printf("inner loop wall: %f sec, cpu total: %f\n",
            my_difftime2(&start, &finish),
            (double)(finish_cpu - start_cpu) / CLOCKS_PER_SEC);
    update_screen(MAJOR, msg, PLACEMENT, FALSE);

#ifdef SPEC
    printf("Total moves attempted: %d.0\n", tot_iter);
#endif

    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE ||
            placer_opts.enable_timing_computations) {
        net_delay = remember_net_delay_original_ptr;
        free_placement_structs(placer_opts.place_cost_type,
                               placer_opts.num_regions, old_region_occ_x,
                               old_region_occ_y, placer_opts);
        free_lookups_and_criticalities(&net_delay, &net_slack);
    }

    /* placement is done - find mst of all nets.
     * creating mst for each net; this gives me an ordering of sinks
     * by which I will direct search (A*) for. */

    /* mst not needed for place_only circuits */
    if (operation != PLACE_ONLY) {
        if (*mst) {
            for (inet = 0; inet < num_nets; ++inet) {
                assert((*mst)[inet]);
                free((*mst)[inet]);
            }
            free(*mst);
        }

        *mst = (t_mst_edge**)my_malloc(sizeof(t_mst_edge*) * num_nets);
        for (inet = 0; inet < num_nets; ++inet) {
            (*mst)[inet] = get_mst_of_net(inet);
        }
    }

    free(x_lookup);
}  /* end of try_place_use_multi_threads() */


/* ====== parallel placement functions ======*/
/* polling barriers */
/* Binary tree based barrier except node 0 and 1.
 * Barrier first waits for its children to arrive, then
 * wait for its parent to release them.
 * Sample structure for a tree with 8 nodes:
 *                0
 *                |
 *                1
 *               / \
 *              2   3
 *             / \  /\
 *            4  5 6  7
 */
void barrier_polling(int kthread_id)
{
    ++(barrier1[kthread_id].entry);

    const int local_entry = (barrier1[kthread_id].entry ) % 2;
    if (NUM_OF_THREADS == 1) {
        return;
    } else if (kthread_id == 0) {
        /* wait for kthread_id 1 finished */
        while (barrier1[1].arrived != local_entry) {};

        barrier1[1].proceed = local_entry;
    } else if (kthread_id == 1) {
        /* wait for children finished */
        if (NUM_OF_THREADS > 3) {
            while (barrier1[2].arrived != local_entry || barrier1[3].arrived != local_entry) {};
        } else if (NUM_OF_THREADS > 2) {
            while (barrier1[2].arrived != local_entry) {};
        }

        /* signal my arrival */
        barrier1[kthread_id].arrived = local_entry;

        while (barrier1[kthread_id].proceed != local_entry) {};

        //release children
        if (NUM_OF_THREADS > 3) {
            barrier1[2].proceed = barrier1[3].proceed = local_entry;
        } else if (NUM_OF_THREADS > 2) {
            barrier1[2].proceed = local_entry;
        }
    } else { /* kthread_id > 1 */
        /* wait for children finished! */
        if (NUM_OF_THREADS > kthread_id * 2 + 1) {
            while (barrier1[kthread_id * 2].arrived != local_entry
                    || barrier1[kthread_id * 2 + 1].arrived != local_entry) {};
        } else if (NUM_OF_THREADS > kthread_id * 2) {
            while (barrier1[kthread_id * 2].arrived != local_entry) {};
        }

        //signal my arrival
        barrier1[kthread_id].arrived = local_entry;

        while (barrier1[kthread_id].proceed != local_entry) {};

        //release children
        if (NUM_OF_THREADS > kthread_id * 2 + 1) {
            barrier1[kthread_id * 2].proceed =
                barrier1[kthread_id * 2 + 1].proceed = local_entry;
        } else if (NUM_OF_THREADS > kthread_id * 2) {
            barrier1[kthread_id * 2].proceed = local_entry;
        }
    } /* end of else */
}  /* end of void barrier_polling(int kthread_id) */

void barrier_polling_reset()
{
    int x = -1;
    for (x = 0; x < NUM_OF_THREADS; ++x) {
        barrier1[x].arrived = 0;
        barrier1[x].proceed = 0;
        barrier1[x].entry = 0;
    }
} /* end of void barrier_polling_reset() */


/* Placement using multi-threads parallely */
static void* try_place_parallel(void* args)
{
    pthread_data_t* input_args = (pthread_data_t*)args;

    int core_affinity = pow(2,
                            input_args->thread_id);
    /*core affinity - locking threads to a particular core*/
    if (pthread_setaffinity_np(pthread_self(),
                               sizeof(core_affinity),
                               &core_affinity) < 0) {
        printf("core affinity error: thread-%d couldn't get Core-%d\n",
               input_args->thread_id, core_affinity);
        exit(-1);
    }

    thread_local_common_paras_t  common_paras;
    int max_pins_per_fb = 0;
    /* tp_initialize() initial all thread_data variables */
    tp_initialize(input_args,
                  &max_pins_per_fb,
                  &common_paras);

    thread_local_data_for_swap_t  swap_data;
    /* allocate space for local data */ 
    tp_alloc_mem(max_pins_per_fb,
                 common_paras.local_placer_opts,
                 &swap_data);

    /* Initial all the 7 local parameters */ 
    tp_init_local(common_paras.local_region_x_boundary,
                  &swap_data);

    const int kthread_id = common_paras.local_thread_id;
    /* build the fan-in tree for each node for timing update */
    find_fanin_parallel(kthread_id);

    if (pthread_mutex_trylock(&global_data_access.mutex) != 0) {
        barrier_polling(kthread_id);
    } else {
    /*another way of storing data, in columns instead of rows. Why? */
        tp_init_localvert_grid();

        barrier_polling(kthread_id);

        pthread_mutex_unlock(&global_data_access.mutex);
    }

    /* program specific parameters */
    int prob = PROB; /* 10 */
    int timing_update_threshold = TIMING_UPDATE; /* 5 */
    /* ensures timing update will occurduring first iteration */
    int freq = timing_update_threshold; 

    /*=========    NOW STARTING MAIN PLACEMENT  =======*/
    int iter = 0;
    int move_counter = -1;
    while (common_paras.local_temper != 0.0) {
        /*  First update parameters that used for controlloing placement! */
        common_paras.local_temper = *(input_args->t);
        int inner_iter_num = *(input_args->inner_iter_num);
        common_paras.local_crit_exponent = *(input_args->crit_exponent);
        common_paras.local_success_ratio = *(input_args->success_rat);
        common_paras.local_range_limit = *(input_args->range_limit);

        for (iter = 0; iter < inner_iter_num; ++iter) {
        /* FIXME, wait for all processors to finish the previous iteration */
            barrier_polling(kthread_id);

            common_paras.local_delay_cost = *(input_args->delay_cost);
            /* place_delay_value = delay_cost / num_connections; */
            common_paras.local_place_delay_value =
                common_paras.local_delay_cost / common_paras.local_num_conns;

            /*parallel timing update. Only run timing update once per
             * 'timing_update_threadhold' iterations */
            if (freq >= timing_update_threshold) {
                /* It first calculate all tnodes arr_time and req_time parallel */
                tp_timing_update_full(kthread_id,
                                      &(common_paras.local_net_delay),
                                      input_args,
                                      &(common_paras.local_max_delay));
                /* then calculate all timing_edge' slack value. */
                tp_compute_net_slack(kthread_id,
                                     &(common_paras.local_net_slack),
                                     &(common_paras.local_timing_cost),
                                     &(common_paras.local_delay_cost),
                                     common_paras.local_crit_exponent,
                                     common_paras.local_max_delay,
                                     input_args);

                /*reset frequency counter*/
                freq = 0;
            } else {
                /*only let one thread do reset update global variable */
                if (pthread_mutex_trylock(&global_data_access.mutex) != 0) {
                    barrier_polling(kthread_id);
                } else {
                    *(input_args->timing_cost) = 0.0;
                    *(input_args->delay_cost) = 0.0;
                    barrier_polling(kthread_id);
                    /* unlock mutex */
                    pthread_mutex_unlock(&global_data_access.mutex);
                }

                /* no criticality update, but still need to recalculate timing */
                tp_timing_calc(kthread_id,
                               &(common_paras.local_timing_cost),
                               &(common_paras.local_delay_cost),
                               input_args);
                ++freq;
            }

            barrier_polling(kthread_id);

            tp_iter_data_update(iter,
                                input_args,
                                &common_paras,
                                &swap_data);

            /* clear local counters. But why? */
            if (iter == 0) {
                common_paras.local_success_sum = 0;
                move_counter = 0;
                common_paras.local_av_cost = 0.0;
                common_paras.local_av_bb_cost = 0.0;
                common_paras.local_av_timing_cost = 0.0;
                common_paras.local_av_delay_cost = 0.0;
                common_paras.local_sum_of_squares = 0.0;
            }

            barrier_polling(kthread_id);
            /* Iterate through each sub-regions of the grid. FIXME, each region
             * had 2x2 sub-regions. */
            int row, col;
            for (row = 0; row < 2; ++row) {
                for (col = 0; col < 2; ++col) {
                    try_place_a_subregion(kthread_id,
                                          row, col,
                                          prob,
                                          &move_counter,
                                          &common_paras,
                                          &swap_data);
                }  /* end of for(col = 0; col < 2; ++col) */
            }  /* end of for(row = 0; row < 2; ++row) */

            /*prepare for bb_box calculation*/
            if (pthread_mutex_trylock(&global_data_access.mutex) != 0) {
                barrier_polling(kthread_id);
            } else {
                *(input_args->bb_cost) = 0.;
                barrier_polling(kthread_id);
                pthread_mutex_unlock(&global_data_access.mutex);
            }

            /*bb cost update
             *parallel calculation followed by serial addition in order to avoid
             *floating point roundoff error
             */
            double wirelength_cost =
                comp_bb_cost_parallel(start_finish_nets[kthread_id].start_sinks,
                                      start_finish_nets[kthread_id].finish_sinks);
            common_paras.local_bb_cost = wirelength_cost;
            partial_results[kthread_id] = wirelength_cost;

            barrier_polling(kthread_id);
            if (kthread_id == 0) {
                wirelength_cost = 0.0;
                int update_count = -1;
                for (update_count = 0; update_count < NUM_OF_THREADS; ++update_count) {
                    wirelength_cost += partial_results[update_count];
                }

                *(input_args->bb_cost) = wirelength_cost;
            }
        } /* end of for(iter = 0; iter < inner_iter_num; ++iter) */

        /*synchronize result*/
        pthread_mutex_lock(&global_data_access.mutex);
        *(input_args->success_sum ) += common_paras.local_success_sum;
        *(input_args->move_lim) += move_counter;
        *(input_args->av_cost) += common_paras.local_av_cost;
        *(input_args->av_bb_cost) += common_paras.local_av_bb_cost;
        *(input_args->av_timing_cost) += common_paras.local_av_timing_cost;
        *(input_args->av_delay_cost) += common_paras.local_av_delay_cost;
        *(input_args->sum_of_squares) += common_paras.local_sum_of_squares;
        pthread_mutex_unlock(&global_data_access.mutex);

        barrier_polling(kthread_id);
        /*data gathering and printing*/
        tp_data_print_to_screen(kthread_id,
                                &common_paras.local_timing_cost,
                                &common_paras.local_delay_cost,
                                &common_paras.local_crit_exponent,
                                common_paras.local_max_delay,
                                input_args,
                                common_paras.local_placer_opts,
                                &common_paras.local_total_cost,
                                &common_paras.local_bb_cost,
                                &common_paras.local_success_sum,
                                &common_paras.local_sum_of_squares,
                                &move_counter,
                                &common_paras.local_total_iter,
                                &common_paras.local_success_ratio,
                                &common_paras.local_av_cost,
                                &common_paras.local_av_bb_cost,
                                &common_paras.local_av_timing_cost,
                                &common_paras.local_av_delay_cost,
                                &common_paras.local_std_dev,
                                &common_paras.local_range_limit,
                                &common_paras.local_final_rlim,
                                &common_paras.local_inverse_delta_rlim,
                                &common_paras.local_temper,
                                &common_paras.local_old_temper,
                                &inner_iter_num,
                                common_paras.local_place_delay_value);
        if (common_paras.local_temper != 0.0) {
            barrier_polling(common_paras.local_thread_id);
        }
    }  /* end of while (t != 0.0) */

    /* free resources */
    int inet = -1;
    const place_algorithm_t kplace_algorithm =
        common_paras.local_placer_opts.place_algorithm;
    if (kplace_algorithm == NET_TIMING_DRIVEN_PLACE ||
            kplace_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        for (inet = 0; inet < num_nets; ++inet) {
            /*add one to the address since it is indexed from 1 not 0 */
            (swap_data.m_local_temp_point_to_point_delay_cost[inet])++;
            free(swap_data.m_local_temp_point_to_point_delay_cost[inet]);

            swap_data.m_local_temp_point_to_point_timing_cost[inet]++;
            free(swap_data.m_local_temp_point_to_point_timing_cost[inet]);
        }

        free(swap_data.m_local_temp_point_to_point_delay_cost);
        swap_data.m_local_temp_point_to_point_delay_cost = NULL;

        free(swap_data.m_local_temp_point_to_point_timing_cost);
        swap_data.m_local_temp_point_to_point_timing_cost = NULL;
    }

    free(swap_data.m_local_block);
    swap_data.m_local_block = NULL;
    free(swap_data.m_local_temp_net_cost);
    swap_data.m_local_temp_net_cost = NULL;
    free(swap_data.m_local_net_cost);
    swap_data.m_local_net_cost = NULL;

    free(swap_data.m_local_bb_coord);
    swap_data.m_local_bb_coord = NULL;
    free(swap_data.m_local_bb_edge);
    swap_data.m_local_bb_edge = NULL;

    free(swap_data.m_nets_to_update);
    swap_data.m_nets_to_update = NULL;

    free(swap_data.m_net_block_moved);
    swap_data.m_net_block_moved = NULL;

    free(swap_data.m_bb_coord_new);
    swap_data.m_bb_coord_new = NULL;

    free(swap_data.m_bb_edge_new);
    swap_data.m_bb_edge_new = NULL;
} /* end of void* try_place_parallel(void* args) */


static void  try_place_a_subregion(const int  kthread_id,
                                   const int  krow,
                                   const int  kcol,
                                   const int  prob,
                                   int*     move_counter,
                                   thread_local_common_paras_t*  common_paras_ptr,
                                   thread_local_data_for_swap_t* swap_data_ptr)
{
    barrier_polling(kthread_id);

    /* update locate Sub-Region from global data, and first
     * update horizontal data, then update vertical data. */
    const int* kregion_x_boundary = common_paras_ptr->local_region_x_boundary;
    const int* kregion_y_boundary = common_paras_ptr->local_region_y_boundary;
    grid_tile_t**  local_grid = swap_data_ptr->m_local_grid;
    local_block_t* local_block = swap_data_ptr->m_local_block;
    tp_local_data_update(kregion_x_boundary,
                         kregion_y_boundary,
                         local_grid,
                         local_block,
                         krow, kcol);

    /*sequentially consider each block within the sub-region*/
    const int hori_start_bound = kregion_x_boundary[krow];
    const int hori_end_bound  =  kregion_x_boundary[krow + 1];
    const int vert_start_bound = kregion_y_boundary[kcol];
    const int vert_end_bound  =  kregion_y_boundary[kcol + 1];
    const boolean kfixed_pins = common_paras_ptr->local_fixed_pins;
    const placer_opts_t kplacer_opts = common_paras_ptr->local_placer_opts;
    const double kt = common_paras_ptr->local_temper;
    int x = 0;
    int y = 0;
    int z = 0;
    for (x = hori_start_bound; x < hori_end_bound; ++x) {
        for (y = vert_start_bound; y < vert_end_bound; ++y) {
            if ((local_grid[x][y].type == EMPTY_TYPE)
                  || (kfixed_pins && local_grid[x][y].type == IO_TYPE)) {
                continue;
            }
            const int kcapacity = local_grid[x][y].type->capacity;
            for (z = 0; z < kcapacity; ++z) {
                //do not consider empty locations, ie - four corners of the grid or fixed pins
                if (local_grid[x][y].blocks[z] == EMPTY) {
                    continue;
                }

                /* PROB_SKIPPED */
                if (my_irand_parallel(100, kthread_id) >= prob) {
                    ++(*move_counter);
                    int return_val = try_swap_parallel(kt, kthread_id,
                                                       x, y, z, local_grid[x][y].blocks[z],
                                                       kplacer_opts.place_cost_type,
                                                       kplacer_opts.place_algorithm,
                                                       kplacer_opts.timing_tradeoff,
                                                       &(common_paras_ptr->local_total_cost),
                                                       &(common_paras_ptr->local_bb_cost),
                                                       &(common_paras_ptr->local_timing_cost),
                                                       &(common_paras_ptr->local_delay_cost),
                                                       common_paras_ptr->local_inverse_prev_bb_cost,
                                                       common_paras_ptr->local_inverse_prev_timing_cost,
                                                       common_paras_ptr->local_range_limit,
                                                       max(kregion_x_boundary[krow] - 2, 0),
                                                       min(kregion_x_boundary[krow + 1] + 1,
                                                           num_grid_columns + 1), /* x direction range */
                                                       max(kregion_y_boundary[kcol] - 2, 0),
                                                       min(kregion_y_boundary[kcol + 1] + 1,
                                                           num_grid_rows + 1), /* y direction range */
                                                       swap_data_ptr);
                    if (return_val == 1) {
                        ++(common_paras_ptr->local_success_sum);
                        const double ktotal_cost = common_paras_ptr->local_total_cost;
                        common_paras_ptr->local_av_cost += ktotal_cost;
                        common_paras_ptr->local_av_bb_cost += common_paras_ptr->local_bb_cost;
                        common_paras_ptr->local_av_timing_cost += common_paras_ptr->local_timing_cost;
                        common_paras_ptr->local_av_delay_cost += common_paras_ptr->local_delay_cost;
                        common_paras_ptr->local_sum_of_squares +=
                            ktotal_cost * ktotal_cost;
                    }
                } /* end of if(my_irand_parallel(100, kthread_id) >= prob) */
            } /* end of for(z = 0; z < local_grid[x][y].type->capacity; ++z) */
        } /* end of for(y = vert_start_bound; y < vert_end_bound; ++y) */
    } /* end of for(x = hori_start_bound; x < hori_end_bound; ++x) */

    /*update global data with local changes*/
    update_from_local_to_global(local_block,
                                local_grid,
                                max(0, kregion_x_boundary[krow] - 2),
                                min(num_grid_columns + 2,
                                    kregion_x_boundary[krow + 1] + 2),
                                max(0, kregion_y_boundary[kcol] - 2),
                                min(num_grid_rows + 2,
                                    kregion_y_boundary[kcol + 1] + 2));
} /* end of static void try_place_a_subregion() */

/* FIXME, Important funtion! It initial all imporant parameters that used for *
 * parallel placement.                                                        */
static void tp_initialize(pthread_data_t*  input_args,
                          int*  max_pins_per_fb,
                          thread_local_common_paras_t*  common_paras_ptr)
{
    /*thread-dependent data*/
    common_paras_ptr->local_thread_id = input_args->thread_id;
    common_paras_ptr->local_y_start = input_args->y_start;
    common_paras_ptr->local_y_end = input_args->y_end;
    common_paras_ptr->local_x_start = input_args->x_start;
    common_paras_ptr->local_x_end = input_args->x_end;

    /*placer-depdent data*/
    common_paras_ptr->local_fixed_pins = input_args->fixed_pins;
    common_paras_ptr->local_temper = *(input_args->t);
    common_paras_ptr->local_placer_opts = input_args->placer_opts;

    /* square based partitioning boundary */
    const int ky_start = common_paras_ptr->local_y_start;
    const int ky_end = common_paras_ptr->local_y_end;
    const int kx_start = common_paras_ptr->local_x_start;
    const int kx_end = common_paras_ptr->local_x_end;
    common_paras_ptr->local_region_y_boundary[0] = ky_start;
    common_paras_ptr->local_region_y_boundary[1] =
        ky_start + ((int)(ky_end - ky_start) / 2);
    common_paras_ptr->local_region_y_boundary[2] = ky_end;

    common_paras_ptr->local_region_x_boundary[0] = kx_start;
    common_paras_ptr->local_region_x_boundary[1] =
        kx_start + ((int)(kx_end - kx_start) / 2);
    common_paras_ptr->local_region_x_boundary[2] = kx_end;
 
    /* other initializations */
    common_paras_ptr->local_first_rlim = (double)max(num_grid_columns,
                                                     num_grid_rows);
    common_paras_ptr->local_range_limit = common_paras_ptr->local_first_rlim;
    common_paras_ptr->local_final_rlim = 1;
    common_paras_ptr->local_inverse_delta_rlim =
        1 / (common_paras_ptr->local_first_rlim - common_paras_ptr->local_final_rlim);

    common_paras_ptr->local_place_delay_value = *(input_args->place_delay_value);
    common_paras_ptr->local_max_delay = *(input_args->max_delay);
    common_paras_ptr->local_num_conns = *(input_args->num_connections);
    common_paras_ptr->local_delay_cost = *(input_args->delay_cost);
    common_paras_ptr->local_net_slack = input_args->net_slack;
    common_paras_ptr->local_net_delay = input_args->net_delay;

    /* max_pins_per_fb initialization*/
    int i = -1;
    for (i = 0; i < num_types; ++i) {
        *max_pins_per_fb = max(*max_pins_per_fb,
                               type_descriptors[i].num_pins);
    }

    /* dynamic workload distribution for timing update initialization*/
    const int knets_assign_to_thread = ceil((double)num_nets / NUM_OF_THREADS);
    const int kthread_id = common_paras_ptr->local_thread_id;
    start_finish_nets[kthread_id].start_edge = kthread_id * knets_assign_to_thread;
    start_finish_nets[kthread_id].finish_edge =
        min(num_nets, (kthread_id + 1) * knets_assign_to_thread);

    start_finish_nets[kthread_id].start_sinks =
        start_finish_nets[kthread_id].start_edge;
    start_finish_nets[kthread_id].finish_sinks =
        start_finish_nets[kthread_id].finish_edge;

    start_finish_nets[kthread_id].edge_partition_size =
        (start_finish_nets[kthread_id].finish_edge
             - start_finish_nets[kthread_id].start_edge);
    start_finish_nets[kthread_id].sink_partition_size =
        (start_finish_nets[kthread_id].finish_sinks
             - start_finish_nets[kthread_id].start_sinks);

    start_finish_nets[kthread_id].counter_sink = 0;
    start_finish_nets[kthread_id].counter_edge = 0;
} /* end of void tp_initialize() */

static void tp_alloc_mem(const int max_pins_per_fb,
                         const placer_opts_t  placer_opts,
                         thread_local_data_for_swap_t*  swap_data_ptr)
{
    /* local_grid[col][row] */
    swap_data_ptr->m_local_grid =
            (grid_tile_t**)alloc_matrix(0, (num_grid_columns + 1),
                                        0, (num_grid_rows + 1),
                                        sizeof(grid_tile_t));

    swap_data_ptr->m_local_block =
            (local_block_t*)malloc(num_blocks * sizeof(local_block_t));

    /* Why did this array allocate 2 * max_pins_per_fb? */
    swap_data_ptr->m_bb_coord_new =
            (bbox_t*)my_malloc(2 * max_pins_per_fb * sizeof(bbox_t));
    swap_data_ptr->m_bb_edge_new =
            (bbox_t*)my_malloc(2 * max_pins_per_fb * sizeof(bbox_t));

    swap_data_ptr->m_nets_to_update =
            (int*)my_malloc(2 * max_pins_per_fb * sizeof(int));
    swap_data_ptr->m_net_block_moved =
            (int*)my_malloc(2 * max_pins_per_fb * sizeof(int));

    swap_data_ptr->m_local_bb_coord = (bbox_t*)my_malloc(num_nets * sizeof(bbox_t));
    swap_data_ptr->m_local_bb_edge = (bbox_t*)my_malloc(num_nets * sizeof(bbox_t));

    swap_data_ptr->m_local_temp_net_cost = (double*)malloc(num_nets * sizeof(double));
    swap_data_ptr->m_local_net_cost = (double*)malloc(num_nets * sizeof(double));
    /* point_to_point_delay_cost & point_to_point_timing_cost are not updated.
     * Hence, no local copies of these are needed, and global data is read
     * during placement. Global data is updated once per iter */
    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        swap_data_ptr->m_local_temp_point_to_point_delay_cost =
            (double**)my_malloc(num_nets * sizeof(double*));
        swap_data_ptr->m_local_temp_point_to_point_timing_cost =
            (double**)my_malloc(num_nets * sizeof(double*));

        int inet = -1;
        for (inet = 0; inet < num_nets; ++inet) {
            (swap_data_ptr->m_local_temp_point_to_point_delay_cost)[inet] =
                    (double*)my_malloc(net[inet].num_sinks * sizeof(double));
            --(swap_data_ptr->m_local_temp_point_to_point_delay_cost)[inet];

            (swap_data_ptr->m_local_temp_point_to_point_timing_cost)[inet] =
                (double*)my_malloc(net[inet].num_sinks * sizeof(double));
            --(swap_data_ptr->m_local_temp_point_to_point_timing_cost)[inet];
        }
    }
} /* end of static void tp_alloc_mem() */


static int count_connections()
{
    /*only count non-global connections */
    int count = 0;
    int inet = -1;
    for (inet = 0; inet < num_nets; inet++) {
        if (net[inet].is_global) {
            continue;
        }
        count += net[inet].num_sinks;
    }

    return count;
}

static void compute_net_pin_index_values()
{
    /*computes net_pin_index array, this array allows us to quickly */
    /*find what pin on the net a block pin corresponds to */

    int inet, netpin, blk, iblk, ipin;
    block_type_ptr type;

    /*initialize values to OPEN */
    for (iblk = 0; iblk < num_blocks; iblk++) {
        type = block[iblk].type;

        for (ipin = 0; ipin < type->num_pins; ipin++) {
            net_pin_index[iblk][ipin] = OPEN;
        }
    }

    for (inet = 0; inet < num_nets; inet++) {

        if (net[inet].is_global) {
            continue;
        }

        for (netpin = 0; netpin <= net[inet].num_sinks; netpin++) {
            blk = net[inet].node_block[netpin];
            net_pin_index[blk][net[inet].node_block_pin[netpin]] =
                netpin;
        }
    }
}

static double
get_std_dev(int n,
            double sum_x_squared,
            double av_x)
{

    /* Returns the standard deviation of data set x.  There are n sample points, *
     * sum_x_squared is the summation over n of x^2 and av_x is the average x.   *
     * All operations are done in double precision, since round off error can be *
     * a problem in the initial temp. std_dev calculation for big circuits.      */

    double std_dev;

    if (n <= 1) {
        std_dev = 0.;
    } else {
        std_dev = (sum_x_squared - n * av_x * av_x) / (double)(n - 1);
    }

    if (std_dev > 0.) {     /* Very small variances sometimes round negative */
        std_dev = sqrt(std_dev);
    } else {
        std_dev = 0.;
    }

    return (std_dev);
}


static void
update_rlim(double* range_limit,
            double success_rat)
{

    /* Update the range limited to keep acceptance prob. near 0.44.  Use *
     * a floating point range_limit to allow gradual transitions at low temps.  */

    double upper_lim;

    *range_limit = (*range_limit) * (1. - 0.60 + success_rat);
    upper_lim = max(num_grid_columns, num_grid_rows);
    *range_limit = min(*range_limit, upper_lim);
    *range_limit = max(*range_limit, 1.);
}


static int
exit_crit(double t,
          double cost,
          annealing_sched_t annealing_sched)
{

    /* Return 1 when the exit criterion is met.                        */

    if (annealing_sched.type == USER_SCHED) {
        if (t < annealing_sched.exit_t) {
            return (1);
        } else {
            return (0);
        }
    }

    /* Automatic annealing schedule */

    if (t < 0.005 * cost / num_nets) {
        return (1);
    } else {
        return (0);
    }
}


static double
starting_t(double* cost_ptr,
           double* bb_cost_ptr,
           double* timing_cost_ptr,
           int place_cost_type,
           double** old_region_occ_x,
           double** old_region_occ_y,
           int num_regions,
           boolean fixed_pins,
           annealing_sched_t annealing_sched,
           int max_moves,
           double range_limit,
           place_algorithm_t place_algorithm,
           double timing_tradeoff,
           double inverse_prev_bb_cost,
           double inverse_prev_timing_cost,
           double* delay_cost_ptr)
{

    /* Finds the starting temperature (hot condition).              */

    int i, num_accepted, move_lim;
    double std_dev, av, sum_of_squares; /* Double important to avoid round off */
    int* x_lookup;

    x_lookup = (int*)my_malloc(num_grid_columns * sizeof(int));

    if (annealing_sched.type == USER_SCHED) {
        return (annealing_sched.init_t);
    }

    move_lim = min(max_moves, num_blocks);

    num_accepted = 0;
    av = 0.;
    sum_of_squares = 0.;

    /* Try one move per block.  Set t high so essentially all accepted. */
    for (i = 0; i < move_lim; i++) {
        if (try_swap(1.e30, cost_ptr, bb_cost_ptr, timing_cost_ptr, range_limit,
                     place_cost_type,
                     old_region_occ_x, old_region_occ_y, num_regions,
                     fixed_pins, place_algorithm, timing_tradeoff,
                     inverse_prev_bb_cost, inverse_prev_timing_cost,
                     delay_cost_ptr, x_lookup) == 1) {
            num_accepted++;
            av += *cost_ptr;
            sum_of_squares += *cost_ptr * (*cost_ptr);
        }
    }

    if (num_accepted != 0) {
        av /= num_accepted;
    } else {
        av = 0.;
    }

    std_dev = get_std_dev(num_accepted, sum_of_squares, av);

#ifdef DEBUG

    if (num_accepted != move_lim) {
        printf
        ("Warning:  Starting t: %d of %d configurations accepted.\n",
         num_accepted, move_lim);
    }

#endif

#ifdef VERBOSE
    printf("std_dev: %g, average cost: %g, starting temp: %g\n",
           std_dev, av, 20. * std_dev);
#endif

    free(x_lookup);

    /* Set the initial temperature to 20 times the standard of deviation */
    /* so that the initial temperature adjusts according to the circuit */
    return (20. * std_dev);
}

/* Picks some block and moves it to another spot.  If this spot is   *
 * occupied, switch the blocks.  Assess the change in cost function  *
 * and accept or reject the move.  If rejected, return 0.  If        *
 * accepted return 1.  Pass back the new value of the cost function. *
 * range_limit is the range limiter.                                                                            */
static int try_swap(double t,
                    double* cost,
                    double* bb_cost,
                    double* timing_cost,
                    double range_limit,
                    int place_cost_type,
                    double** old_region_occ_x,
                    double** old_region_occ_y,
                    int num_regions,
                    boolean fixed_pins,
                    place_algorithm_t place_algorithm,
                    double timing_tradeoff,
                    double inverse_prev_bb_cost,
                    double inverse_prev_timing_cost,
                    double* delay_cost,
                    int* x_lookup)
{
    int i, k, inet;
    int max_pins_per_fb = 0;
    for (i = 0; i < num_types; ++i) {
        max_pins_per_fb =
            max(max_pins_per_fb, type_descriptors[i].num_pins);
    }

    /* Allocate the local bb_coordinate storage, etc. only once. */
    static bbox_t* bb_coord_new = NULL;
    static bbox_t* bb_edge_new = NULL;
    static int* nets_to_update = NULL;
    static int* net_block_moved = NULL;
    if (bb_coord_new == NULL) {
        bb_coord_new = (bbox_t*)my_malloc(2 * max_pins_per_fb * sizeof(bbox_t));
        bb_edge_new = (bbox_t*)my_malloc(2 * max_pins_per_fb * sizeof(bbox_t));
        nets_to_update = (int*)my_malloc(2 * max_pins_per_fb * sizeof(int));
        net_block_moved = (int*)my_malloc(2 * max_pins_per_fb * sizeof(int));
    }

    int from_block = my_irand(num_blocks - 1);
    /* If the pins are fixed we never move them from their initial    *
     * random locations.  The code below could be made more efficient *
     * by using the fact that pins appear first in the block list,    *
     * but this shouldn't cause any significant slowdown and won't be *
     * broken if I ever change the parser so that the pins aren't     *
     * necessarily at the start of the block list.                    */
    if (fixed_pins == TRUE) {
        while (block[from_block].type == IO_TYPE) {
            from_block = my_irand(num_blocks - 1);
        }
    }

    int x_from = block[from_block].x;
    int y_from = block[from_block].y;
    int z_from = block[from_block].z;
    int x_to = 0;
    int y_to = 0;
    if (!find_to(x_from,
                 y_from,
                 block[from_block].type,
                 range_limit,
                 x_lookup,
                 &x_to,
                 &y_to)) {
        return FALSE;
    }

    /* Make the switch in order to make computing the new bounding *
     * box simpler.  If the cost increase is too high, switch them *
     * back.  (block data structures switched, clbs not switched   *
     * until success of move is determined.)                       */
    int z_to = 0;
    if (grid[x_to][y_to].type->capacity > 1) {
        z_to = my_irand(grid[x_to][y_to].type->capacity - 1);
    }

    int to_block = EMPTY;
    if (grid[x_to][y_to].blocks[z_to] == EMPTY) {
        /* Moving to an empty location */
        to_block = EMPTY;
        block[from_block].x = x_to;
        block[from_block].y = y_to;
        block[from_block].z = z_to;
    } else {
        /* Swapping two blocks */
        to_block = grid[x_to][y_to].blocks[z_to];
        block[to_block].x = x_from;
        block[to_block].y = y_from;
        block[to_block].z = z_from;

        block[from_block].x = x_to;
        block[from_block].y = y_to;
        block[from_block].z = z_to;
    }

    /* Now update the cost function.  May have to do major optimizations *
     * here later.                                                       */
    int num_of_pins = block[from_block].type->num_pins;
    int num_nets_affected = find_affected_nets(nets_to_update,
                                               net_block_moved,
                                               from_block,
                                               to_block,
                                               num_of_pins);
    if (place_cost_type == NONLINEAR_CONG) {
        save_region_occ(old_region_occ_x,
                        old_region_occ_y,
                        num_regions);
    }

    /* I'm using negative values of temp_net_cost as a flag, so DO NOT   *
     * use cost functions that can go negative.                          */
    double delta_c = 0.0;  /* Change in cost due to this swap. */
    double bb_delta_c = 0.0;
    int bb_index = 0; /* Index of new bounding box. */
    for (k = 0; k < num_nets_affected; ++k) {
        inet = nets_to_update[k];
        /* If we swapped two blocks connected to the same net, its bounding box *
         * doesn't change.                                                      */
        if (net_block_moved[k] == FROM_AND_TO) {
            continue;
        }

        if (net[inet].num_sinks < SMALL_NET) {
            get_non_updateable_bb(inet, &bb_coord_new[bb_index]);
        } else {
            if (net_block_moved[k] == FROM) {
                update_bb(inet, &bb_coord_new[bb_index],
                          &bb_edge_new[bb_index], x_from, y_from,
                          x_to, y_to);
            } else {
                update_bb(inet, &bb_coord_new[bb_index],
                          &bb_edge_new[bb_index], x_to, y_to, x_from,
                          y_from);
            }
        }

        if (place_cost_type != NONLINEAR_CONG) {
            temp_net_cost[inet] =
                get_net_cost(inet, &bb_coord_new[bb_index]);
            bb_delta_c += temp_net_cost[inet] - net_cost[inet];
        } else {
            /* Rip up, then replace with new bb. */
            update_region_occ(inet, &bb_coords[inet], -1,
                              num_regions);
            update_region_occ(inet, &bb_coord_new[bb_index], 1,
                              num_regions);
        }

        ++bb_index;
    }  /* end of for(k = 0; k < num_nets_affected; k++) */

    double newcost = 0.0;
    if (place_cost_type == NONLINEAR_CONG) {
        newcost = nonlinear_cong_cost(num_regions);
        bb_delta_c = newcost - *bb_cost;
    }

    double timing_delta_c = 0.0;
    double delay_delta_c = 0.0;
    if (place_algorithm == NET_TIMING_DRIVEN_PLACE ||
           place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        /*in this case we redefine delta_c as a combination of timing and bb.  *
         *additionally, we normalize all values, therefore delta_c is in       *
         *relation to 1*/
        comp_delta_td_cost(from_block,
                           to_block,
                           num_of_pins,
                           &timing_delta_c,
                           &delay_delta_c);

        delta_c = (1 - timing_tradeoff) * bb_delta_c * inverse_prev_bb_cost
                     + timing_tradeoff * timing_delta_c * inverse_prev_timing_cost;
    } else {
        delta_c = bb_delta_c;
    }

    /* 1->move accepted, 0->rejected. */
    int keep_switch = assess_swap(delta_c,
                                  t);
    if (keep_switch) {
        *cost = *cost + delta_c;
        *bb_cost = *bb_cost + bb_delta_c;

        if (place_algorithm == NET_TIMING_DRIVEN_PLACE ||
                place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
            /*update the point_to_point_timing_cost and point_to_point_delay_cost
             * values from the temporary values */
            *timing_cost = *timing_cost + timing_delta_c;
            *delay_cost = *delay_cost + delay_delta_c;
            update_td_cost(from_block,
                           to_block,
                           num_of_pins);
        }

        /* update net cost functions and reset flags. */
        bb_index = 0;
        for (k = 0; k < num_nets_affected; ++k) {
            inet = nets_to_update[k];
            /* If we swapped two blocks connected to the same net, its bounding box *
             * doesn't change.                                                      */
            if (net_block_moved[k] == FROM_AND_TO) {
                temp_net_cost[inet] = -1;
                continue;
            }

            bb_coords[inet] = bb_coord_new[bb_index];
            if (net[inet].num_sinks >= SMALL_NET) {
                bb_num_on_edges[inet] = bb_edge_new[bb_index];
            }

            ++bb_index;

            net_cost[inet] = temp_net_cost[inet];
            temp_net_cost[inet] = -1;
        }

        /* Update fb data structures since we kept the move. */
        /* Swap physical location */
        grid[x_to][y_to].blocks[z_to] = from_block;
        grid[x_from][y_from].blocks[z_from] = to_block;

        if (EMPTY == to_block) {
            /* Moved to an empty location */
            grid[x_to][y_to].usage++;
            grid[x_from][y_from].usage--;
        }
    } else {
        /* Move was rejected.  */
        /* Reset the net cost function flags first. */
        for (k = 0; k < num_nets_affected; ++k) {
            inet = nets_to_update[k];
            temp_net_cost[inet] = -1;
        }

        /* Restore the block data structures to their state before the move. */
        block[from_block].x = x_from;
        block[from_block].y = y_from;
        block[from_block].z = z_from;
        if (to_block != EMPTY) {
            block[to_block].x = x_to;
            block[to_block].y = y_to;
            block[to_block].z = z_to;
        }

        /* Restore the region occupancies to their state before the move. */
        if (place_cost_type == NONLINEAR_CONG) {
            restore_region_occ(old_region_occ_x,
                               old_region_occ_y,
                               num_regions);
        }
    }

    return  keep_switch;
}  /* end of static int try_swap(double t, ) */


static void
save_region_occ(double** old_region_occ_x,
                double** old_region_occ_y,
                int num_regions)
{
    /* Saves the old occupancies of the placement subregions in case the  *
     * current move is not accepted.  Used only for NONLINEAR_CONG.       */
    int i, j;
    for (i = 0; i < num_regions; i++) {
        for (j = 0; j < num_regions; j++) {
            old_region_occ_x[i][j] = place_region_x[i][j].occupancy;
            old_region_occ_y[i][j] = place_region_y[i][j].occupancy;
        }
    }
}


static void
restore_region_occ(double** old_region_occ_x,
                   double** old_region_occ_y,
                   int num_regions)
{
    /* Restores the old occupancies of the placement subregions when the  *
     * current move is not accepted.  Used only for NONLINEAR_CONG.       */
    int i, j;
    for (i = 0; i < num_regions; i++) {
        for (j = 0; j < num_regions; j++) {
            place_region_x[i][j].occupancy = old_region_occ_x[i][j];
            place_region_y[i][j].occupancy = old_region_occ_y[i][j];
        }
    }
}


/* Puts a list of all the nets connected to from_block and to_block into *
 * nets_to_update. Returns the number of affected nets. Net_block_moved *
 * is either FROM, TO or FROM_AND_TO -- the block connected to this net  *
 * that has moved.                                                       */
static int find_affected_nets(int* nets_to_update,
                              int* net_block_moved,
                              int from_block,
                              int to_block,
                              int num_of_pins)
{
    int inet, count;
    int affected_index = 0;
    int k = -1;
    for (k = 0; k < num_of_pins; ++k) {
        inet = block[from_block].nets[k];
        if (inet == OPEN || TRUE == net[inet].is_global) {
            continue;
        }
        /* This is here in case the same block connects to a net twice. */
        if (temp_net_cost[inet] > 0.) {
            continue;
        }

        nets_to_update[affected_index] = inet;
        net_block_moved[affected_index] = FROM;
        affected_index++;
        temp_net_cost[inet] = 1.;   /* Flag to say we've marked this net. */
    }

    if (to_block != EMPTY) {
        for (k = 0; k < num_of_pins; k++) {
            inet = block[to_block].nets[k];
            if (inet == OPEN) {
                continue;
            }

            if (net[inet].is_global) {
                continue;
            }

            if (temp_net_cost[inet] > 0.) {
                /* Net already marked. */
                for (count = 0; count < affected_index; count++) {
                    if (nets_to_update[count] == inet) {
                        if (net_block_moved[count] == FROM)
                            net_block_moved[count] =
                                FROM_AND_TO;

                        break;
                    }
                }

#ifdef DEBUG

                if (count > affected_index) {
                    printf
                    ("Error in find_affected_nets -- count = %d,"
                     " affected index = %d.\n", count,
                     affected_index);
                    exit(1);
                }

#endif
            }

            else {
                /* Net not marked yet. */

                nets_to_update[affected_index] = inet;
                net_block_moved[affected_index] = TO;
                affected_index++;
                temp_net_cost[inet] = 1.;   /* Flag means we've  marked net. */
            }
        }
    }

    return (affected_index);
}


static boolean
find_to(int x_from,
        int y_from,
        block_type_ptr type,
        double range_limit,
        int* x_lookup,
        int* x_to,
        int* y_to)
{

    /* Returns the point to which I want to swap, properly range limited.
     * range_limit must always be between 1 and num_grid_columns (inclusive) for this routine
     * to work.  Assumes that a column only contains blocks of the same type.
     */

    int x_rel, y_rel, iside, iplace, rlx, rly, min_x, max_x, min_y, max_y;
    int num_col_same_type, i, j;

    rlx = min(num_grid_columns, range_limit);    /* Only needed when num_grid_columns < num_grid_rows. */
    rly = min(num_grid_rows, range_limit);   /* Added rly for aspect_ratio != 1 case. */

    min_x = max(1, x_from - rlx);
    max_x = min(num_grid_columns, x_from + rlx);
    min_y = max(1, y_from - rly);
    max_y = min(num_grid_rows, y_from + rly);

    num_col_same_type = 0;
    j = 0;

    if (type != IO_TYPE) {
        for (i = min_x; i <= max_x; i++) {
            if (grid[i][1].type == type) {
                num_col_same_type++;
                x_lookup[j] = i;
                j++;
            }
        }

        assert(num_col_same_type != 0);

        if (num_col_same_type == 1 &&
                ((((max_y - min_y) / type->height) - 1) <= 0
                 || type->height > (num_grid_rows / 2))) {
            return FALSE;
        }
    }

#ifdef DEBUG

    if (rlx < 1 || rlx > num_grid_columns) {
        printf("Error in find_to: rlx = %d\n", rlx);
        exit(1);
    }

#endif

    do {
        /* Until (x_to, y_to) different from (x_from, y_from) */
        if (type == IO_TYPE) {
            /* io_block to be moved. */
            if (rlx >= num_grid_columns) {
                iside = my_irand(3);

                /*                              *
                 *       +-----1----+           *
                 *       |          |           *
                 *       |          |           *
                 *       0          2           *
                 *       |          |           *
                 *       |          |           *
                 *       +-----3----+           *
                 *                              */
                switch (iside) {
                    case 0:
                        iplace = my_irand(num_grid_rows - 1) + 1;
                        *x_to = 0;
                        *y_to = iplace;
                        break;

                    case 1:
                        iplace = my_irand(num_grid_columns - 1) + 1;
                        *x_to = iplace;
                        *y_to = num_grid_rows + 1;
                        break;

                    case 2:
                        iplace = my_irand(num_grid_rows - 1) + 1;
                        *x_to = num_grid_columns + 1;
                        *y_to = iplace;
                        break;

                    case 3:
                        iplace = my_irand(num_grid_columns - 1) + 1;
                        *x_to = iplace;
                        *y_to = 0;
                        break;

                    default:
                        printf
                        ("Error in find_to.  Unexpected io swap location.\n");
                        exit(1);
                }
            } else {
                /* rlx is less than whole chip */
                if (x_from == 0) {
                    iplace = my_irand(2 * rly);
                    *y_to = y_from - rly + iplace;
                    *x_to = x_from;

                    if (*y_to > num_grid_rows) {
                        *y_to = num_grid_rows + 1;
                        *x_to = my_irand(rlx - 1) + 1;
                    } else if (*y_to < 1) {
                        *y_to = 0;
                        *x_to = my_irand(rlx - 1) + 1;
                    }
                } else if (x_from == num_grid_columns + 1) {
                    iplace = my_irand(2 * rly);
                    *y_to = y_from - rly + iplace;
                    *x_to = x_from;

                    if (*y_to > num_grid_rows) {
                        *y_to = num_grid_rows + 1;
                        *x_to = num_grid_columns - my_irand(rlx - 1);
                    } else if (*y_to < 1) {
                        *y_to = 0;
                        *x_to = num_grid_columns - my_irand(rlx - 1);
                    }
                } else if (y_from == 0) {
                    iplace = my_irand(2 * rlx);
                    *x_to = x_from - rlx + iplace;
                    *y_to = y_from;

                    if (*x_to > num_grid_columns) {
                        *x_to = num_grid_columns + 1;
                        *y_to = my_irand(rly - 1) + 1;
                    } else if (*x_to < 1) {
                        *x_to = 0;
                        *y_to = my_irand(rly - 1) + 1;
                    }
                } else {
                    /* *y_from == num_grid_rows + 1 */
                    iplace = my_irand(2 * rlx);
                    *x_to = x_from - rlx + iplace;
                    *y_to = y_from;

                    if (*x_to > num_grid_columns) {
                        *x_to = num_grid_columns + 1;
                        *y_to = num_grid_rows - my_irand(rly - 1);
                    } else if (*x_to < 1) {
                        *x_to = 0;
                        *y_to = num_grid_rows - my_irand(rly - 1);
                    }
                }
            }   /* End rlx if */
        }       /* end type if */
        else {
            x_rel = my_irand(num_col_same_type - 1);
            y_rel =
                my_irand(max
                         (0, ((max_y - min_y) / type->height) - 1));
            *x_to = x_lookup[x_rel];
            *y_to = min_y + y_rel * type->height;
            *y_to = (*y_to) - grid[*x_to][*y_to].offset;    /* align it */
            assert(*x_to >= 1 && *x_to <= num_grid_columns);
            assert(*y_to >= 1 && *y_to <= num_grid_rows);
        }
    } while ((x_from == *x_to) && (y_from == *y_to));

#ifdef DEBUG

    if (*x_to < 0 || *x_to > num_grid_columns + 1 || *y_to < 0 || *y_to > num_grid_rows + 1) {
        printf("Error in routine find_to:  (x_to,y_to) = (%d,%d)\n",
               *x_to, *y_to);
        exit(1);
    }

#endif
    assert(type == grid[*x_to][*y_to].type);
    return TRUE;
}


static int
assess_swap(double delta_c,
            double t)
{

    /* Returns: 1->move accepted, 0->rejected. */

    int accept;
    double prob_fac, fnum;

    if (delta_c <= 0) {

#ifdef SPEC         /* Reduce variation in final solution due to round off */
        fnum = my_frand();
        //fnum = randfloat();
#endif

        accept = 1;
        return (accept);
    }

    if (t == 0.) {
        return (0);
    }

    fnum = my_frand();
    //fnum = randfloat();
    prob_fac = exp(-delta_c / t);

    if (prob_fac > fnum) {
        accept = 1;
    } else {
        accept = 0;
    }

    return (accept);
}

static double comp_td_point_to_point_delay(int inet,
                                          int ipin)
{
    /*returns the Tdel of one point to point connection */
    double delay_source_to_sink = 0.0;

    const int source_block = net[inet].node_block[0];
    const int sink_block = net[inet].node_block[ipin];
    const block_type_ptr sink_type = block[sink_block].type;
    const block_type_ptr source_type = block[source_block].type;
    assert(source_type != NULL && sink_type != NULL);

    const int delta_x = abs(block[sink_block].x - block[source_block].x);
    const int delta_y = abs(block[sink_block].y - block[source_block].y);

    /* TODO low priority: Could be merged into one look-up table */
    /* Note: This heuristic is terrible on Quality of Results.
     * A much better heuristic is to create a more comprehensive lookup table but
     * it's too late in the release cycle to do this.  Pushing until the next release */
    if (source_type == IO_TYPE) {
        if (sink_type == IO_TYPE) {
            delay_source_to_sink = delta_io_to_io[delta_x][delta_y];
        } else {
            delay_source_to_sink = delta_io_to_fb[delta_x][delta_y];
        }
    } else {
        if (sink_type == IO_TYPE) {
            delay_source_to_sink = delta_fb_to_io[delta_x][delta_y];
        } else {
            delay_source_to_sink = delta_fb_to_fb[delta_x][delta_y];
        }
    }

    if (delay_source_to_sink < 0) {
        printf("Error in comp_td_point_to_point_delay in place.c, bad delay_source_to_sink value\n");
        exit(1);
    }

    if (delay_source_to_sink < 0.) {
        printf("Error in comp_td_point_to_point_delay in place.c, Tdel is less than 0\n");
        exit(1);
    }

    return delay_source_to_sink;
}  /* end of static double comp_td_point_to_point_delay(int inet,) */


static void update_td_cost(int from_block,
                           int to_block,
                           int num_of_pins)
{
    /*update the point_to_point_timing_cost values from the temporary */
    /*values for all connections that have changed */
    int blkpin, net_pin, inet, ipin;
    for (blkpin = 0; blkpin < num_of_pins; blkpin++) {
        inet = block[from_block].nets[blkpin];
        if (inet == OPEN) {
            continue;
        }
        if (net[inet].is_global) {
            continue;
        }

        net_pin = net_pin_index[from_block][blkpin];
        if (net_pin != 0) {
            /*the following "if" prevents the value from being updated twice */
            if (net[inet].node_block[0] != to_block
                    && net[inet].node_block[0] != from_block) {

                point_to_point_delay_cost[inet][net_pin] =
                    temp_point_to_point_delay_cost[inet][net_pin];
                temp_point_to_point_delay_cost[inet][net_pin] =
                    -1;

                point_to_point_timing_cost[inet][net_pin] =
                    temp_point_to_point_timing_cost[inet]
                    [net_pin];
                temp_point_to_point_timing_cost[inet][net_pin] =
                    -1;
            }
        } else {
            /*this net is being driven by a moved block, recompute */
            /*all point to point connections on this net. */
            for (ipin = 1; ipin <= net[inet].num_sinks; ipin++) {

                point_to_point_delay_cost[inet][ipin] =
                    temp_point_to_point_delay_cost[inet][ipin];
                temp_point_to_point_delay_cost[inet][ipin] = -1;

                point_to_point_timing_cost[inet][ipin] =
                    temp_point_to_point_timing_cost[inet][ipin];
                temp_point_to_point_timing_cost[inet][ipin] = -1;
            }
        }
    }

    if (to_block != EMPTY) {
        for (blkpin = 0; blkpin < num_of_pins; blkpin++) {

            inet = block[to_block].nets[blkpin];

            if (inet == OPEN) {
                continue;
            }

            if (net[inet].is_global) {
                continue;
            }

            net_pin = net_pin_index[to_block][blkpin];

            if (net_pin != 0) {

                /*the following "if" prevents the value from being updated 2x */
                if (net[inet].node_block[0] != to_block
                        && net[inet].node_block[0] != from_block) {

                    point_to_point_delay_cost[inet][net_pin] =
                        temp_point_to_point_delay_cost[inet]
                        [net_pin];
                    temp_point_to_point_delay_cost[inet]
                    [net_pin] = -1;

                    point_to_point_timing_cost[inet][net_pin]
                        =
                            temp_point_to_point_timing_cost[inet]
                            [net_pin];
                    temp_point_to_point_timing_cost[inet]
                    [net_pin] = -1;
                }
            } else {
                /*this net is being driven by a moved block, recompute */
                /*all point to point connections on this net. */
                for (ipin = 1; ipin <= net[inet].num_sinks; ipin++) {

                    point_to_point_delay_cost[inet][ipin] =
                        temp_point_to_point_delay_cost[inet]
                        [ipin];
                    temp_point_to_point_delay_cost[inet][ipin]
                        = -1;

                    point_to_point_timing_cost[inet][ipin] =
                        temp_point_to_point_timing_cost[inet]
                        [ipin];
                    temp_point_to_point_timing_cost[inet]
                    [ipin] = -1;
                }
            }
        }
    }
}

static void comp_delta_td_cost(int from_block,
                               int to_block,
                               int num_of_pins,
                               double* delta_timing,
                               double* delta_delay)
{
    /*a net that is being driven by a moved block must have all of its  */
    /*sink timing costs recomputed. A net that is driving a moved block */
    /*must only have the timing cost on the connection driving the input_args */
    /*pin computed */
    int inet, k, ipin;

    double delta_timing_cost = 0.0;
    double delta_delay_cost = 0.0;
    for (k = 0; k < num_of_pins; ++k) {
        const int inet = block[from_block].nets[k];
        if (inet == OPEN || net[inet].is_global) {
            continue;
        }

        const int net_pin = net_pin_index[from_block][k];
        if (net_pin != 0) {
            /*this net is driving a moved block               */

            /*if this net is being driven by a block that has moved, we do not  */
            /*need to compute the change in the timing cost (here) since it will */
            /*be computed in the fanout of the net on  the driving block, also  */
            /*computing it here would double count the change, and mess up the  */
            /*delta_timing_cost value */
            if (net[inet].node_block[0] != to_block
                    && net[inet].node_block[0] != from_block) {
                double temp_delay = comp_td_point_to_point_delay(inet,
                                                                net_pin);

                temp_point_to_point_delay_cost[inet][net_pin] = temp_delay;
                temp_point_to_point_timing_cost[inet][net_pin] =
                        timing_place_crit[inet][net_pin] * temp_delay;

                delta_delay_cost +=
                    (temp_point_to_point_delay_cost[inet][net_pin]
                        - point_to_point_delay_cost[inet][net_pin]);

                delta_timing_cost +=
                    (temp_point_to_point_timing_cost[inet][net_pin]
                        - point_to_point_timing_cost[inet][net_pin]);
            }
        } else {
            /*this net is being driven by a moved block, recompute */
            /*all point to point connections on this net. */
            for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                double temp_delay = comp_td_point_to_point_delay(inet,
                                                                ipin);
                temp_point_to_point_delay_cost[inet][ipin] = temp_delay;
                temp_point_to_point_timing_cost[inet][ipin] =
                        timing_place_crit[inet][ipin] * temp_delay;

                delta_delay_cost +=
                    temp_point_to_point_delay_cost[inet][ipin]
                      - point_to_point_delay_cost[inet][ipin];

                delta_timing_cost +=
                    temp_point_to_point_timing_cost[inet][ipin]
                      - point_to_point_timing_cost[inet][ipin];
            }  /* end of for(ipin = 1; ipin <= net[inet].num_sinks; ++ipin) */
        }  /* end of else(net_pin == 0) */
    }  /* end of for (k = 0; k < num_of_pins; ++k) */

    if (to_block != EMPTY) {
        for (k = 0; k < num_of_pins; ++k) {
            inet = block[to_block].nets[k];
            if (inet == OPEN || net[inet].is_global) {
                continue;
            }

            const int net_pin = net_pin_index[to_block][k];
            if (net_pin != 0) {
                /*this net is driving a moved block */

                /*if this net is being driven by a block that has moved, we do not */
                /*need to compute the change in the timing cost (here) since it was */
                /*computed in the fanout of the net on  the driving block, also    */
                /*computing it here would double count the change, and mess up the */
                /*delta_timing_cost value */
                if (net[inet].node_block[0] != to_block
                        && net[inet].node_block[0] != from_block) {
                    double temp_delay = comp_td_point_to_point_delay(inet,
                                                                    net_pin);
                    temp_point_to_point_delay_cost[inet][net_pin] = temp_delay;
                    temp_point_to_point_timing_cost[inet][net_pin] =
                        timing_place_crit[inet][net_pin] * temp_delay;

                    delta_delay_cost +=
                        temp_point_to_point_delay_cost[inet][net_pin]
                          - point_to_point_delay_cost[inet][net_pin];
                    delta_timing_cost +=
                        temp_point_to_point_timing_cost[inet][net_pin]
                          - point_to_point_timing_cost[inet][net_pin];
                }
            } else {
                /*this net is being driven by a moved block, recompute */
                /*all point to point connections on this net. */
                for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                    double temp_delay = comp_td_point_to_point_delay(inet,
                                                                    ipin);
                    temp_point_to_point_delay_cost[inet][ipin] = temp_delay;
                    temp_point_to_point_timing_cost[inet][ipin] =
                        timing_place_crit[inet][ipin] * temp_delay;

                    delta_delay_cost +=
                        temp_point_to_point_delay_cost[inet][ipin]
                          - point_to_point_delay_cost[inet][ipin];
                    delta_timing_cost +=
                        temp_point_to_point_timing_cost[inet][ipin]
                          - point_to_point_timing_cost[inet][ipin];
                }
            }  /* end of else(net_pin == 0), this net_pin was a driver_pin*/
        } /* end of for(k = 0; k < num_of_pins; ++k) */
    }  /* end of if(to_block != EMPTY) */

    *delta_timing = delta_timing_cost;
    *delta_delay = delta_delay_cost;
}

/*computes the cost (from scratch) due to the delays and criticalities*
 *on all point to point connections, we define the timing cost of     *
 *each connection as criticality*Tdel */
static void comp_td_costs(double* timing_cost,
                          double* connection_delay_sum)
{
    double loc_timing_cost = 0.0;
    double loc_connection_delay_sum = 0.0;
    int inet = 0;
    for (inet = 0; inet < num_nets; ++inet) {
        /* for each net ... */
        if (net[inet].is_global == FALSE) {
            /* Do only if not global. */
            int ipin = 0;
            for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                double temp_delay_cost = comp_td_point_to_point_delay(inet,
                                                                     ipin);
                double temp_timing_cost =
                         temp_delay_cost * timing_place_crit[inet][ipin];

                loc_connection_delay_sum += temp_delay_cost;
                point_to_point_delay_cost[inet][ipin] = temp_delay_cost;
                temp_point_to_point_delay_cost[inet][ipin] = -1; /*undefined */

                point_to_point_timing_cost[inet][ipin] = temp_timing_cost;
                temp_point_to_point_timing_cost[inet][ipin] = -1;   /*undefined */
                loc_timing_cost += temp_timing_cost;
            }
        }
    }

    *timing_cost = loc_timing_cost;
    *connection_delay_sum = loc_connection_delay_sum;
}  /* end of static void comp_td_costs(double* timing_cost, ) */

/* Finds the cost from scratch.  Done only when the placement   *
 * has been radically changed (i.e. after initial placement).   *
 * Otherwise find the cost change incrementally.  If method     *
 * check is NORMAL, we find bounding boxes that are updateable  *
 * for the larger nets.  If method is CHECK, all bounding boxes *
 * are found via the non_updateable_bb routine, to provide a    *
 * cost which can be used to check the correctness of the       *
 * other routine.                                               */
static double comp_bb_cost(int method,
                           int place_cost_type,
                           int num_regions)
{
    double cost = 0;
    double expected_wirelength = 0.0;

    /* Initialize occupancies to zero if regions are being used. */
    int i, j, k;
    if (place_cost_type == NONLINEAR_CONG) {
        for (i = 0; i < num_regions; i++) {
            for (j = 0; j < num_regions; j++) {
                place_region_x[i][j].occupancy = 0.;
                place_region_y[i][j].occupancy = 0.;
            }
        }
    }

    for (k = 0; k < num_nets; k++) {
        /* for each net ... */
        if (net[k].is_global == FALSE) {
            /* Do only if not global. */

            /* Small nets don't use incremental updating on their bounding boxes, *
             * so they can use a fast bounding box calculator.                    */

            if (net[k].num_sinks >= SMALL_NET && method == NORMAL) {
                get_bb_from_scratch(k, &bb_coords[k],
                                    &bb_num_on_edges[k]);
            } else {
                get_non_updateable_bb(k, &bb_coords[k]);
            }

            if (place_cost_type != NONLINEAR_CONG) {
                net_cost[k] = get_net_cost(k, &bb_coords[k]);
                cost += net_cost[k];

                if (method == CHECK)
                    expected_wirelength +=
                        get_net_wirelength_estimate(k,
                                                    &bb_coords
                                                    [k]);
            } else {
                /* Must be nonlinear_cong case. */
                update_region_occ(k, &bb_coords[k], 1,
                                  num_regions);
            }
        }
    }

    if (place_cost_type == NONLINEAR_CONG) {
        cost = nonlinear_cong_cost(num_regions);
    }

    //printf("BB estimate of min-dist (placement) wirelength is ;%.0f\n",expected_wirelength);
    if (method == CHECK)
    {
        return (cost);
    }
}

/* This routine computes the cost of a placement when the NONLINEAR_CONG *
 * option is selected.  It assumes that the occupancies of all the       *
 * placement subregions have been properly updated, and simply           *
 * computes the cost due to these occupancies by summing over all        *
 * subregions.  This will be inefficient for moves that don't affect     *
 * many subregions (i.e. small moves late in placement), esp. when there *
 * are a lot of subregions.  May recode later to update only affected    *
 * subregions.                                                           */
static double nonlinear_cong_cost(int num_regions)
{
    double cost = 0.0;
    double tmp = 0.0;
    int i, j;
    for (i = 0; i < num_regions; i++) {
        for (j = 0; j < num_regions; j++) {

            /* Many different cost metrics possible.  1st try:  */
            if (place_region_x[i][j].occupancy < place_region_x[i][j].capacity) {
                cost += place_region_x[i][j].occupancy *
                        place_region_x[i][j].inv_capacity;
            } else {
                /* Overused region -- penalize. */
                tmp = place_region_x[i][j].occupancy *
                      place_region_x[i][j].inv_capacity;
                cost += tmp * tmp;
            }

            if (place_region_y[i][j].occupancy <
                    place_region_y[i][j].capacity) {
                cost += place_region_y[i][j].occupancy *
                        place_region_y[i][j].inv_capacity;
            } else {
                /* Overused region -- penalize. */

                tmp = place_region_y[i][j].occupancy *
                      place_region_y[i][j].inv_capacity;
                cost += tmp * tmp;
            }

        }
    }

    return (cost);
}


/* Called only when the place_cost_type is NONLINEAR_CONG.  If add_or_sub *
 * is 1, this uses the new net bounding box to increase the occupancy     *
 * of some regions.  If add_or_sub = - 1, it decreases the occupancy      *
 * by that due to this bounding box.                                      */
static void update_region_occ(int inet,
                              bbox_t* coords,
                              int add_or_sub,
                              int num_regions)
{

    double net_xmin, net_xmax, net_ymin, net_ymax, crossing;
    double inv_region_len, inv_region_height;
    double inv_bb_len, inv_bb_height;
    double overlap_xlow, overlap_xhigh, overlap_ylow, overlap_yhigh;
    double y_overlap, x_overlap, x_occupancy, y_occupancy;
    int imin, imax, jmin, jmax, i, j;

    if (net[inet].num_sinks >= 50) {
        crossing = 2.7933 + 0.02616 * ((net[inet].num_sinks + 1) - 50);
    } else {
        crossing = cross_count[net[inet].num_sinks];
    }

    net_xmin = coords->xmin - 0.5;
    net_xmax = coords->xmax + 0.5;
    net_ymin = coords->ymin - 0.5;
    net_ymax = coords->ymax + 0.5;

    /* I could precompute the two values below.  Should consider this. */

    inv_region_len = (double)num_regions / (double)num_grid_columns;
    inv_region_height = (double)num_regions / (double)num_grid_rows;

    /* Get integer coordinates defining the rectangular area in which the *
     * subregions have to be updated.  Formula is as follows:  subtract   *
     * 0.5 from net_xmin, etc. to get numbers from 0 to num_grid_columns or num_grid_rows;         *
     * divide by num_grid_columns or num_grid_rows to scale between 0 and 1; multiply by           *
     * num_regions to scale between 0 and num_regions; and truncate to    *
     * get the final answer.                                              */

    imin = (int)(net_xmin - 0.5) * inv_region_len;
    imax = (int)(net_xmax - 0.5) * inv_region_len;
    imax = min(imax, num_regions - 1);  /* Watch for weird roundoff */

    jmin = (int)(net_ymin - 0.5) * inv_region_height;
    jmax = (int)(net_ymax - 0.5) * inv_region_height;
    jmax = min(jmax, num_regions - 1);  /* Watch for weird roundoff */

    inv_bb_len = 1. / (net_xmax - net_xmin);
    inv_bb_height = 1. / (net_ymax - net_ymin);

    /* See RISA paper (ICCAD '94, pp. 690 - 695) for a description of why *
     * I use exactly this cost function.                                  */

    for (i = imin; i <= imax; i++) {
        for (j = jmin; j <= jmax; j++) {
            overlap_xlow = max(place_region_bounds_x[i], net_xmin);
            overlap_xhigh =
                min(place_region_bounds_x[i + 1], net_xmax);
            overlap_ylow = max(place_region_bounds_y[j], net_ymin);
            overlap_yhigh =
                min(place_region_bounds_y[j + 1], net_ymax);

            x_overlap = overlap_xhigh - overlap_xlow;
            y_overlap = overlap_yhigh - overlap_ylow;

#ifdef DEBUG

            if (x_overlap < -0.001) {
                printf
                ("Error in update_region_occ:  x_overlap < 0"
                 "\n inet = %d, overlap = %g\n", inet,
                 x_overlap);
            }

            if (y_overlap < -0.001) {
                printf
                ("Error in update_region_occ:  y_overlap < 0"
                 "\n inet = %d, overlap = %g\n", inet,
                 y_overlap);
            }

#endif


            x_occupancy =
                crossing * y_overlap * x_overlap * inv_bb_height *
                inv_region_len;
            y_occupancy =
                crossing * x_overlap * y_overlap * inv_bb_len *
                inv_region_height;

            place_region_x[i][j].occupancy +=
                add_or_sub * x_occupancy;
            place_region_y[i][j].occupancy +=
                add_or_sub * y_occupancy;
        }
    }
}


static void free_place_regions(int num_regions)
{
    /* Frees the place_regions data structures needed by the NONLINEAR_CONG *
     * cost function.                                                       */
    free_matrix(place_region_x, 0, num_regions - 1,
                0, sizeof(place_region_t));

    free_matrix(place_region_y, 0, num_regions - 1,
                0, sizeof(place_region_t));

    free(place_region_bounds_x);
    free(place_region_bounds_y);
}


static void free_placement_structs(int place_cost_type,
                                   int num_regions,
                                   double** old_region_occ_x,
                                   double** old_region_occ_y,
                                   placer_opts_t placer_opts)
{
    /* Frees the major structures needed by the placer (and not needed       *
     * elsewhere).   */
    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE ||
            placer_opts.enable_timing_computations) {
        int inet = -1;
        for (inet = 0; inet < num_nets; ++inet) {
            /*add one to the address since it is indexed from 1 not 0 */

            point_to_point_delay_cost[inet]++;
            free(point_to_point_delay_cost[inet]);

            point_to_point_timing_cost[inet]++;
            free(point_to_point_timing_cost[inet]);

            temp_point_to_point_delay_cost[inet]++;
            free(temp_point_to_point_delay_cost[inet]);

            temp_point_to_point_timing_cost[inet]++;
            free(temp_point_to_point_timing_cost[inet]);
        }

        free(point_to_point_delay_cost);
        free(temp_point_to_point_delay_cost);

        free(point_to_point_timing_cost);
        free(temp_point_to_point_timing_cost);

        free_matrix(net_pin_index,
                    0,
                    num_blocks - 1,
                    0,
                    sizeof(int));
    }


    free(net_cost);
    free(temp_net_cost);
    free(bb_num_on_edges);
    free(bb_coords);

    net_cost = NULL;        /* Defensive coding. */
    temp_net_cost = NULL;
    bb_num_on_edges = NULL;
    bb_coords = NULL;

    free_unique_pin_list();

    if (place_cost_type == NONLINEAR_CONG) {
        free_place_regions(num_regions);
        free_matrix(old_region_occ_x, 0, num_regions - 1, 0,
                    sizeof(double));
        free_matrix(old_region_occ_y, 0, num_regions - 1, 0,
                    sizeof(double));
    } else if (place_cost_type == LINEAR_CONG) {
        free_fast_cost_update_structs();
    } else {
        /* NO OPERATION. */
    }
}  /* end of static void free_placement_structs(int place_cost_type,) */


static void alloc_and_load_placement_structs(int place_cost_type,
                                             int num_regions,
                                             double place_cost_exp,
                                             double** *old_region_occ_x,
                                             double** *old_region_occ_y,
                                             placer_opts_t placer_opts)
{
    /* Allocates the major structures needed only by the placer, primarily for *
     * computing costs quickly and such.                                       */
    int inet, ipin, i;
    int max_pins_per_fb = 0;
    for (i = 0; i < num_types; ++i) {
        max_pins_per_fb = max(max_pins_per_fb,
                              type_descriptors[i].num_pins);
    }

    /* Allocate the structures associated with timing-driven_placement   */
    /* [0..num_nets-1][1..num_pins-1], pin 0 was driver pin of net. */
    if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE ||
            placer_opts.enable_timing_computations) {
        point_to_point_delay_cost =
                         (double**)my_malloc(num_nets * sizeof(double*));
        temp_point_to_point_delay_cost =
                         (double**)my_malloc(num_nets * sizeof(double*));

        point_to_point_timing_cost =
                         (double**)my_malloc(num_nets * sizeof(double*));
        temp_point_to_point_timing_cost =
                         (double**)my_malloc(num_nets * sizeof(double*));

        net_pin_index = (int**)alloc_matrix(0, num_blocks - 1,
                                            0, max_pins_per_fb - 1,
                                            sizeof(int));
        for (inet = 0; inet < num_nets; ++inet) {
            /* in the following, subract one so index starts at *
             * 1 instead of 0 */
            point_to_point_delay_cost[inet] = 
                    (double*)my_malloc(net[inet].num_sinks * sizeof(double));
            --(point_to_point_delay_cost[inet]);

            temp_point_to_point_delay_cost[inet] =
                (double*)my_malloc(net[inet].num_sinks * sizeof(double));
            --(temp_point_to_point_delay_cost[inet]);

            point_to_point_timing_cost[inet] =
                (double*)my_malloc(net[inet].num_sinks * sizeof(double));
            --(point_to_point_timing_cost[inet]);

            temp_point_to_point_timing_cost[inet] =
                (double*)my_malloc(net[inet].num_sinks * sizeof(double));
            --(temp_point_to_point_timing_cost[inet]);
        }

        for (inet = 0; inet < num_nets; ++inet) {
            for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                point_to_point_delay_cost[inet][ipin] = 0;
                temp_point_to_point_delay_cost[inet][ipin] = 0;
            }
        }
    }

    net_cost = (double*)my_malloc(num_nets * sizeof(double));
    temp_net_cost = (double*)my_malloc(num_nets * sizeof(double));
    /* Used to store costs for moves not yet made and to indicate when a net's   *
     * cost has been recomputed. temp_net_cost[inet] < 0 means net's cost hasn't *
     * been recomputed.                                                        */
    for (inet = 0; inet < num_nets; ++inet) {
        temp_net_cost[inet] = -1.0;
    }

    bb_coords = (bbox_t*)my_malloc(num_nets * sizeof(bbox_t));
    bb_num_on_edges = (bbox_t*)my_malloc(num_nets * sizeof(bbox_t));

    /* Get a list of pins with no duplicates. */
    alloc_and_load_unique_pin_list();

    /* Allocate storage for subregion data, if needed. */
    if (place_cost_type == NONLINEAR_CONG) {
        alloc_place_regions(num_regions);
        load_place_regions(num_regions);
        *old_region_occ_x = (double**)alloc_matrix(0, num_regions - 1,
                                                  0, num_regions - 1,
                                                  sizeof(double));
        *old_region_occ_y = (double**)alloc_matrix(0, num_regions - 1,
                                                  0, num_regions - 1,
                                                  sizeof(double));
    } else {
        /* Shouldn't use them; crash hard if I do!   */
        *old_region_occ_x = NULL;
        *old_region_occ_y = NULL;
    }

    if (place_cost_type == LINEAR_CONG) {
        alloc_and_load_for_fast_cost_update(place_cost_exp);
    }
}

/* Allocates memory for the regional occupancy, cost, etc. counts *
 * kept when we're using the NONLINEAR_CONG placement cost        *
 * function.                                                      */
static void alloc_place_regions(int num_regions)
{
    place_region_x =
        (place_region_t**)alloc_matrix(0, num_regions - 1, 0,
                                       num_regions - 1,
                                       sizeof(place_region_t));

    place_region_y =
        (place_region_t**)alloc_matrix(0, num_regions - 1, 0,
                                       num_regions - 1,
                                       sizeof(place_region_t));

    place_region_bounds_x = (double*)my_malloc((num_regions + 1) *
                                              sizeof(double));

    place_region_bounds_y = (double*)my_malloc((num_regions + 1) *
                                              sizeof(double));
}

/* Loads the capacity values in each direction for each of the placement *
 * regions.  The chip is divided into a num_regions x num_regions array. */
static void load_place_regions(int num_regions)
{

    int i, j, low_block, high_block, rnum;
    double low_lim, high_lim, capacity, fac, block_capacity;
    double len_fac, height_fac;

    /* First load up horizontal channel capacities.  */
    for (j = 0; j < num_regions; j++) {
        capacity = 0.;
        low_lim = (double)j / (double)num_regions * num_grid_rows + 1.;
        high_lim = (double)(j + 1) / (double)num_regions * num_grid_rows;

        low_block = floor(low_lim);
        low_block = max(1, low_block);  /* Watch for weird roundoff effects. */
        high_block = ceil(high_lim);
        high_block = min(high_block, num_grid_rows);

        block_capacity = (chan_width_x[low_block - 1] +
                          chan_width_x[low_block]) / 2.;

        if (low_block == 1) {
            block_capacity += chan_width_x[0] / 2.;
        }

        fac = 1. - (low_lim - low_block);
        capacity += fac * block_capacity;

        for (rnum = low_block + 1; rnum < high_block; rnum++) {
            block_capacity =
                (chan_width_x[rnum - 1] + chan_width_x[rnum]) / 2.;
            capacity += block_capacity;
        }

        block_capacity = (chan_width_x[high_block - 1] +
                          chan_width_x[high_block]) / 2.;

        if (high_block == num_grid_rows) {
            block_capacity += chan_width_x[num_grid_rows] / 2.;
        }

        fac = 1. - (high_block - high_lim);
        capacity += fac * block_capacity;

        for (i = 0; i < num_regions; i++) {
            place_region_x[i][j].capacity = capacity;
            place_region_x[i][j].inv_capacity = 1. / capacity;
            place_region_x[i][j].occupancy = 0.;
            place_region_x[i][j].cost = 0.;
        }
    }

    /* Now load vertical channel capacities.  */

    for (i = 0; i < num_regions; i++) {
        capacity = 0.;
        low_lim = (double)i / (double)num_regions * num_grid_columns + 1.;
        high_lim = (double)(i + 1) / (double)num_regions * num_grid_columns;

        low_block = floor(low_lim);
        low_block = max(1, low_block);  /* Watch for weird roundoff effects. */
        high_block = ceil(high_lim);
        high_block = min(high_block, num_grid_columns);

        block_capacity = (chan_width_y[low_block - 1] +
                          chan_width_y[low_block]) / 2.;

        if (low_block == 1) {
            block_capacity += chan_width_y[0] / 2.;
        }

        fac = 1. - (low_lim - low_block);
        capacity += fac * block_capacity;

        for (rnum = low_block + 1; rnum < high_block; rnum++) {
            block_capacity =
                (chan_width_y[rnum - 1] + chan_width_y[rnum]) / 2.;
            capacity += block_capacity;
        }

        block_capacity = (chan_width_y[high_block - 1] +
                          chan_width_y[high_block]) / 2.;

        if (high_block == num_grid_columns) {
            block_capacity += chan_width_y[num_grid_columns] / 2.;
        }

        fac = 1. - (high_block - high_lim);
        capacity += fac * block_capacity;

        for (j = 0; j < num_regions; j++) {
            place_region_y[i][j].capacity = capacity;
            place_region_y[i][j].inv_capacity = 1. / capacity;
            place_region_y[i][j].occupancy = 0.;
            place_region_y[i][j].cost = 0.;
        }
    }

    /* Finally set up the arrays indicating the limits of each of the *
     * placement subregions.                                          */

    len_fac = (double)num_grid_columns / (double)num_regions;
    height_fac = (double)num_grid_rows / (double)num_regions;

    place_region_bounds_x[0] = 0.5;
    place_region_bounds_y[0] = 0.5;

    for (i = 1; i <= num_regions; i++) {
        place_region_bounds_x[i] = place_region_bounds_x[i - 1] + len_fac;
        place_region_bounds_y[i] =
            place_region_bounds_y[i - 1] + height_fac;
    }
}


/* Frees the unique pin list structures.                               */
static void free_unique_pin_list(void)
{
    int any_dup, inet;
    any_dup = 0;

    for (inet = 0; inet < num_nets; inet++) {
        if (duplicate_pins[inet] != 0) {
            free(unique_pin_list[inet]);
            any_dup = 1;
        }
    }

    if (any_dup != 0) {
        free(unique_pin_list);
    }

    free(duplicate_pins);
}


/* This routine looks for (multiple pins going to the same block) in the *
 * pinlist of each net. If it finds any, it marks that net as having  *
 * duplicate pins, and creates a new pinlist with no duplicates. This *
 * is then used by the updatable bounding box calculation routine for  *
 * efficiency.                                                         */
static void alloc_and_load_unique_pin_list(void)
{
    duplicate_pins = (int*)my_calloc(num_nets,
                                     sizeof(int));
    /* [0..num_blocks-1]: number of times a block is listed in the pinlist *
     * of a net. Temp. storage. */
    int* times_listed = (int*)my_calloc(num_blocks,
                                        sizeof(int));
    int any_dups = 0;

    int inet = -1;
    for (inet = 0; inet < num_nets; ++inet) {
        int num_dup = 0;

        int ipin = -1;
        for (ipin = 0; ipin <= net[inet].num_sinks; ++ipin) {
            /* For a block, it may had more than 1 pin connect to a net. */
            int block_num = net[inet].node_block[ipin];
            ++times_listed[block_num];
            if (times_listed[block_num] > 1) {
                ++num_dup;
            }
        }

        if (num_dup > 0) {
            /* Duplicates found. Make unique pin list. */
            duplicate_pins[inet] = num_dup;
            if (any_dups == 0) {
                /* This is the first duplicate found */
                unique_pin_list = (int**)my_calloc(num_nets,
                                                   sizeof(int*));
                any_dups = 1;
            }

            unique_pin_list[inet] = (int*)my_malloc((net[inet].num_sinks + 1
                                                      - num_dup) * sizeof(int));
            int offset = 0;
            for (ipin = 0; ipin <= net[inet].num_sinks; ++ipin) {
                int block_num = net[inet].node_block[ipin];
                /* If a block pin had added to unique_pin_list, then set it 0.
                 * Otherwise, it must add more than 1 time. */
                if (times_listed[block_num] != 0) {
                    times_listed[block_num] = 0;
                    unique_pin_list[inet][offset] = block_num;
                    ++offset;
                }
            }
        } else {
            /* No duplicates found. Reset times_listed. */
            for (ipin = 0; ipin <= net[inet].num_sinks; ++ipin) {
                int block_num = net[inet].node_block[ipin];
                times_listed[block_num] = 0;
            }
        }
    } /* end of for(inet = 0; inet < num_nets; ++inet) */

    free((void*)times_listed);
} /* end of static void alloc_and_load_unique_pin_list(void) */

/* This routine finds the bounding box of each net from scratch (i.e.    *
 * from only the block location information).  It updates both the       *
 * coordinate and number of blocks on each tedge information.  It        *
 * should only be called when the bounding box information is not valid. */
static void
get_bb_from_scratch(int inet,
                    bbox_t* coords,
                    bbox_t* num_on_edges)
{
    int ipin, block_num, x, y, xmin, xmax, ymin, ymax;
    int xmin_edge, xmax_edge, ymin_edge, ymax_edge;
    int n_pins;
    int* plist;

    /* I need a list of blocks to which this net connects, with no block listed *
     * more than once, in order to get a proper count of the number on the tedge *
     * of the bounding box.                                                     */
    if (duplicate_pins[inet] == 0) {
        plist = net[inet].node_block;
        n_pins = net[inet].num_sinks + 1;
    } else {
        plist = unique_pin_list[inet];
        n_pins = (net[inet].num_sinks + 1) - duplicate_pins[inet];
    }

    x = block[plist[0]].x;
    y = block[plist[0]].y;

    x = max(min(x, num_grid_columns), 1);
    y = max(min(y, num_grid_rows), 1);

    xmin = x;
    ymin = y;
    xmax = x;
    ymax = y;
    xmin_edge = 1;
    ymin_edge = 1;
    xmax_edge = 1;
    ymax_edge = 1;

    for (ipin = 1; ipin < n_pins; ipin++) {

        block_num = plist[ipin];
        x = block[block_num].x;
        y = block[block_num].y;

        /* Code below counts IO blocks as being within the 1..num_grid_columns, 1..num_grid_rows clb array. *
         * This is because channels do not go out of the 0..num_grid_columns, 0..num_grid_rows range, and   *
         * I always take all channels impinging on the bounding box to be within   *
         * that bounding box.  Hence, this "movement" of IO blocks does not affect *
         * the which channels are included within the bounding box, and it         *
         * simplifies the code a lot.                                              */
        x = max(min(x, num_grid_columns), 1);
        y = max(min(y, num_grid_rows), 1);

        if (x == xmin) {
            xmin_edge++;
        }

        if (x == xmax) {
            /* Recall that xmin could equal xmax -- don't use else */
            xmax_edge++;
        } else if (x < xmin) {
            xmin = x;
            xmin_edge = 1;
        } else if (x > xmax) {
            xmax = x;
            xmax_edge = 1;
        }

        if (y == ymin) {
            ymin_edge++;
        }

        if (y == ymax) {
            ymax_edge++;
        } else if (y < ymin) {
            ymin = y;
            ymin_edge = 1;
        } else if (y > ymax) {
            ymax = y;
            ymax_edge = 1;
        }
    }

    /* Copy the coordinates and number on edges information into the proper   *
     * structures.                                                            */

    coords->xmin = xmin;
    coords->xmax = xmax;
    coords->ymin = ymin;
    coords->ymax = ymax;

    num_on_edges->xmin = xmin_edge;
    num_on_edges->xmax = xmax_edge;
    num_on_edges->ymin = ymin_edge;
    num_on_edges->ymax = ymax_edge;
}


/* WMF: Finds the estimate of wirelength due to one net by looking at   *
 * its coordinate bounding box.                                         */
static double
get_net_wirelength_estimate(int inet,
                            bbox_t* bbptr)
{

    double ncost, crossing;
    /* Get the expected "crossing count" of a net, based on its number *
     * of pins.  Extrapolate for very large nets.                      */

    if (((net[inet].num_sinks + 1) > 50) && ((net[inet].num_sinks + 1) < 85)) {
        crossing = 2.7933 + 0.02616 * ((net[inet].num_sinks + 1) - 50);
    } else if ((net[inet].num_sinks + 1) >= 85) {
        crossing =
            2.7933 + 0.011 * (net[inet].num_sinks + 1) -
            0.0000018 * (net[inet].num_sinks + 1) * (net[inet].num_sinks +
                                                     1);
    } else {
        crossing = cross_count[(net[inet].num_sinks + 1) - 1];
    }

    /* Could insert a check for xmin == xmax.  In that case, assume  *
     * connection will be made with no bends and hence no x-cost.    *
     * Same thing for y-cost.                                        */

    /* Cost = wire length along channel * cross_count / average      *
     * channel capacity.   Do this for x, then y direction and add.  */

    ncost = (bbptr->xmax - bbptr->xmin + 1) * crossing;

    ncost += (bbptr->ymax - bbptr->ymin + 1) * crossing;

    return (ncost);
}


static double get_net_cost(int inet,
                          bbox_t* bbox_ptr)
{
    /* Finds the cost due to one net by looking at its coordinate bounding  *
     * box.                                                                 */

    /* Get the expected "crossing count" of a net, based on its number *
     * of pins.  Extrapolate for very large nets.                      */
    double crossing = 0.0;
    if ((net[inet].num_sinks + 1) > 50) {
        /* crossing = 3.0; Old value  */
        crossing = 2.7933 + 0.02616 * ((net[inet].num_sinks + 1) - 50);
    } else {
        crossing = cross_count[(net[inet].num_sinks + 1) - 1];
    }

    /* Could insert a check for xmin == xmax.  In that case, assume  *
     * connection will be made with no bends and hence no x-cost.    *
     * Same thing for y-cost.                                        */

    /* Cost = wirelength along channel * cross_count / average_   *
     * channel_capacity. Do this for x, then y direction and add. */
    const int xmin = bbox_ptr->xmin;
    const int xmax = bbox_ptr->xmax;
    const int ymin = bbox_ptr->ymin;
    const int ymax = bbox_ptr->ymax;
    double net_cost = (xmax - xmin + 1) * crossing * chanx_place_cost_fac[ymax][ymin-1];
    net_cost += (ymax - ymin + 1) * crossing * chany_place_cost_fac[xmax][xmin-1];

    return net_cost;
}  /* end of static double get_net_cost(int inet, ) */

/* Finds the bounding-box of a net and stores its coordinates in the *
 * bb_coord_new data structure. This routine should only be called   *
 * for small nets, since it does not determine enough information for *
 * the bounding-box to be updated incrementally later.                *
 * Currently assumes channels on both sides of the CLBs forming the  *
 * edges of the bounding-box can be used. Essentially, I am assuming *
 * the pins always lie on the outside of the bounding box.           */
static void get_non_updateable_bb(int inet,
                                  bbox_t* bb_coord_new)
{
    /* Only for small nets */
    int x = block[net[inet].node_block[0]].x;
    int y = block[net[inet].node_block[0]].y;

    /* get the biggest bounding-box boundary among all nets of current block */
    int xmin = x;
    int ymin = y;
    int xmax = x;
    int ymax = y;
    int k = 0;
    for (k = 1; k < (net[inet].num_sinks + 1); ++k) {
        x = block[net[inet].node_block[k]].x;
        y = block[net[inet].node_block[k]].y;

        if (x < xmin) {
            xmin = x;
        } else if (x > xmax) {
            xmax = x;
        }

        if (y < ymin) {
            ymin = y;
        } else if (y > ymax) {
            ymax = y;
        }
    }  /* end of for(k = 1;... ) */

    /* Now I've found the coordinates of the bounding-box. There are no channels *
     * beyond num_grid_columns and num_grid_rows, so I want to clip to that. As well,*
     * since I'll always include the channel immediately below and the channel   *
     * immediately to the left of the bounding-box, I want to clip to 1 in both *
     * directions as well(since minimum channel index is 0). See route.c for a  *
     * channel diagram.                        */
    bb_coord_new->xmin = max(min(xmin, num_grid_columns), 1);
    bb_coord_new->ymin = max(min(ymin, num_grid_rows), 1);
    bb_coord_new->xmax = max(min(xmax, num_grid_columns), 1);
    bb_coord_new->ymax = max(min(ymax, num_grid_rows), 1);
}  /* end of static void get_non_updateable_bb(int inet,) */

/* Updates the bounding box of a net by storing its coordinates in    *
 * the bb_coord_new data structure and the number of blocks on each   *
 * tedge in the bb_edge_new data structure.  This routine should only *
 * be called for large nets, since it has some overhead relative to   *
 * just doing a brute force bounding box calculation.  The bounding   *
 * box coordinate and tedge information for inet must be valid before *
 * this routine is called.                                            *
 * Currently assumes channels on both sides of the CLBs forming the   *
 * edges of the bounding box can be used.  Essentially, I am assuming *
 * the pins always lie on the outside of the bounding box.            */
static void update_bb(int inet,
                      bbox_t* bb_coord_new,
                      bbox_t* bb_edge_new,
                      int xold, int yold,
                      int xnew, int ynew)
{
    /* IO blocks are considered to be one cell in for simplicity.*/
    xnew = max(min(xnew, num_grid_columns), 1);
    ynew = max(min(ynew, num_grid_rows), 1);
    xold = max(min(xold, num_grid_columns), 1);
    yold = max(min(yold, num_grid_rows), 1);

    /* Check if I can update the bounding box incrementally. */
    if (xnew < xold) {
        /*------    Move To Left.    ------*/
        /* Update the xmax fields for coordinates and number of edges first. */
        if (xold == bb_coords[inet].xmax) {
            /* Old position at xmax. */
            if (bb_num_on_edges[inet].xmax == 1) {
                get_bb_from_scratch(inet,
                                    bb_coord_new,
                                    bb_edge_new);
                return;
            } else {
                bb_edge_new->xmax = bb_num_on_edges[inet].xmax - 1;
                bb_coord_new->xmax = bb_coords[inet].xmax;
            }
        } else {
            /* Move to left, old postion was not at xmax. */
            bb_coord_new->xmax = bb_coords[inet].xmax;
            bb_edge_new->xmax = bb_num_on_edges[inet].xmax;
        }

        /* Now do the xmin fields for coordinates and number of edges. */
        if (xnew < bb_coords[inet].xmin) {
            /* Moved past xmin */
            bb_coord_new->xmin = xnew;
            bb_edge_new->xmin = 1;
        } else if (xnew == bb_coords[inet].xmin) {
            /* Moved to xmin */
            bb_coord_new->xmin = xnew;
            bb_edge_new->xmin = bb_num_on_edges[inet].xmin + 1;
        } else {
            /* Xmin unchanged. */
            bb_coord_new->xmin = bb_coords[inet].xmin;
            bb_edge_new->xmin = bb_num_on_edges[inet].xmin;
        }
    /* End of move to left case. */
    } else if (xnew > xold) {
        /* Move to right. */
        /* Update the xmin fields for coordinates and number of edges first. */
        if (xold == bb_coords[inet].xmin) {
            /* Old position at xmin. */
            if (bb_num_on_edges[inet].xmin == 1) {
                get_bb_from_scratch(inet,
                                    bb_coord_new,
                                    bb_edge_new);
                return;
            } else {
                bb_edge_new->xmin = bb_num_on_edges[inet].xmin - 1;
                bb_coord_new->xmin = bb_coords[inet].xmin;
            }
        } else {
            /* Move to right, old position was not at xmin. */
            bb_coord_new->xmin = bb_coords[inet].xmin;
            bb_edge_new->xmin = bb_num_on_edges[inet].xmin;
        }

        /* Now do the xmax fields for coordinates and number of edges. */
        if (xnew > bb_coords[inet].xmax) {
            /* Moved past xmax. */
            bb_coord_new->xmax = xnew;
            bb_edge_new->xmax = 1;
        } else if (xnew == bb_coords[inet].xmax) {
            /* Moved to xmax */
            bb_coord_new->xmax = xnew;
            bb_edge_new->xmax = bb_num_on_edges[inet].xmax + 1;
        } else {
            /* Xmax unchanged. */
            bb_coord_new->xmax = bb_coords[inet].xmax;
            bb_edge_new->xmax = bb_num_on_edges[inet].xmax;
        }
    /* End of move to right case. */
    } else {
        /* xnew == xold -- no x motion. */
        bb_coord_new->xmin = bb_coords[inet].xmin;
        bb_coord_new->xmax = bb_coords[inet].xmax;
        bb_edge_new->xmin = bb_num_on_edges[inet].xmin;
        bb_edge_new->xmax = bb_num_on_edges[inet].xmax;
    }

    /* Now account for the y-direction motion. It similar with x-direction. */
    if (ynew < yold) {
        /* Move down. */
        /* Update the ymax fields for coordinates and number of edges first. */
        if (yold == bb_coords[inet].ymax) {
            /* Old position at ymax. */
            if (bb_num_on_edges[inet].ymax == 1) {
                get_bb_from_scratch(inet,
                                    bb_coord_new,
                                    bb_edge_new);
                return;
            } else {
                bb_edge_new->ymax = bb_num_on_edges[inet].ymax - 1;
                bb_coord_new->ymax = bb_coords[inet].ymax;
            }
        } else {
            /* Move down, old postion was not at ymax. */
            bb_coord_new->ymax = bb_coords[inet].ymax;
            bb_edge_new->ymax = bb_num_on_edges[inet].ymax;
        }

        /* Now do the ymin fields for coordinates and number of edges. */
        if (ynew < bb_coords[inet].ymin) {
            /* Moved past ymin */
            bb_coord_new->ymin = ynew;
            bb_edge_new->ymin = 1;
        } else if (ynew == bb_coords[inet].ymin) {
            /* Moved to ymin */
            bb_coord_new->ymin = ynew;
            bb_edge_new->ymin = bb_num_on_edges[inet].ymin + 1;
        } else {
            /* ymin unchanged. */
            bb_coord_new->ymin = bb_coords[inet].ymin;
            bb_edge_new->ymin = bb_num_on_edges[inet].ymin;
        }
    /* End of move down case. */
    } else if (ynew > yold) {
        /* Moved up. */
        /* Update the ymin fields for coordinates and number of edges first. */
        if (yold == bb_coords[inet].ymin) {
            /* Old position at ymin. */
            if (bb_num_on_edges[inet].ymin == 1) {
                get_bb_from_scratch(inet,
                                    bb_coord_new,
                                    bb_edge_new);
                return;
            } else {
                bb_edge_new->ymin = bb_num_on_edges[inet].ymin - 1;
                bb_coord_new->ymin = bb_coords[inet].ymin;
            }
        } else {
            /* Moved up, old position was not at ymin. */
            bb_coord_new->ymin = bb_coords[inet].ymin;
            bb_edge_new->ymin = bb_num_on_edges[inet].ymin;
        }

        /* Now do the ymax fields for coordinates and number of edges. */
        if (ynew > bb_coords[inet].ymax) {
            /* Moved past ymax. */
            bb_coord_new->ymax = ynew;
            bb_edge_new->ymax = 1;
        } else if (ynew == bb_coords[inet].ymax) {
            /* Moved to ymax */
            bb_coord_new->ymax = ynew;
            bb_edge_new->ymax = bb_num_on_edges[inet].ymax + 1;
        } else {
            /* ymax unchanged. */
            bb_coord_new->ymax = bb_coords[inet].ymax;
            bb_edge_new->ymax = bb_num_on_edges[inet].ymax;
        }
    /* End of move up case. */
    } else {
        /* ynew == yold -- no y motion. */
        bb_coord_new->ymin = bb_coords[inet].ymin;
        bb_coord_new->ymax = bb_coords[inet].ymax;
        bb_edge_new->ymin = bb_num_on_edges[inet].ymin;
        bb_edge_new->ymax = bb_num_on_edges[inet].ymax;
    }
}  /* end of static void update_bb(int inet, ), Not parallel */


static void initial_placement(pad_loc_t pad_loc_type,
                              char* pad_loc_file)
{
    /* Randomly places the blocks to create an initial placement.     */
    struct s_pos {
        int x;
        int y;
        int z;
    }** pos;         /* [0..num_types-1][0..type_tsize - 1] */

    pos = (struct s_pos**)my_malloc(num_types * sizeof(struct s_pos*));
    int* count = (int*)my_calloc(num_types, sizeof(int));
    int* index = (int*)my_calloc(num_types, sizeof(int));

    /* Initialize all occupancy to zero. */
    int i, j, k, iblk, choice, type_index, x, y, z;
    for (i = 0; i <= num_grid_columns + 1; i++) {
        for (j = 0; j <= num_grid_rows + 1; j++) {
            grid[i][j].usage = 0;

            for (k = 0; k < grid[i][j].type->capacity; k++) {
                grid[i][j].blocks[k] = EMPTY;

                if (grid[i][j].offset == 0) {
                    ++count[grid[i][j].type->index];
                }
            }
        }
    }

    for (i = 0; i < num_types; i++) {
        pos[i] = (struct s_pos*)my_malloc(count[i] * sizeof(struct s_pos));
    }

    for (i = 0; i <= num_grid_columns + 1; ++i) {
        for (j = 0; j <= num_grid_rows + 1; ++j) {
            for (k = 0; k < grid[i][j].type->capacity; ++k) {
                if (grid[i][j].offset == 0) {
                    type_index = grid[i][j].type->index;
                    pos[type_index][index[type_index]].x = i;
                    pos[type_index][index[type_index]].y = j;
                    pos[type_index][index[type_index]].z = k;
                    ++index[type_index];
                }
            }
        }
    }

    for (iblk = 0; iblk < num_blocks; iblk++) {
        /* Don't do IOs if the user specifies IOs */
        if (!(block[iblk].type == IO_TYPE && pad_loc_type == USER)) {
            type_index = block[iblk].type->index;
            assert(count[type_index] > 0);
            /* choice >= 0 && choice <= count[type_index] -1 */
            choice = my_irand(count[type_index] - 1);
            x = pos[type_index][choice].x;
            y = pos[type_index][choice].y;
            z = pos[type_index][choice].z;
            grid[x][y].blocks[z] = iblk;
            ++grid[x][y].usage;

            /* Ensure randomizer doesn't pick this block again */
            pos[type_index][choice] = pos[type_index][count[type_index] - 1]; /* overwrite used block position */
            --count[type_index];
        }
    }

    if (pad_loc_type == USER) {
        read_user_pad_loc(pad_loc_file);
    }

    /* All the blocks are placed now.  Make the block array agree with the    *
     * clb array.                                                             */
    for (i = 0; i <= (num_grid_columns + 1); ++i) {
        for (j = 0; j <= (num_grid_rows + 1); ++j) {
            for (k = 0; k < grid[i][j].type->capacity; ++k) {
                assert(grid[i][j].blocks != NULL);

                iblk = grid[i][j].blocks[k];
                if (iblk != EMPTY) {
                    block[iblk].x = i;
                    block[iblk].y = j;
                    block[iblk].z = k;
                }
            }
        }
    }

#ifdef VERBOSE
    printf("At end of initial_placement.\n");
    dump_clbs();
#endif

    for (i = 0; i < num_types; ++i) {
        free(pos[i]);
    }

    free(pos);          /* Free the mapping list */
    free(index);
    free(count);
}


/* Frees the structures used to speed up evaluation of the nonlinear   *
 * congestion cost function.                                           */
static void free_fast_cost_update_structs(void)
{
    int i;
    for (i = 0; i <= num_grid_rows; i++) {
        free(chanx_place_cost_fac[i]);
    }

    free(chanx_place_cost_fac);

    for (i = 0; i <= num_grid_columns; i++) {
        free(chany_place_cost_fac[i]);
    }

    free(chany_place_cost_fac);
} 

/* Allocates and loads the chanx_place_cost_fac and chany_place_cost_fac *
 * arrays with the inverse of the average number of tracks per channel   *
 * between [subhigh] and [sublow].  This is only useful for the cost     *
 * function that takes the length of the net bounding box in each        *
 * dimension divided by the average number of tracks in that direction.  *
 * For other cost functions, you don't have to bother calling this       *
 * routine; when using the cost function described above, however, you   *
 * must always call this routine after you call init_chan and before     *
 * you do any placement cost determination.  The place_cost_exp factor   *
 * specifies to what power the width of the channel should be taken --   *
 * larger numbers make narrower channels more expensive.                 */
static void alloc_and_load_for_fast_cost_update(double place_cost_exp)
{    
    /* Access arrays below as chan?_place_cost_fac[subhigh][sublow].  Since   *
     * subhigh must be greater than or equal to sublow, we only need to       *
     * allocate storage for the lower half of a matrix.                       */
    int low, high, i;
    chanx_place_cost_fac = (double**)my_malloc((num_grid_rows + 1) * sizeof(double*));
    for (i = 0; i <= num_grid_rows; i++) {
        chanx_place_cost_fac[i] = (double*)my_malloc((i + 1) * sizeof(double));
    }

    chany_place_cost_fac = (double**)my_malloc((num_grid_columns + 1) * sizeof(double*));
    for (i = 0; i <= num_grid_columns; i++) {
        chany_place_cost_fac[i] = (double*)my_malloc((i + 1) * sizeof(double));
    }


    /* First compute the number of tracks between channel high and channel *
     * low, inclusive, in an efficient manner.                             */
    chanx_place_cost_fac[0][0] = chan_width_x[0];
    for (high = 1; high <= num_grid_rows; high++) {
        chanx_place_cost_fac[high][high] = chan_width_x[high];

        for (low = 0; low < high; low++) {
            chanx_place_cost_fac[high][low] =
                chanx_place_cost_fac[high - 1][low] +
                chan_width_x[high];
        }
    }

    /* Now compute the inverse of the average number of tracks per channel *
     * between high and low.  The cost function divides by the average     *
     * number of tracks per channel, so by storing the inverse I convert   *
     * this to a faster multiplication.  Take this final number to the     *
     * place_cost_exp power -- numbers other than one mean this is no      *
     * longer a simple "average number of tracks"; it is some power of     *
     * that, allowing greater penalization of narrow channels.             */
    for (high = 0; high <= num_grid_rows; high++)
        for (low = 0; low <= high; low++) {
            chanx_place_cost_fac[high][low] = (high - low + 1.) /
                                              chanx_place_cost_fac[high][low];
            chanx_place_cost_fac[high][low] =
                pow((double)chanx_place_cost_fac[high][low],
                    (double)place_cost_exp);
        }

    /* Now do the same thing for the y-directed channels.  First get the  *
     * number of tracks between channel high and channel low, inclusive.  */
    chany_place_cost_fac[0][0] = chan_width_y[0];
    for (high = 1; high <= num_grid_columns; high++) {
        chany_place_cost_fac[high][high] = chan_width_y[high];

        for (low = 0; low < high; low++) {
            chany_place_cost_fac[high][low] =
                chany_place_cost_fac[high - 1][low] +
                chan_width_y[high];
        }
    }

    /* Now compute the inverse of the average number of tracks per channel *
     * between high and low.  Take to specified power.                     */

    for (high = 0; high <= num_grid_columns; high++)
        for (low = 0; low <= high; low++) {
            chany_place_cost_fac[high][low] = (high - low + 1.) /
                                              chany_place_cost_fac[high][low];
            chany_place_cost_fac[high][low] =
                pow((double)chany_place_cost_fac[high][low],
                    (double)place_cost_exp);
        }
}

/* Checks that the placement has not confused our data structures. *
 * i.e. the clb and block structures agree about the locations of  *
 * every block, blocks are in legal spots, etc.  Also recomputes   *
 * the final placement cost from scratch and makes sure it is      *
 * within roundoff of what we think the cost is.                   */
static double check_place(double bb_cost,
                         double timing_cost,
                         int place_cost_type,
                         int num_regions,
                         place_algorithm_t place_algorithm,
                         double delay_cost)
{
    static int* bdone;
    int i, j, k, error = 0, block_num;
    int usage_check;

    double bb_cost_check = comp_bb_cost(CHECK, place_cost_type, num_regions);
    printf("bb_cost recomputed from scratch is %g.\n", bb_cost_check);
    if (fabs(bb_cost_check - bb_cost) > bb_cost * ERROR_TOL) {
        printf("Error:  bb_cost_check: %g and bb_cost: %g differ in check_place.\n",
               bb_cost_check, bb_cost);
        //error++;
    }

    /*
    //timing checking is not done here, since it will be recomputed - ie, it will be wrong for sure
    if(place_algorithm == NET_TIMING_DRIVEN_PLACE ||
    place_algorithm == PATH_TIMING_DRIVEN_PLACE)
    {
    comp_td_costs(&timing_cost_check, &delay_cost_check);
    printf("timing_cost recomputed from scratch is %g. \n",
    timing_cost_check);
    if(fabs(timing_cost_check - timing_cost) >
    timing_cost * ERROR_TOL)
    {
    printf("Error:  timing_cost_check: %g and timing_cost: "
    "%g differ in check_place.\n",
    timing_cost_check, timing_cost);
    //error++;
    }
    printf("delay_cost recomputed from scratch is %g. \n",
    delay_cost_check);
    if(fabs(delay_cost_check - delay_cost) > delay_cost * ERROR_TOL)
    {
    printf("Error:  delay_cost_check: %g and delay_cost: "
    "%g differ in check_place.\n",
    delay_cost_check, delay_cost);
    //error++;
    }
    }
     */
    bdone = (int*)my_malloc(num_blocks * sizeof(int));

    for (i = 0; i < num_blocks; i++) {
        bdone[i] = 0;
    }

    /* Step through grid array. Check it against block array. */
    for (i = 0; i <= (num_grid_columns + 1); i++)
        for (j = 0; j <= (num_grid_rows + 1); j++) {
            if (grid[i][j].usage > grid[i][j].type->capacity) {
                printf
                ("Error:  block at grid location (%d,%d) overused. "
                 "Usage is %d\n", i, j, grid[i][j].usage);
                error++;
            }

            usage_check = 0;

            for (k = 0; k < grid[i][j].type->capacity; k++) {
                block_num = grid[i][j].blocks[k];

                if (EMPTY == block_num) {
                    continue;
                }

                if (block[block_num].type != grid[i][j].type) {
                    printf
                    ("Error:  block %d type does not match grid location (%d,%d) type.\n",
                     block_num, i, j);
                    error++;
                }

                if ((block[block_num].x != i) || (block[block_num].y != j)) {
                    printf
                    ("Error:  block %d (%d,%d) location conflicts with grid(%d,%d)"
                     "data.\n", block_num, block[block_num].x, block[block_num].y, i, j);
                    error++;
                }

                ++usage_check;
                bdone[block_num]++;
            }

            if (usage_check != grid[i][j].usage) {
                printf
                ("Error:  Location (%d,%d) usage is %d, but has actual usage %d.\n",
                 i, j, grid[i][j].usage, usage_check);
            }
        }

    /* Check that every block exists in the grid and block arrays somewhere. */
    for (i = 0; i < num_blocks; i++)
        if (bdone[i] != 1) {
            printf
            ("Error:  block %d listed %d times in data structures.\n",
             i, bdone[i]);
            error++;
        }

    free(bdone);

    if (error == 0) {
        printf
        ("\nCompleted placement consistency check successfully.\n\n");
#ifdef PRINT_REL_POS_DISTR
        print_relative_pos_distr();
#endif
    } else {
        printf
        ("\nCompleted placement consistency check, %d Errors found.\n\n",
         error);
        printf("Aborting program.\n");
        //exit(1);
    }

    return bb_cost_check;
}

/*  Parallel Placement functions Startings....     */
/*partition nets evenly based on number of edges each net is connected to(that is sub-connections). */
static void balance_two_consecutive_threads_edge(int kthread_id)
{
    ++(start_finish_nets[kthread_id].counter_edge);
    unsigned long work_in_this_region = start_finish_nets[kthread_id].edges_in_this_partition;
    unsigned long work_in_next_region = start_finish_nets[kthread_id + 1].edges_in_this_partition;

    /*if next partition has more work (edges)
     *adjust the boundary according to % of difference */
    int net_shift = 0;
    double exceed_ratio = 0.0;
    int edge_partition_size = 0;
    if (work_in_this_region < work_in_next_region) {
        exceed_ratio = (work_in_next_region - work_in_this_region) / work_in_this_region;
        edge_partition_size = start_finish_nets[kthread_id + 1].edge_partition_size;
        net_shift = (int)(min(1, exceed_ratio) * edge_partition_size * 0.25);

        if (net_shift == 0 && (start_finish_nets[kthread_id + 1].edge_partition_size >= 2)) {
            net_shift = 1;
        }

        start_finish_nets[kthread_id].finish_edge += net_shift;
        start_finish_nets[kthread_id + 1].start_edge += net_shift;
    } else if (work_in_this_region > work_in_next_region) {
        /* if this partition has more edges
         * adjust the boundary according to % of difference */
        exceed_ratio = (work_in_this_region - work_in_next_region) / work_in_next_region;
        edge_partition_size = start_finish_nets[kthread_id].edge_partition_size;
        net_shift = (int)(min(1, exceed_ratio) * edge_partition_size * 0.25);

        if (net_shift == 0 && start_finish_nets[kthread_id].edge_partition_size >= 2) {
            net_shift = 1;
        }

        start_finish_nets[kthread_id].finish_edge -= net_shift;
        start_finish_nets[kthread_id + 1].start_edge -= net_shift;
    } else {
        /* No operations */
    }

    assert(start_finish_nets[kthread_id].finish_edge
             > start_finish_nets[kthread_id].start_edge);
} /* end of void balance_two_consecutive_threads_edge(int kthread_id)  */

/*partition nets evenly based on number of sinks each net is connected to*/
static void balance_two_consecutive_threads_sinks(int kthread_id)
{
    ++(start_finish_nets[kthread_id].counter_sink);
    unsigned long work_in_this_region = start_finish_nets[kthread_id].sinks_in_this_partition;
    unsigned long work_in_next_region = start_finish_nets[kthread_id + 1].sinks_in_this_partition;

    /*if next partition has more edges
     *adjust the boundary according to % of difference */
    int net_shift = 0;
    double exceed_ratio = 0.0;
    int sink_partition_size = 0;
    if (work_in_this_region < work_in_next_region) {
        exceed_ratio = (work_in_next_region - work_in_this_region) / work_in_this_region;
        sink_partition_size = start_finish_nets[kthread_id + 1].sink_partition_size;
        net_shift = (int)(min(1, exceed_ratio) * sink_partition_size * 0.25);

        start_finish_nets[kthread_id].finish_sinks += net_shift;
        start_finish_nets[kthread_id + 1].start_sinks += net_shift;
    } else if (work_in_this_region > work_in_next_region) {
        /*if this partition has more edges
         *adjust the boundary according to % of difference */
        exceed_ratio = (work_in_this_region - work_in_next_region) / work_in_next_region;
        sink_partition_size = start_finish_nets[kthread_id].sink_partition_size;
        net_shift = (int)(min(1, exceed_ratio) * sink_partition_size * 0.25);

        start_finish_nets[kthread_id].finish_sinks -= net_shift;
        start_finish_nets[kthread_id + 1].start_sinks -= net_shift;
    } else {
        /* No operations */
    }

    assert (start_finish_nets[kthread_id].finish_sinks
              > start_finish_nets[kthread_id].start_sinks);
}  /* end of void balance_two_consecutive_threads_sinks(int kthread_id) */


/* Copy memory from bb_num_on_edges(or bb_coords) to *
 * local_bb_edge(local_bb_coord) */
static void tp_init_local(const int*  region_x_boundary,
                          thread_local_data_for_swap_t*  swap_data_ptr)
{
    memcpy(swap_data_ptr->m_local_bb_edge,
           bb_num_on_edges,
           num_nets * sizeof(bbox_t));
    memcpy(swap_data_ptr->m_local_bb_coord,
           bb_coords,
           num_nets * sizeof(bbox_t));

    int x, y, z;
    for (x = 0; x < num_nets; ++x) {
        swap_data_ptr->m_local_temp_net_cost[x] = temp_net_cost[x];
        swap_data_ptr->m_local_net_cost[x] = net_cost[x];
    }

    /* extend sub-region[start_boundary-2..end_boundary_x+2] */
    grid_tile_t** local_grid = swap_data_ptr->m_local_grid;
    for (x = region_x_boundary[0] - 2;
            x < region_x_boundary[2] + 2 && x < (num_grid_columns + 2); ++x) {
        /*  takes care of proc #1, where start_end_boundary = 0 */
        if (x < 0) {
            x = 0;
        }
        /* Why not set y as Extend-SubRegion[region_y_boundary-2..region_y_boundary+2] */
        for (y = 0; y <= (num_grid_rows + 1); ++y) {
            local_grid[x][y].type = grid[x][y].type;
            local_grid[x][y].usage = grid[x][y].usage;
            local_grid[x][y].offset = grid[x][y].offset;

            local_grid[x][y].blocks =
                (int*)my_malloc(sizeof(int) * grid[x][y].type->capacity);

            for (z = 0; z < grid[x][y].type->capacity; ++z) {
                local_grid[x][y].blocks[z] = grid[x][y].blocks[z];
            }
        }
    }

    local_block_t* local_block = swap_data_ptr->m_local_block;
    for (x = 0; x < num_blocks; ++x) {
        local_block[x].x = block[x].x;
        local_block[x].y = block[x].y;
        local_block[x].z = block[x].z;
    }
}  /* end of void tp_init_local(bbox_t* local_bb_edge,...) */


static void tp_init_localvert_grid()
{
    /* grid[col][row] was column-based, but localvert_grid[col][row] was
     * row-based, it vertical to grid. */
    localvert_grid = (grid_tile_t**)alloc_matrix(0, (num_grid_columns + 1),
                                                 0, (num_grid_rows + 1),
                                                 sizeof(grid_tile_t));
    int col, row, cap;
    for (col = 0; col < (num_grid_columns + 2); ++col) {
        for (row = 0; row < (num_grid_rows + 2); ++row) {
            localvert_grid[row][col].type = grid[col][row].type;
            localvert_grid[row][col].usage = grid[col][row].usage;
            localvert_grid[row][col].offset = grid[col][row].offset;

            localvert_grid[row][col].blocks =
                (int*)my_malloc(sizeof(int) * grid[col][row].type->capacity);

            for (cap = 0; cap < grid[col][row].type->capacity; ++cap) {
                localvert_grid[row][col].blocks[cap] = grid[col][row].blocks[cap];
            }
        }
    }
} /* end of void tp_init_localvert_grid() */


/* FIXME, important for update timing parallel */ 
static void tp_timing_update_full(int kthread_id,
                                  double*** net_delay,
                                  pthread_data_t* input_args,
                                  double*   max_delay)
{
    /*dynamic workload setup*/
    const int counter_sink = start_finish_nets[kthread_id].counter_sink;
    if (kthread_id != NUM_OF_THREADS - 1
          && counter_sink < (num_nets / PARITION_UPDATE)) {
        const int partition_size = (start_finish_nets[kthread_id].finish_sinks
                                       - start_finish_nets[kthread_id].start_sinks);
            start_finish_nets[kthread_id].sink_partition_size = partition_size;
        }

    /* loads the net_delay, and return the total number of sinks visited for dynamic workload */
    const int cur_thread_start_sinks = start_finish_nets[kthread_id].start_sinks;
    const int cur_thread_finish_sinks = start_finish_nets[kthread_id].finish_sinks;
    start_finish_nets[kthread_id].sinks_in_this_partition =
                   load_timing_graph_net_delays_parallel(*net_delay,
                                                         cur_thread_start_sinks,
                                                         cur_thread_finish_sinks);

    barrier_polling(kthread_id);
    /* update dynamic workload distribution*/
    if (kthread_id != NUM_OF_THREADS - 1 && counter_sink < num_nets / PARITION_UPDATE) {
        balance_two_consecutive_threads_sinks(kthread_id);
    }
    /*  reset global variables  */
    if (kthread_id == 0) {
        *(input_args->timing_cost) = 0;
        *(input_args->delay_cost) = 0;
        *(input_args->max_delay) = 0.0;
    }

    /*  Load-net-slack starting....   */
    /* part 1 *
     * Reset all arrival times to -ve infinity. Can't just set to zero or the *
     * constant propagation(constant generators work at -ve infinity) won't  *
     * work.                                                                 */
    int tnodes_assign_to_thread = ceil((double)num_tnodes / NUM_OF_THREADS);
    int start_node = kthread_id * tnodes_assign_to_thread;
    int finish_node = min((kthread_id + 1) * tnodes_assign_to_thread, num_tnodes);
    int i = 0;
    for (i = start_node; i < finish_node; ++i) {
        tnode[i].arr_time = T_CONSTANT_GENERATOR; /* -1000 */
    }

    barrier_polling(kthread_id);
    /* Part 2 *
     * reset arrivial time for all nodes at level 0. */
    int tnodes_at_level0 = tnodes_at_level[0].nelem;
    tnodes_assign_to_thread = ceil((double)tnodes_at_level0 / NUM_OF_THREADS);
    start_node = kthread_id * tnodes_assign_to_thread;
    finish_node = min((kthread_id + 1) * tnodes_assign_to_thread,
                       tnodes_at_level0);
    for (i = start_node; i < finish_node; ++i) {
        int node_index = tnodes_at_level[0].list[i];
        tnode[node_index].arr_time = 0.0;
    }
    barrier_polling(kthread_id);

    /* Part 3, compute all tnodes arrival_time parallely! And find out the *
     * Critical_delay. The functions are parallelized on a per-level basis *
     * each processor will receive a min. of 100 nodes to work with.       *
     * If there are not enough nodes to distrubute to all processors,
     * processor with larger kthread_id will remain idle.       */
    *max_delay = 0.0;
    int ilevel = 0;
    for (ilevel = 1; ilevel < num_tnode_levels; ++ilevel) {
        barrier_polling(kthread_id);
        int num_at_level = tnodes_at_level[ilevel].nelem;
        int num_of_thread_used = ceil((double)num_at_level / 100);
        tnodes_assign_to_thread = ceil((double)num_at_level /
                                       min(NUM_OF_THREADS, num_of_thread_used));

        /* for a big enough partition */
        if (kthread_id < num_of_thread_used) {
            start_node = kthread_id * tnodes_assign_to_thread;
            finish_node = min((kthread_id + 1) * tnodes_assign_to_thread,
                               num_at_level);
            /* the following line does the actual work.
               rest of the stuff above is for the workload distrubtion */
            *max_delay = max(*max_delay,
                             calc_tnodes_arr_time_parallel(start_node,
                                                           finish_node,
                                                           ilevel));
        }
    }  /* end of for(ilevel = 1; ilevel < num_tnode_levels; ++ilevel) */

    /* the MAX_DELAY value is written back */
    pthread_mutex_lock(&global_data_access.mutex);
    *(input_args->max_delay) = max(*(input_args->max_delay), *max_delay);
    pthread_mutex_unlock(&global_data_access.mutex);

    /* Part 4, compute all tnodes required_time parallely!  *
     * same concept as part 3, but for a different function */
    barrier_polling(kthread_id);
    *max_delay = *(input_args->max_delay);
    for (ilevel = num_tnode_levels - 1; ilevel >= 0; --ilevel) {
        int num_at_level = tnodes_at_level[ilevel].nelem;
        int num_of_thread_used = ceil((double)num_at_level / 100);
        tnodes_assign_to_thread = ceil((double)num_at_level /
                                       min(NUM_OF_THREADS, num_of_thread_used));

        /* for a big enough partition */
        if (kthread_id < num_of_thread_used) {
            start_node = kthread_id * tnodes_assign_to_thread;
            finish_node = min((kthread_id + 1) * tnodes_assign_to_thread,
                               num_at_level);
            calc_tnodes_req_time_parallel(*max_delay,
                                          start_node,
                                          finish_node,
                                          ilevel);
        }
        barrier_polling(kthread_id);
    }  /* end of for(ilevel = num_tnode_levels - 1; ilevel >= 0; --ilevel)*/
}  /* end of void tp_timing_update_full(int kthread_id,...)  */

/* FIXME, If I want to compute net_slack, first I must calcuatel all edge's delay value,
 * then calcuate all vertexes' arr_time and req_time. Last compute slack       */
static void tp_compute_net_slack(int kthread_id,
                                 double***  net_slack,
                                 double*  timing_cost,
                                 double*  delay_cost,
                                 double  crit_exponent,
                                 double  max_delay,
                                 pthread_data_t* input_args)
{
    /* part 5
     * compute net slacks */
    const int thread_start_edge = start_finish_nets[kthread_id].start_edge;
    const int thread_finish_edge = start_finish_nets[kthread_id].finish_edge;
    const int partition_size = thread_finish_edge - thread_start_edge;

    const int thread_counter_edge = start_finish_nets[kthread_id].counter_edge;
    if (thread_counter_edge < (num_nets / PARITION_UPDATE)) {
        start_finish_nets[kthread_id].edge_partition_size = partition_size;
    }

    /* Compute [thread_start_edge, thread_finish_edge] nets' slack OK! */
    start_finish_nets[kthread_id].edges_in_this_partition =
                                 compute_net_slacks_parallel(*net_slack,
                                                             thread_start_edge,
                                                             thread_finish_edge);
    barrier_polling(kthread_id);
    /* dynamic workload distribution *
     * equalize partition */
    if (kthread_id != NUM_OF_THREADS - 1
          && thread_counter_edge < (num_nets / PARITION_UPDATE)) {
        /* update and balance the thread and (thread + 1)'s start_edge and
         * finish_edge */ 
        balance_two_consecutive_threads_edge(kthread_id);
    }

    const int thread_start_sinks = start_finish_nets[kthread_id].start_sinks;
    const int thread_finish_sinks = start_finish_nets[kthread_id].finish_sinks;
    const int thread_counter_sink = start_finish_nets[kthread_id].counter_sink;
    if (kthread_id != NUM_OF_THREADS - 1
          && thread_counter_sink < (num_nets / PARITION_UPDATE)) {
        start_finish_nets[kthread_id].sink_partition_size =
                              thread_finish_sinks - thread_start_sinks;
    }

    /* Part 6, compute timing_driven cost Parallel */
    start_finish_nets[kthread_id].sinks_in_this_partition =
                            compute_td_costs_parallel_with_update_crit(timing_cost,
                                                                       delay_cost,
                                                                       *net_slack,
                                                                       max_delay,
                                                                       crit_exponent,
                                                                       thread_start_sinks,
                                                                       thread_finish_sinks);

    /* write back the partial timing and Tdel cost to the global variable */
    /* Why did author using these 2 following variables? */
    partial_results[kthread_id] = *timing_cost;
    partial_results2[kthread_id] = *delay_cost;

    barrier_polling(kthread_id);
    /* master thread sums up the partial values */
    if (kthread_id == 0) {
        *timing_cost = 0.0;
        *delay_cost = 0.0;
        int thread_idx = 0;
        for (thread_idx = 0; thread_idx < NUM_OF_THREADS; ++thread_idx) {
            *timing_cost += partial_results[thread_idx];
            *delay_cost += partial_results2[thread_idx];
        }

        *(input_args->timing_cost) = *timing_cost;
        *(input_args->delay_cost) = *delay_cost;
    }  /* end of if(kthread_id == 0) */

    barrier_polling(kthread_id);
    /* Update partition size  */
    if (kthread_id != NUM_OF_THREADS - 1
          && thread_counter_sink < num_nets / PARITION_UPDATE) {
        balance_two_consecutive_threads_sinks(kthread_id);
    }
}  /* end of void tp_compute_net_slack(int kthread_id,..) */

static void tp_timing_calc(int kthread_id,
                           double* timing_cost,
                           double* delay_cost,
                           pthread_data_t* input_args)
{
    const int thread_counter_sink = start_finish_nets[kthread_id].counter_sink;
    const int thread_finish_sinks = start_finish_nets[kthread_id].finish_sinks;
    const int thread_start_sinks = start_finish_nets[kthread_id].start_sinks;
    if (kthread_id != NUM_OF_THREADS - 1 && thread_counter_sink < num_nets / PARITION_UPDATE) {
        const int partition_size = thread_finish_sinks - thread_start_sinks;
        start_finish_nets[kthread_id].sink_partition_size = partition_size;
    }

    start_finish_nets[kthread_id].sinks_in_this_partition =
                    compute_td_costs_parallel_without_update_crit(timing_cost,
                                                                  delay_cost,
                                                                  thread_start_sinks,
                                                                  thread_finish_sinks);
    pthread_mutex_lock(&global_data_access.mutex);
    *(input_args->timing_cost) += *timing_cost;
    *(input_args->delay_cost) += *delay_cost;
    pthread_mutex_unlock(&global_data_access.mutex);
    barrier_polling(kthread_id);

    if (kthread_id != NUM_OF_THREADS - 1 && thread_counter_sink < num_nets / PARITION_UPDATE) {
        balance_two_consecutive_threads_sinks(kthread_id);
    }
}  /* end of static void tp_timing_calc(int kthread_id,..) */


static void tp_iter_data_update(const int kiter,
                                pthread_data_t* input_args,
                                thread_local_common_paras_t*  common_paras_ptr,
                                thread_local_data_for_swap_t* swap_data_ptr)
{
    /* Regions in first row of the region grid */
    const int  kthread_id = common_paras_ptr->local_thread_id;
    const int* kregion_x_boundary = common_paras_ptr->local_region_x_boundary;
    const int* kregion_y_boundary = common_paras_ptr->local_region_y_boundary;
    const place_algorithm_t kplace_algorithm =
        common_paras_ptr->local_placer_opts.place_algorithm;
    if (kthread_id / sqrt(NUM_OF_THREADS) == 0) {
    /* Attention, current NUM_OF_THREADS == 4, so when kthread_id = 0 or 1, *
     * the Regions was in first row. Later, I will use Intel Xeon 4-core   *
     * 8-threads CPU, I'd like to use 8 threads parallel. When I set 4x2,  *
     * that is 4 cols in each row, and totally 2 rows. You can see x_partion*
     * and y_partition in place.c. Sqrt(8) = 3.464, when kthread_id = 0,1,2,3,
     * the regions was in first row. But when I set 2x4 for 8 threads(4 rows
     * and 2 cols), this may be error! */
        if (kplace_algorithm == NET_TIMING_DRIVEN_PLACE ||
                kplace_algorithm == PATH_TIMING_DRIVEN_PLACE) {
            *(input_args->cost) = common_paras_ptr->local_total_cost = 1.0;
        }
        /* all average cost initial as 0.0 */
        if (kiter == 0 && kthread_id == 0) {
            *(input_args->av_cost) = 0.0;
            *(input_args->av_bb_cost) = 0.0;
            *(input_args->av_timing_cost) = 0.0;
            *(input_args->av_delay_cost) = 0.0;
            *(input_args->sum_of_squares) = 0.0;
            *(input_args->success_sum) = 0.0;
            *(input_args->move_lim) = 0;
        }
    } else {
    /* update local data from global for top rows of each private region */
        update_from_global_to_local_grid_only(swap_data_ptr->m_local_grid,
                                              max(0,  kregion_x_boundary[0] - 2),
                                              kregion_x_boundary[0] + 2,
                                              max(kregion_y_boundary[0] - 2 , 0),
                                              min(kregion_y_boundary[2] + 2,
                                                  num_grid_rows + 2));
    }

    /* global grid to local block grid update */
    int x = -1;
    local_block_t* local_block = swap_data_ptr->m_local_block;
    for (x = 0; x < num_blocks; ++x) {
        local_block[x].x = block[x].x;
        local_block[x].y = block[x].y;
        local_block[x].z = block[x].z;
    }

    /* timing update. */
    if (kplace_algorithm == NET_TIMING_DRIVEN_PLACE ||
            kplace_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        common_paras_ptr->local_total_cost = 1.0;

        common_paras_ptr->local_timing_cost = *(input_args->timing_cost);
        common_paras_ptr->local_inverse_prev_timing_cost =
            1 / common_paras_ptr->local_timing_cost;
        common_paras_ptr->local_delay_cost = *(input_args->delay_cost);
    }

    common_paras_ptr->local_bb_cost = *(input_args->bb_cost);
    common_paras_ptr->local_inverse_prev_bb_cost =
        1 / common_paras_ptr->local_bb_cost;

    memcpy(swap_data_ptr->m_local_bb_edge,
           bb_num_on_edges,
           num_nets * sizeof(bbox_t));
    memcpy(swap_data_ptr->m_local_bb_coord,
           bb_coords,
           num_nets * sizeof(bbox_t));
    memcpy(swap_data_ptr->m_local_temp_net_cost,
           temp_net_cost,
           num_nets * sizeof(double));
    memcpy(swap_data_ptr->m_local_net_cost,
           net_cost,
           num_nets * sizeof(double));
} /* end of void tp_iter_data_update(int kthread_id... ) */

/* Update local sub-region(in a Region) horizontal and vertical seperately */
static void tp_local_data_update(const int* region_x_boundary,
                                 const int* region_y_boundary,
                                 grid_tile_t**  local_grid,
                                 local_block_t* local_block,
                                 const int krow,
                                 const int kcol)
{
    /* update local grid information, first was horizontal direction, 
     * then was vertical direction. */
    if (krow == 0 && kcol == 0) {
        /* update top of the strip(horizontal direction) */
        update_from_global_to_local_hori(local_block,
                                         local_grid,
                                         max(2, region_x_boundary[0] - 2),
                                         region_x_boundary[0] + 2,
                                         max(0, region_y_boundary[0] + 2),
                                         region_y_boundary[1] + 2);
        /* update left side of the strip(vertical direction) */
        update_from_global_to_local_vert(local_block,
                                         local_grid,
                                         max(0, region_x_boundary[0] - 2),
                                         region_x_boundary[1] + 2,
                                         max(2, region_y_boundary[0] - 2),
                                         region_y_boundary[0] + 2);
    } else if (krow == 0 && kcol == 1) {
        /* update top of the strip(horizontal direction) */
        update_from_global_to_local_hori(local_block,
                                         local_grid,
                                         max(2, region_x_boundary[0] - 2),
                                         region_x_boundary[0] + 2,
                                         region_y_boundary[1] + 2,
                                         min(num_grid_rows + 2,
                                             region_y_boundary[2] + 2));
        /* update right side of the strip(vertical direction) */
        update_from_global_to_local_vert(local_block,
                                         local_grid,
                                         max(0, region_x_boundary[0] - 2),
                                         region_x_boundary[1] + 2,
                                         region_y_boundary[2] - 2,
                                         min(num_grid_rows - 2,
                                             region_y_boundary[2] + 2));
    } else if (krow == 1 && kcol == 0) {
        /* update bottom of the strip(horizontal directions) */
        update_from_global_to_local_hori(local_block,
                                         local_grid,
                                         region_x_boundary[2] - 2,
                                         min(num_grid_columns - 2,
                                             region_x_boundary[2] + 2),
                                         max(0, region_y_boundary[0] + 2),
                                         region_y_boundary[1] + 2);
        /* update left side strip(vertical directions) */
        update_from_global_to_local_vert(local_block,
                                         local_grid,
                                         region_x_boundary[1] - 2, 
                                         min(num_grid_columns + 2,
                                             region_x_boundary[2] + 2),
                                         max(0, region_y_boundary[0] - 2),
                                         region_y_boundary[0] + 2);
    } else if (krow == 1 && kcol == 1) {
        /* update bottom of the strip(horizontal direction) */
        update_from_global_to_local_hori(local_block,
                                         local_grid,
                                         region_x_boundary[2] - 2,
                                         min(num_grid_columns - 2,
                                             region_x_boundary[2] + 2),
                                         region_y_boundary[1] + 2,
                                         min(num_grid_columns + 2,
                                             region_y_boundary[2] + 2));
        /* update right side strip(vertical direction) */
        update_from_global_to_local_vert(local_block,
                                         local_grid,
                                         region_x_boundary[1] - 2,
                                         min(num_grid_columns + 2,
                                             region_x_boundary[2] + 2),
                                         region_y_boundary[2] - 2,
                                         min(num_grid_rows - 2,
                                             region_y_boundary[2] + 2));
    } else {
        printf("incorrect row,col combination in inner loop: (%d, %d)", krow, kcol);
    }
}  /* end of void tp_local_data_update(int* region_x_boundary) */


/* Update the temperature according to the annealing schedule selected. */
static void update_t_parallel(double* t,
                              int* inner_iter_num,
                              double range_limit,
                              double success_rat,
                              annealing_sched_t annealing_sched)
{
    /*  double fac; */
    if (annealing_sched.type == USER_SCHED) {
        *t = annealing_sched.alpha_t * (*t);
    }
    /* Old standard deviation based stuff is below.  This bogs down horribly
     * for big circuits (alu4 and especially bigkey_mod). */
    /* #define LAMBDA .7  */
    /* ------------------------------------ */
#if 0
    else if (std_dev == 0.) {
        *t = 0.;
    } else {
        fac = exp(-LAMBDA * (*t) / std_dev);
        fac = max(0.5, fac);
        *t = (*t) * fac;
    }

#endif
    /* ------------------------------------- */
    else {
        /* AUTO_SCHED */
        *inner_iter_num = annealing_sched.inner_num;
        if (success_rat > 0.96) {
            *t = (*t) * 0.5;
        } else if (success_rat > 0.8) {
            *t = (*t) * 0.9;
        } else if (success_rat > 0.15 || range_limit > 1.) {
            *t = (*t) * 0.9;
            *inner_iter_num = annealing_sched.inner_num / 4;
        } else {
            *t = (*t) * 0.6;
            *inner_iter_num = annealing_sched.inner_num / 20;
        }
    }
}

/* Puts a list of all the nets connected to from_block and to_block into  *
 * nets_to_update.  Returns the number of affected nets.  Net_block_moved *
 * is either FROM, TO or FROM_AND_TO -- the block connected to this net   *
 * that has moved.                                                        */
static int find_affected_nets_parallel(int* nets_to_update,
                                       int* net_block_moved,
                                       int from_block,
                                       int to_block,
                                       int num_of_pins,
                                       double* local_temp_net_cost)
{
    int affected_index = 0;
    int k, inet, count;
    for (k = 0; k < num_of_pins; ++k) {
        inet = block[from_block].nets[k];
        if (OPEN == inet || TRUE == net[inet].is_global) {
            continue;
        }
        /* This is here in case the same block connects to a net twice. */
        if (local_temp_net_cost[inet] > 0.0) {
            continue;
        }

        nets_to_update[affected_index] = inet;
        net_block_moved[affected_index] = FROM;
        ++affected_index;
        local_temp_net_cost[inet] = 1.0; /* Flag to say we've marked this net. */
    }  /* end of for (k = 0; k < num_of_pins; ++k) */

    if (to_block != EMPTY) {
        for (k = 0; k < num_of_pins; ++k) {
            inet = block[to_block].nets[k];
            if (OPEN == inet || TRUE == net[inet].is_global) {
                continue;
            }

            if (local_temp_net_cost[inet] > 0.0) {
                /* Net already marked. */
                for (count = 0; count < affected_index; ++count) {
                    if (nets_to_update[count] == inet) {
                        if (net_block_moved[count] == FROM) {
                            net_block_moved[count] = FROM_AND_TO;
                        }
                        break;
                    }
                }

#ifdef DEBUG
                if (count > affected_index) {
                    printf("Error in find_affected_nets -- count = %d,"
                     " affected index = %d.\n", count,
                     affected_index);
                    exit(1);
                }
#endif
            } else {
                /* Net not marked yet. */
                nets_to_update[affected_index] = inet;
                net_block_moved[affected_index] = TO;
                ++affected_index;
                local_temp_net_cost[inet] = 1.; /* Flag means we've  marked net. */
            }
        }
    }  /* end of if(to_block != NULL) */

    return affected_index;
}  /* end of static int find_affected_nets_parallel(int* nets_to_update,) */

static int assess_swap_parallel(double delta_c,
                                double t,
                                int local_seed)
{
    /* Returns: 1->move accepted, 0->rejected. */
    int accept = -1;
    if (delta_c <= 0) {
#ifdef SPEC  /* Reduce variation in final solution due to round off */
        fnum = my_frand_parallel(local_seed);
#endif
        accept = 1;
        return (accept);
    }

    if (t == 0.) {
        return (0);
    }

    double fnum = my_frand_parallel(local_seed);
    double prob_fac = exp(-delta_c / t);
    if (prob_fac > fnum) {
        accept = 1;
    } else {
        accept = 0;
    }

    return (accept);
}  /* end of int assess_swap_parallel(double delta_c,) */

static boolean find_to_block_parallel(int x_from,
                                      int y_from,
                                      block_type_ptr type,
                                      int* x_to,
                                      int* y_to,
                                      int kthread_id,
                                      int xmin, int xmax,
                                      int ymin, int ymax)
{
    /* Returns the location to which I want to swap within the range (xmin,ymin) & (xmax,ymax)*/
    int test;
    do {
        /* Until (x_to, y_to) different from (x_from, y_from) */
        if (type == IO_TYPE) {
            /* io_block to be moved. */
            /*left bottom corner*/
            if (xmin == 0 && ymin == 0) {
                /*the target is either located on x = 0 or y = 0*/
                if (my_irand_parallel(1, kthread_id)) {
                    /*target is located on x = 0*/
                    *x_to = 0;
                    test = my_irand_parallel(ymax - ymin, kthread_id);
                    *y_to = max(1, test);
                } else {
                    /*target is located on y = 0*/
                    test = my_irand_parallel(xmax - xmin, kthread_id);
                    *x_to = max(1, test);
                    *y_to = 0;
                }
            } else if (xmax == num_grid_columns + 1 && ymin == 0) {
                /* For Bottom-Right corner*/
                /*the target is either located on x = 0 or y = 0*/
                if (my_irand_parallel(1, kthread_id)) {
                    /*target is located on x = num_grid_columns + 1*/
                    *x_to = num_grid_columns + 1;
                    test = my_irand_parallel(ymax - ymin, kthread_id);
                    *y_to = max(1, test);
                } else {
                    /*target is located on y = 0*/
                    test = xmin + my_irand_parallel(xmax - xmin, kthread_id);
                    *x_to = min(num_grid_columns, test);
                    *y_to = 0;
                }
            } else if (xmin == 0 && ymax == num_grid_rows + 1) {
                /* For Top-Left corner*/
                /*the target is either located on x = 0 or y = 0*/
                if (my_irand_parallel(1, kthread_id)) {
                    /*target is located on x = 0*/
                    *x_to = 0;
                    test = ymin + my_irand_parallel(ymax - ymin, kthread_id);
                    *y_to = min(num_grid_rows, test);
                } else {
                    /*target is located on y = num_grid_rows+1*/
                    test = my_irand_parallel(xmax - xmin, kthread_id);
                    *x_to = max(1, test);
                    *y_to = num_grid_rows + 1;
                }
            } else if (xmax == num_grid_columns + 1 && ymax == num_grid_rows + 1) {
                /* For Top-Right corner*/
                /*the target is either located on x = 0 or y = 0*/
                if (my_irand_parallel(1, kthread_id)) {
                    /*target is located on x = num_grid_columns + 1*/
                    *x_to = num_grid_columns + 1;
                    test = ymin + my_irand_parallel(ymax - ymin, kthread_id);
                    *y_to = min(num_grid_rows, test);
                } else {
                    /*target is located on y = 0*/
                    test = xmin + my_irand_parallel(xmax - xmin, kthread_id);
                    *x_to = min(num_grid_columns, test);
                    *y_to = num_grid_rows + 1;
                }
            } else if (xmin == 0 || xmax == num_grid_columns + 1) {
                *x_to = x_from;
                *y_to = ymin + my_irand_parallel(ymax - ymin, kthread_id);
            } else if (ymin == 0 || ymax == num_grid_rows + 1) {
                *x_to = xmin + my_irand_parallel(xmax - xmin, kthread_id);
                *y_to = y_from;
            } else {
                printf("came to a wrong loop\n");
                exit(-1);
            }

            assert(type == grid[*x_to][*y_to].type);
        } else {
            /* For other types except IO_TYPES */
            /* generate a {x_offset, y_offset} pairs that between (xmin, xmax) *
             * (ymin, ymax). */
            int x_offset = my_irand_parallel(xmax - xmin, kthread_id);
            int y_offset = my_irand_parallel(ymax - ymin, kthread_id);

            *x_to = x_offset + xmin;
            *y_to = y_offset + ymin;

            /* make sure the destination is not an IO */
            *x_to = max(1, *x_to);
            *x_to = min(num_grid_columns, *x_to);
            *y_to = max(1, *y_to);
            *y_to = min(num_grid_rows, *y_to);

            assert(*x_to >= 1 && *x_to <= num_grid_columns);
            assert(*y_to >= 1 && *y_to <= num_grid_rows);
        }
    } while ((x_from == *x_to) && (y_from == *y_to));

    assert(type == grid[*x_to][*y_to].type);
    return TRUE;
}  /* end of static boolean find_to_block_parallel(int x_from,) */

static double comp_td_point_to_point_delay_parallel(int inet,
                                                    int ipin,
                                                    local_block_t* local_block)
{
    /*returns the Tdel of one point to point connection */
    int source_block = net[inet].node_block[0];
    int sink_block = net[inet].node_block[ipin];

    block_type_ptr source_type = block[source_block].type;
    block_type_ptr sink_type = block[sink_block].type;
    assert(source_type != NULL && sink_type != NULL);

    int delta_x = abs(local_block[sink_block].x - local_block[source_block].x);
    int delta_y = abs(local_block[sink_block].y - local_block[source_block].y);

    /* TODO low priority: Could be merged into one look-up table */
    /* Note: This heuristic is terrible on Quality of Results.
     * A much better heuristic is to create a more comprehensive lookup table but
     * it's too late in the release cycle to do this.  Pushing until the next release */
    double delay_source_to_sink = 0.0;
    if (source_type == IO_TYPE) {
        if (sink_type == IO_TYPE) {
            delay_source_to_sink = delta_io_to_io[delta_x][delta_y];
        } else {
            delay_source_to_sink = delta_io_to_fb[delta_x][delta_y];
        }
    } else {
        if (sink_type == IO_TYPE) {
            delay_source_to_sink = delta_fb_to_io[delta_x][delta_y];
        } else {
            delay_source_to_sink = delta_fb_to_fb[delta_x][delta_y];
        }
    }

    if (delay_source_to_sink < 0) {
        printf
        ("Error in comp_td_point_to_point_delay in place.c, bad delay_source_to_sink value\n");
        exit(1);
    }

    if (delay_source_to_sink < 0.) {
        printf
        ("Error in comp_td_point_to_point_delay in place.c, Tdel is less than 0\n");
        exit(1);
    }

    return (delay_source_to_sink);
}


/*a net that is being driven by a moved block must have all of its  */
/*sink timing costs recomputed. A net that is driving a moved block */
/*must only have the timing cost on the connection driving the input_args */
/*pin computed */
static void comp_delta_td_cost_parallel(int from_block,
                                        int to_block,
                                        int num_of_pins,
                                        double* delta_timing,
                                        double* delta_delay,
                                        local_block_t* local_block,
                                        double** local_temp_point_to_point_timing_cost,
                                        double** local_temp_point_to_point_delay_cost)
{
    double delta_timing_cost = 0.0;
    double delta_delay_cost = 0.0;
    double temp_delay = 0.0;
    /* int inet, pin_index, net_pin, ipin; */
    int pin_index = -1;
    int ipin = 0;
    for (pin_index = 0; pin_index < num_of_pins; ++pin_index) {
        const int inet = block[from_block].nets[pin_index];
        if (OPEN == inet || TRUE == net[inet].is_global) {
            continue;
        }

        const int net_pin = net_pin_index[from_block][pin_index];
        if (net_pin != 0) {
        /* This net_pin is not a driver_pin, so this net drived a moved block. *
         * if this net is being driven by a block that has moved, we do not  *
         * need to compute the change in the timing cost (here) since it will *
         * be computed in the fanout of the net on the driving block, also  *
         * computing it here would double count the change, and mess up the  *
         * delta_timing_cost value */
            if (net[inet].node_block[0] != to_block
                    && net[inet].node_block[0] != from_block) {
            /* In this situation, the drivering node of this net neither *
             * from_block nor to_block */
                temp_delay = comp_td_point_to_point_delay_parallel(inet,
                                                                   net_pin,
                                                                   local_block);

                local_temp_point_to_point_delay_cost[inet][net_pin] = temp_delay;
                local_temp_point_to_point_timing_cost[inet][net_pin] =
                                  timing_place_crit[inet][net_pin] * temp_delay;

                delta_delay_cost +=
                            (local_temp_point_to_point_delay_cost[inet][net_pin]
                               - point_to_point_delay_cost[inet][net_pin]);
                delta_timing_cost +=
                    (local_temp_point_to_point_timing_cost[inet][net_pin]
                       - point_to_point_timing_cost[inet][net_pin]);
            }
        } else { /* net_pin == 0 */
        /* this net is being driven by a moved block, recompute *
         * all point to point connections on this net. */
            for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                temp_delay = comp_td_point_to_point_delay_parallel(inet,
                                                                   ipin,
                                                                   local_block);
                local_temp_point_to_point_delay_cost[inet][ipin] = temp_delay;
                local_temp_point_to_point_timing_cost[inet][ipin] =
                                    timing_place_crit[inet][ipin] * temp_delay;

                delta_delay_cost +=
                    (local_temp_point_to_point_delay_cost[inet][ipin]
                       - point_to_point_delay_cost[inet][ipin]);
                delta_timing_cost +=
                    (local_temp_point_to_point_timing_cost[inet][ipin]
                       - point_to_point_timing_cost[inet][ipin]);
            }
        }
    }  /* end of for (pin_index = 0; pin_index < num_of_pins; ++pin_index) */

    /* end of from_block, then consider the to_block. It was similar with former program. */
    if (to_block != EMPTY) {
        for (pin_index = 0; pin_index < num_of_pins; ++pin_index) {
            const int inet = block[to_block].nets[pin_index];
            if (OPEN == inet || TRUE == net[inet].is_global) {
                continue;
            }
            /* just considering the nets connected to to_block */
            const int net_pin = net_pin_index[to_block][pin_index];
            if (net_pin != 0) {
            /* This net is driving moved to_block! *
             * If this net is being driven by a block that has moved, we do not *
             * need to compute the change in the timing cost (here) since it was *
             * computed in the fanout of the net on the driving block, also    *
             * computing it here would double count the change, and mess up the *
             * delta_timing_cost value */
                if (net[inet].node_block[0] != to_block
                        && net[inet].node_block[0] != from_block) {
                    temp_delay = comp_td_point_to_point_delay_parallel(inet,
                                                                       net_pin,
                                                                       local_block);

                    local_temp_point_to_point_delay_cost[inet][net_pin] = temp_delay;
                    local_temp_point_to_point_timing_cost[inet][net_pin] =
                                timing_place_crit[inet][net_pin] * temp_delay;

                    delta_delay_cost +=
                        (local_temp_point_to_point_delay_cost[inet][net_pin]
                           - point_to_point_delay_cost[inet][net_pin]);
                    delta_timing_cost +=
                        (local_temp_point_to_point_timing_cost[inet][net_pin]
                           - point_to_point_timing_cost[inet][net_pin]);
                }
            } else {  /* net_pin was driver_pin, so to_block was a driver block */
                /*this net is being driven by a moved block, recompute */
                /*all point to point connections on this net. */
                for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                    temp_delay = comp_td_point_to_point_delay_parallel(inet,
                                                                       ipin,
                                                                       local_block);

                    local_temp_point_to_point_delay_cost[inet][ipin] = temp_delay;
                    local_temp_point_to_point_timing_cost[inet][ipin] =
                                timing_place_crit[inet][ipin] * temp_delay;

                    delta_delay_cost +=
                        (local_temp_point_to_point_delay_cost[inet][ipin]
                           - point_to_point_delay_cost[inet][ipin]);
                    delta_timing_cost +=
                        (local_temp_point_to_point_timing_cost[inet][ipin]
                          - point_to_point_timing_cost[inet][ipin]);
                }
            }  /* end of else(net_pin was a driver_pin) */
        }  /* end of for (pin_index = 0; pin_index < num_of_pins; ++pin_index) */
    }  /* end of if (to_block != EMPTY) */

    *delta_timing = delta_timing_cost;
    *delta_delay = delta_delay_cost;
}  /* end of void comp_delta_td_cost_parallel(...) */

/*computes the cost (from scratch) due to the delays and criticalities*
 *on all point to point connections, we define the timing cost of     *
 *each connection as criticality * Tdel */
static unsigned long compute_td_costs_parallel_without_update_crit(double* timing_cost,
                                             double* connection_delay_sum,
                                             int start_net,
                                             int finish_net)
{
    double loc_timing_cost = 0.0;
    double loc_connection_delay_sum = 0.0;
    unsigned long local_work = 0;

    int inet = 0;
    for (inet = start_net; inet < finish_net; ++inet) {
        /* for each net ... */
        if (net[inet].is_global == FALSE) {
            /* Do only if not global. */
            const int sinks = net[inet].num_sinks;
            local_work += sinks;

            int ipin = 0;
            for (ipin = 1; ipin <= net[inet].num_sinks; ++ipin) {
                double temp_delay_cost = comp_td_point_to_point_delay(inet,
                                                                     ipin);
                double temp_timing_cost =
                    temp_delay_cost * timing_place_crit[inet][ipin];

                loc_connection_delay_sum += temp_delay_cost;
                point_to_point_delay_cost[inet][ipin] = temp_delay_cost;
                temp_point_to_point_delay_cost[inet][ipin] = -1; /*undefined */

                point_to_point_timing_cost[inet][ipin] = temp_timing_cost;
                temp_point_to_point_timing_cost[inet][ipin] = -1;   /*undefined */
                loc_timing_cost += temp_timing_cost;
            } /* for(ipin = 1; ipin <= net[inet].num_sinks; ++ipin) */
        }  /* end of if(net[inet].is_global == FALSE) */
    }  /* end of for(inet = start_net; inet < finish_net; ++inet) */

    *timing_cost = loc_timing_cost;
    *connection_delay_sum = loc_connection_delay_sum;

    return local_work;
}  /* end of static unsigned long compute_td_costs_parallel_without_update_crit(double* timing_cost,) */

/*computes the cost(from scratch) due to the delays and criticalities*
 *on all point to point connections, we define the timing cost of     *
 *each connection as criticality * Tdel */
static unsigned long compute_td_costs_parallel_with_update_crit(double* timing_cost,
                                                                double* connection_delay_sum,
                                                                double** net_slack,
                                                                double max_delay,
                                                                double crit_exponent,
                                                                int start_net,
                                                                int finish_net)
{   /* Just like compute_timing_driven_cost in VPR4.3_double */
    double loc_timing_cost = 0.0;
    double loc_connection_delay_sum = 0.0;
    unsigned long local_work = 0;
    int inet = 0;
    /* for each net between [start_net, finish_net)... */
    for (inet = start_net; inet < finish_net; ++inet) {
        if (!net[inet].is_global) {
            /* Do only if not global. */
            const int sinks = net[inet].num_sinks;
            local_work += sinks;

            int ipin = 0;
            for (ipin = 1; ipin <= sinks; ++ipin) {
                double pin_crit = max(1 - net_slack[inet][ipin] / max_delay,
                                     0.0);
                timing_place_crit[inet][ipin] = pow(pin_crit,
                                                    crit_exponent);
                double temp_delay_cost = comp_td_point_to_point_delay(inet,
                                                                     ipin);
                double temp_timing_cost =
                        temp_delay_cost * timing_place_crit[inet][ipin];

                loc_connection_delay_sum += temp_delay_cost;
                point_to_point_delay_cost[inet][ipin] = temp_delay_cost;
                temp_point_to_point_delay_cost[inet][ipin] = -1;

                point_to_point_timing_cost[inet][ipin] = temp_timing_cost;
                temp_point_to_point_timing_cost[inet][ipin] = -1;
                loc_timing_cost += temp_timing_cost;
            } /* end of for (ipin = 1; ipin <= sinks; ++ipin) */
        } /* if (!net[inet].is_global) */
    } /* end of for (inet = start_net; inet < finish_net; ++inet) */

    *timing_cost = loc_timing_cost;
    *connection_delay_sum = loc_connection_delay_sum;

    return local_work;
}  /* end of static unsigned long compute_td_costs_parallel_with_update_crit(double* timing_cost,) */

/* Finds the cost from scratch.  Done only when the placement   *
 * has been radically changed (i.e. after initial placement).   *
 * Otherwise find the cost change incrementally.  If method     *
 * check is NORMAL, we find bounding boxes that are updateable  *
 * for the larger nets.  If method is CHECK, all bounding boxes *
 * are found via the non_updateable_bb routine, to provide a    *
 * cost which can be used to check the correctness of the       *
 * other routine.                                               */
static double comp_bb_cost_parallel(int start_net,
                                    int end_net)
{
    double cost = 0;
    int k = 0;
    for (k = start_net; k < end_net; ++k) {
        /* for each net ... */
        if (net[k].is_global == FALSE) {
            /* Do only if not global. */
            /* Small nets don't use incremental updating on their bounding boxes, *
             * so they can use a fast bounding box calculator.                    */
            get_non_updateable_bb(k, &bb_coords[k]);

            net_cost[k] = get_net_cost(k, &bb_coords[k]);
            cost += net_cost[k];
        }
    }
    return cost;
}  /* end of static double comp_bb_cost_parallel(int start_net,) */

/* This routine finds the bounding box of each net from scratch (i.e.    *
 * from only the block location information).  It updates both the       *
 * coordinate and number of blocks on each tedge information.  It        *
 * should only be called when the bounding box information is not valid. */
static void get_bb_from_scratch_parallel(int inet,
                                         bbox_t* coords,
                                         bbox_t* num_on_edges,
                                         local_block_t* local_block)
{
    /* I need a list of blocks to which this net connects, with no block listed *
     * more than once, in order to get a proper count of the number on the tedge *
     * of the bounding box.                                                     */
    int n_pins = -1;
    int* plist = NULL;
    if (duplicate_pins[inet] == 0) {
        plist = net[inet].node_block;
        n_pins = net[inet].num_sinks + 1;
    } else {
        plist = unique_pin_list[inet];
        n_pins = (net[inet].num_sinks + 1) - duplicate_pins[inet];
    }

    int x = local_block[plist[0]].x;
    int y = local_block[plist[0]].y;
    x = max(min(x, num_grid_columns), 1);
    y = max(min(y, num_grid_rows), 1);

    int xmin = x;
    int ymin = y;
    int xmax = x;
    int ymax = y;
    int xmin_edge = 1;
    int ymin_edge = 1;
    int xmax_edge = 1;
    int ymax_edge = 1;

    int ipin = -1;
    int block_num = 0;
    for (ipin = 1; ipin < n_pins; ++ipin) {
        block_num = plist[ipin];
        x = local_block[block_num].x;
        y = local_block[block_num].y;

        /* Code below counts IO blocks as being within the 1..num_grid_columns,
         * 1..num_grid_rows clb array. *
         * This is because channels do not go out of the 0..num_grid_columns,  *
         * 0..num_grid_rows range, and I always take all channels impinging on *
         * the bounding box to be within that bounding box. Hence, this      *
         * "movement" of IO blocks does not affect the which channels are    *
         * included within the bounding box, and it simplifies the code a lot.*/
        x = max(min(x, num_grid_columns), 1);
        y = max(min(y, num_grid_rows), 1);

        if (x == xmin) {
            ++xmin_edge;
        }
        if (x == xmax) {
            /* Recall that xmin could equal xmax -- don't use else */
            ++xmax_edge;
        } else if (x < xmin) {
            xmin = x;
            xmin_edge = 1;
        } else if (x > xmax) {
            xmax = x;
            xmax_edge = 1;
        }

        if (y == ymin) {
            ++ymin_edge;
        }
        if (y == ymax) {
            ++ymax_edge;
        } else if (y < ymin) {
            ymin = y;
            ymin_edge = 1;
        } else if (y > ymax) {
            ymax = y;
            ymax_edge = 1;
        }
    } /*end of for (ipin = 1; ipin < num_pins; ++ipin) */

    /* Copy the coordinates and number on edges information into the proper *
     * structures.                                                          */
    coords->xmin = xmin;
    coords->xmax = xmax;
    coords->ymin = ymin;
    coords->ymax = ymax;

    num_on_edges->xmin = xmin_edge;
    num_on_edges->xmax = xmax_edge;
    num_on_edges->ymin = ymin_edge;
    num_on_edges->ymax = ymax_edge;
}  /* end of void get_bb_from_scratch_parallel(int inet,...)  */

/* Finds the bounding-box of a net and stores its coordinates in the *
 * bb_coord_new data structure. This routine should only be called   *
 * for small nets, since it does not determine enough information for *
 * the bounding-box to be updated incrementally later.                *
 * Currently assumes channels on both sides of the CLBs forming the  *
 * edges of the bounding box can be used. Essentially, I am assuming *
 * the pins always lie on the outside of the bounding box.           */
static void get_non_updateable_bb_parallel(int inet,
                                    bbox_t* bb_coord_new_ptr,
                                    local_block_t* local_block_ptr)
{
    /* current (x,y) was coordinate of inet's driver block */
    int x = local_block_ptr[net[inet].node_block[0]].x;
    int y = local_block_ptr[net[inet].node_block[0]].y;

    /* [xmin, xmax] and [ymin, ymax] record the bounding-box boundary of inet */
    int xmin = x;
    int ymin = y;
    int xmax = x;
    int ymax = y;

    int k = 0;
    for (k = 1; k < (net[inet].num_sinks + 1); ++k) {
        x = local_block_ptr[net[inet].node_block[k]].x;
        y = local_block_ptr[net[inet].node_block[k]].y;

        if (x < xmin) {
            xmin = x;
        } else if (x > xmax) {
            xmax = x;
        }

        if (y < ymin) {
            ymin = y;
        } else if (y > ymax) {
            ymax = y;
        }
    }

    /* Now I've found the coordinates of the bounding box. There are no  *
     * channels beyond num_grid_columns and num_grid_rows, so I want to clip *
     * to that. As well, since I'll always include the channel immediately *
     * below and the channel immediately to the left of the bounding box, *
     * I want to clip to 1 in both directions as well(since minimum channel *
     * index is 0). See route.c for a channel diagram.                     */
    bb_coord_new_ptr->xmin = max(min(xmin, num_grid_columns), 1);
    bb_coord_new_ptr->ymin = max(min(ymin, num_grid_rows), 1);
    bb_coord_new_ptr->xmax = max(min(xmax, num_grid_columns), 1);
    bb_coord_new_ptr->ymax = max(min(ymax, num_grid_rows), 1);
}  /* end of static void get_non_updateable_bb_parallel(int inet, ) */


/* Updates the bounding-box of a net by storing its coordinates in   *
 * the bb_coord_new data structure and the number of blocks on each  *
 * tedge in the bb_edge_new data structure. This routine should only *
 * be called for large nets(sinks > 4), since it has some overhead   *
 * relative to just doing a brute-force bounding box calculation. The *
 * bounding-box coordinate and tedge information for inet must be     *
 * valid before this routine is called.                               *
 * Currently assumes channels on both sides of the CLBs forming the  *
 * edges of the bounding-box can be used. Essentially, I am assuming *
 * the pins always lie on the outside of the bounding box.           */
static void update_bb_parallel(int inet,
                               bbox_t* local_bb_coord_ptr,
                               bbox_t* local_bb_edge_ptr,
                               bbox_t* bb_coord_new_ptr,
                               bbox_t* bb_edge_new_ptr,
                               int xold, int yold,
                               int xnew, int ynew,
                               local_block_t* local_block)
{
    /* IO blocks are considered to be one cell in for simplicity. */
    xnew = max(min(xnew, num_grid_columns), 1);
    ynew = max(min(ynew, num_grid_rows), 1);
    xold = max(min(xold, num_grid_columns), 1);
    yold = max(min(yold, num_grid_rows), 1);

    /*======   First update x-coordinate, then update y-coordinate like ======*
     *======     x-coordinate                                           ======*/
    /* Check if I can update the bounding box incrementally. */
    if (xnew < xold) {
        /* (1) Move to left. */
        /* Update the xmax fields for coordinates and number of edges first. */
        if (xold == local_bb_coord_ptr[inet].xmax) {
            /* Old position at xmax. */
            if (local_bb_edge_ptr[inet].xmax == 1) {
                get_bb_from_scratch_parallel(inet,
                                             bb_coord_new_ptr,
                                             bb_edge_new_ptr,
                                             local_block);
                return;
            } else {
                bb_edge_new_ptr->xmax = local_bb_edge_ptr[inet].xmax - 1;
                bb_coord_new_ptr->xmax = local_bb_coord_ptr[inet].xmax;
            }
        } else {
            /* Move to left, old postion was not at xmax. */
            bb_coord_new_ptr->xmax = local_bb_coord_ptr[inet].xmax;
            bb_edge_new_ptr->xmax = local_bb_edge_ptr[inet].xmax;
        }

        /* Now do the xmin fields for coordinates and number of edges. */
        if (xnew < local_bb_coord_ptr[inet].xmin) {
            /* Moved past xmin */
            bb_coord_new_ptr->xmin = xnew;
            bb_edge_new_ptr->xmin = 1;
        } else if (xnew == local_bb_coord_ptr[inet].xmin) {
            /* Moved to xmin */
            bb_coord_new_ptr->xmin = xnew;
            bb_edge_new_ptr->xmin = local_bb_edge_ptr[inet].xmin + 1;
        } else {
            /* Xmin unchanged. */
            bb_coord_new_ptr->xmin = local_bb_coord_ptr[inet].xmin;
            bb_edge_new_ptr->xmin = local_bb_edge_ptr[inet].xmin;
        }
        /* End of move to left case. */
    } else if (xnew > xold) {
        /* Move to right. */
        /* Update the xmin fields for coordinates and number of edges first. */
        if (xold == local_bb_coord_ptr[inet].xmin) {
            /* Old position at xmin. */
            if (local_bb_edge_ptr[inet].xmin == 1) {
                get_bb_from_scratch_parallel(inet,
                                             bb_coord_new_ptr,
                                             bb_edge_new_ptr,
                                             local_block);
                return;
            } else {
                bb_edge_new_ptr->xmin = local_bb_edge_ptr[inet].xmin - 1;
                bb_coord_new_ptr->xmin = local_bb_coord_ptr[inet].xmin;
            }
        } else {
            /* Move to right, old position was not at xmin. */
            bb_coord_new_ptr->xmin = local_bb_coord_ptr[inet].xmin;
            bb_edge_new_ptr->xmin = local_bb_edge_ptr[inet].xmin;
        }

        /* Now do the xmax fields for coordinates and number of edges. */
        if (xnew > local_bb_coord_ptr[inet].xmax) {
            /* Moved past xmax. */
            bb_coord_new_ptr->xmax = xnew;
            bb_edge_new_ptr->xmax = 1;
        } else if (xnew == local_bb_coord_ptr[inet].xmax) {
            /* Moved to xmax */
            bb_coord_new_ptr->xmax = xnew;
            bb_edge_new_ptr->xmax = local_bb_edge_ptr[inet].xmax + 1;
        } else {
            /* Xmax unchanged. */
            bb_coord_new_ptr->xmax = local_bb_coord_ptr[inet].xmax;
            bb_edge_new_ptr->xmax = local_bb_edge_ptr[inet].xmax;
        }
    /* End of move to right case. */
    } else {
        /* xnew == xold -- no x motion. */
        bb_coord_new_ptr->xmin = local_bb_coord_ptr[inet].xmin;
        bb_coord_new_ptr->xmax = local_bb_coord_ptr[inet].xmax;
        bb_edge_new_ptr->xmin = local_bb_edge_ptr[inet].xmin;
        bb_edge_new_ptr->xmax = local_bb_edge_ptr[inet].xmax;
    }

    /* Now account for the y-direction motion. */
    if (ynew < yold) {
        /* Move down. */
        /* Update the ymax fields for coordinates and number of edges first. */
        if (yold == local_bb_coord_ptr[inet].ymax) {
            /* Old position at ymax. */
            if (local_bb_edge_ptr[inet].ymax == 1) {
                get_bb_from_scratch_parallel(inet,
                                             bb_coord_new_ptr,
                                             bb_edge_new_ptr,
                                             local_block);
                return;
            } else {
                bb_edge_new_ptr->ymax = local_bb_edge_ptr[inet].ymax - 1;
                bb_coord_new_ptr->ymax = local_bb_coord_ptr[inet].ymax;
            }
        } else {
            /* Move down, old postion was not at ymax. */
            bb_coord_new_ptr->ymax = local_bb_coord_ptr[inet].ymax;
            bb_edge_new_ptr->ymax = local_bb_edge_ptr[inet].ymax;
        }

        /* Now do the ymin fields for coordinates and number of edges. */
        if (ynew < local_bb_coord_ptr[inet].ymin) {
            /* Moved past ymin */
            bb_coord_new_ptr->ymin = ynew;
            bb_edge_new_ptr->ymin = 1;
        } else if (ynew == local_bb_coord_ptr[inet].ymin) {
            /* Moved to ymin */
            bb_coord_new_ptr->ymin = ynew;
            bb_edge_new_ptr->ymin = local_bb_edge_ptr[inet].ymin + 1;
        } else {
            /* ymin unchanged. */
            bb_coord_new_ptr->ymin = local_bb_coord_ptr[inet].ymin;
            bb_edge_new_ptr->ymin = local_bb_edge_ptr[inet].ymin;
        }
    /* End of Move Down case. */
    } else if (ynew > yold) {
        /* Moved up. */
        /* Update the ymin fields for coordinates and number of edges first. */
        if (yold == local_bb_coord_ptr[inet].ymin) {
            /* Old position at ymin. */
            if (local_bb_edge_ptr[inet].ymin == 1) {
                get_bb_from_scratch_parallel(inet,
                                             bb_coord_new_ptr,
                                             bb_edge_new_ptr,
                                             local_block);
                return;
            } else {
                bb_edge_new_ptr->ymin = local_bb_edge_ptr[inet].ymin - 1;
                bb_coord_new_ptr->ymin = local_bb_coord_ptr[inet].ymin;
            }
        } else {
            /* Moved up, old position was not at ymin. */
            bb_coord_new_ptr->ymin = local_bb_coord_ptr[inet].ymin;
            bb_edge_new_ptr->ymin = local_bb_edge_ptr[inet].ymin;
        }

        /* Now do the ymax fields for coordinates and number of edges. */
        if (ynew > local_bb_coord_ptr[inet].ymax) {
            /* Moved past ymax. */
            bb_coord_new_ptr->ymax = ynew;
            bb_edge_new_ptr->ymax = 1;
        } else if (ynew == local_bb_coord_ptr[inet].ymax) {
            /* Moved to ymax */
            bb_coord_new_ptr->ymax = ynew;
            bb_edge_new_ptr->ymax = local_bb_edge_ptr[inet].ymax + 1;
        } else {
            /* ymax unchanged. */
            bb_coord_new_ptr->ymax = local_bb_coord_ptr[inet].ymax;
            bb_edge_new_ptr->ymax = local_bb_edge_ptr[inet].ymax;
        }
        /* End of move up case. */
    } else {
        /* ynew == yold -- no y motion. */
        bb_coord_new_ptr->ymin = local_bb_coord_ptr[inet].ymin;
        bb_coord_new_ptr->ymax = local_bb_coord_ptr[inet].ymax;
        bb_edge_new_ptr->ymin = local_bb_edge_ptr[inet].ymin;
        bb_edge_new_ptr->ymax = local_bb_edge_ptr[inet].ymax;
    }
}  /* end of static void update_bb_parallel(int inet,...) */

/* I insisted that the most important parameter was the kthread_id */
/* try swap a pair of blocks in each Extend-SubRegion by each thread parallely */
static int  try_swap_parallel(const double t,
                              const int  kthread_id,
                              const int  x_from,
                              const int  y_from,
                              const int  z_from,
                              const int  from_block,
                              const int  place_cost_type,
                              const place_algorithm_t place_algorithm,
                              const double timing_tradeoff,
                              double*  cost,
                              double*  bb_cost,
                              double*  timing_cost,
                              double*  delay_cost,
                              double  inverse_prev_bb_cost,
                              double  inverse_prev_timing_cost,
                              double  range_limit,
                              int xMin, int xMax,
                              int yMin, int yMax,
                              thread_local_data_for_swap_t* swap_data_ptr)
{
    /* the flow of try_swap_parallel was that: (1) find to_block, it will swap *
     * with to_block */
    /*constrain the swap region*/
    int limit = (int)range_limit;
    limit = min(10, range_limit);
    xMin = max(xMin, x_from - limit);
    xMax = min(xMax, x_from + limit);
    yMin = max(yMin, y_from - limit);
    yMax = min(yMax, y_from + limit);

    /* First, find a block within the current swap region*/
    double delay_delta_c = 0.0;
    int x_to = 0;
    int y_to = 0;
    find_to_block_parallel(x_from, y_from, block[from_block].type,
                           &x_to, &y_to, kthread_id,
                           xMin, xMax, yMin, yMax);
    /* Ensure that from_block and to_block are same type */
    grid_tile_t** local_grid = swap_data_ptr->m_local_grid;
    assert(local_grid[x_from][y_from].type == local_grid[x_to][y_to].type);

    /* then find a location for to_grid tile */
    int z_to = 0;
    if (local_grid[x_to][y_to].type->capacity > 1) {
        int to_grid_capacity = local_grid[x_to][y_to].type->capacity - 1;
        z_to = my_irand_parallel(to_grid_capacity, kthread_id);
    }

    /* Second, swap the from_block and to_block location */
    /* The Swap a pair of locations */
    local_block_t*  local_block = swap_data_ptr->m_local_block;
    int to_block = EMPTY;
    if (local_grid[x_to][y_to].blocks[z_to] == EMPTY) {
        /* Moving from_block to an empty location */
        to_block = EMPTY;
        local_block[from_block].x = x_to;
        local_block[from_block].y = y_to;
        local_block[from_block].z = z_to;
    } else {
        /* Swapping two non-empty location */
        to_block = local_grid[x_to][y_to].blocks[z_to];
        local_block[to_block].x = x_from;
        local_block[to_block].y = y_from;
        local_block[to_block].z = z_from;

        local_block[from_block].x = x_to;
        local_block[from_block].y = y_to;
        local_block[from_block].z = z_to;
    }

    /*------------------------  Third, compute the swap cost ----------------*/
    /* Change in cost due to this swap. */
    double bb_delta_c = 0;
    int num_of_pins = block[from_block].type->num_pins;
    /* (3.1). I must found out the affected nets, it restored in *
     * int* nets_to_update array. */
    int* nets_to_update = swap_data_ptr->m_nets_to_update;
    int* net_block_moved = swap_data_ptr->m_net_block_moved;
    double* local_temp_net_cost = swap_data_ptr->m_local_temp_net_cost;
    const int knum_nets_affected = find_affected_nets_parallel(nets_to_update,
                                                               net_block_moved,
                                                               from_block,
                                                               to_block,
                                                               num_of_pins,
                                                               local_temp_net_cost);
    /* (3.2) compute wirelength-cost parallel */
    int bb_index = 0;
    int k = -1;
    bbox_t*  bb_coord_new = swap_data_ptr->m_bb_coord_new;
    bbox_t*  bb_edge_new = swap_data_ptr->m_bb_edge_new;

    bbox_t*  local_bb_coord = swap_data_ptr->m_local_bb_coord;
    bbox_t*  local_bb_edge = swap_data_ptr->m_local_bb_edge;

    double*  local_net_cost = swap_data_ptr->m_local_net_cost;
    for (k = 0; k < knum_nets_affected; ++k) {
        int inet = nets_to_update[k];
        if (net_block_moved[k] == FROM_AND_TO) {
            continue;
        }

        if (net[inet].num_sinks < SMALL_NET) {
            get_non_updateable_bb_parallel(inet,
                                           &bb_coord_new[bb_index],
                                           local_block);
        } else {
            if (net_block_moved[k] == FROM) {
                update_bb_parallel(inet,
                                   local_bb_coord,
                                   local_bb_edge,
                                   &bb_coord_new[bb_index],
                                   &bb_edge_new[bb_index],
                                   x_from, y_from,
                                   x_to, y_to,
                                   local_block);
            } else {
                update_bb_parallel(inet,
                                   local_bb_coord,
                                   local_bb_edge,
                                   &bb_coord_new[bb_index],
                                   &bb_edge_new[bb_index],
                                   x_to, y_to,
                                   x_from, y_from,
                                   local_block);
            }
        }  /* end of else(num_sinks >= SMALL_NET) */

        if (place_cost_type != NONLINEAR_CONG) {
            local_temp_net_cost[inet] = get_net_cost(inet,
                                                     &bb_coord_new[bb_index]);
            bb_delta_c += (local_temp_net_cost[inet] - local_net_cost[inet]);
        } else {
            printf("can't do nonlinear_cong\n");
            exit(-1);
        }

        ++bb_index;
    } /* end of for (.k = 0; .k < .num_nets_affected; ++.k) */

    /* (3.3) Calcualte the timing_cost parallel */
    double timing_delta_c = 0.0;
    double delta_c = 0.0;
    double**  local_temp_point_to_point_timing_cost =
        swap_data_ptr->m_local_temp_point_to_point_timing_cost;
    double**  local_temp_point_to_point_delay_cost =
        swap_data_ptr->m_local_temp_point_to_point_delay_cost;
    if (place_algorithm == NET_TIMING_DRIVEN_PLACE ||
            place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
        /* When compute timing_driven_cost(no matteer parallel or sequential), 
         * It should distinguish the driver_pin and not driver_pins */
        comp_delta_td_cost_parallel(from_block,
                                    to_block,
                                    num_of_pins,
                                    &timing_delta_c,
                                    &delay_delta_c,
                                    local_block,
                                    local_temp_point_to_point_timing_cost,
                                    local_temp_point_to_point_delay_cost);

        const double bbox_cost = bb_delta_c * inverse_prev_bb_cost;
        const double timing_cost = timing_delta_c * inverse_prev_timing_cost;
        delta_c = (1 - timing_tradeoff) * bbox_cost + timing_tradeoff * timing_cost;
    } else {
        delta_c = bb_delta_c;
    }

    /* Forth, After calcuate placement cost, then assess swap *
     * 1->move accepted, 0->rejected.                  */
    int keep_switch = assess_swap_parallel(delta_c,
                                           t,
                                           kthread_id);
    if (keep_switch) {
        /* Swap successful! first update cost value. */
        *cost += delta_c;
        *bb_cost +=  bb_delta_c;
        if (place_algorithm == NET_TIMING_DRIVEN_PLACE ||
                place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
            *timing_cost += timing_delta_c;
            *delay_cost += delay_delta_c;
        }

        /* then update affected nets */
        bb_index = 0;
        for (k = 0; k < knum_nets_affected; ++k) {
            int inet = nets_to_update[k];
            if (net_block_moved[k] == FROM_AND_TO) {
                local_temp_net_cost[inet] = -1;
                continue;
            }

            local_bb_coord[inet] = bb_coord_new[bb_index];
            if (net[inet].num_sinks >= SMALL_NET) {
                local_bb_edge[inet] = bb_edge_new[bb_index];
            }
            ++bb_index;

            local_net_cost[inet] = local_temp_net_cost[inet];
            local_temp_net_cost[inet] = -1;
        }

        local_grid[x_to][y_to].blocks[z_to] = from_block;
        local_grid[x_from][y_from].blocks[z_from] = to_block;

        if (EMPTY == to_block) {
            ++(local_grid[x_to][y_to].usage);
            --(local_grid[x_from][y_from].usage);
        }
    } else {
        for (k = 0; k < knum_nets_affected; ++k) {
            int inet = nets_to_update[k];
            local_temp_net_cost[inet] = -1;
        }

        local_block[from_block].x = x_from;
        local_block[from_block].y = y_from;
        local_block[from_block].z = z_from;

        if (to_block != EMPTY) {
            local_block[to_block].x = x_to;
            local_block[to_block].y = y_to;
            local_block[to_block].z = z_to;
        }
    }  /* end of else(swap failure!) */

    return keep_switch;
} /* end of static void try_swap_parallel(double t,) */

/*update the global variables from local variables*/
static void update_from_local_to_global(local_block_t* local_block,
                                        grid_tile_t** local_grid,
                                        int x_start, int x_end,
                                        int y_start, int y_end)
{
    int x, y, z, block_moved;
    for (x = x_start; x < x_end; ++x) {
        for (y = y_start; y < y_end; ++y) {
            if (local_grid[x][y].type != EMPTY_TYPE) {
                if (y_end - y < 5 || y - y_start < 5) {
                    localvert_grid[y][x].usage = local_grid[x][y].usage;
                }

                grid[x][y].usage = local_grid[x][y].usage;
                for (z = 0; z < grid[x][y].type->capacity; ++z) {
                    if (local_grid[x][y].blocks[z] != grid[x][y].blocks[z]) {
                        /*block has been moved*/
                        block_moved = local_grid[x][y].blocks[z];

                        //if (y_end - y < 5 && y_end != num_grid_rows + 1 || y - y_start < 5 && y_start != 0)
                        if (y_end - y < 5 || y - y_start < 5) {
                            localvert_grid[y][x].blocks[z] = block_moved;
                        }

                        grid[x][y].blocks[z] = block_moved;

                        /*if the location becomes empty, don't worry about the rest*/
                        if (block_moved == -1) {
                            continue;
                        }

                        block[block_moved].x = local_block[block_moved].x;
                        block[block_moved].y = local_block[block_moved].y;
                        block[block_moved].z = local_block[block_moved].z;
                    }
                }
            }
        }
    }
}  /* end of void update_from_local_to_global(local_block_t* local_block,) */

static void update_from_global_to_local_hori(local_block_t* local_block,
                                             grid_tile_t** local_grid,
                                             int x_start, int x_end,
                                             int y_start, int y_end)
{
    int x, y, z, block_moved;
    for (x = x_start; x < x_end; ++x) {
        for (y = y_start; y < y_end; ++y) {
            if (grid[x][y].type != EMPTY_TYPE) {
                local_grid[x][y].usage = grid[x][y].usage;

                for (z = 0; z < grid[x][y].type->capacity; ++z) {
                    if (grid[x][y].blocks[z] != local_grid[x][y].blocks[z]) {
                        //block has been moved
                        block_moved = grid[x][y].blocks[z];
                        local_grid[x][y].blocks[z] = block_moved;

                        if (block_moved == -1) {
                            continue;
                        }

                        local_block[block_moved].x = block[block_moved].x;
                        local_block[block_moved].y = block[block_moved].y;
                        local_block[block_moved].z = block[block_moved].z;
                    }
                }
            }
        }
    }
} /* end of static void update_from_global_to_local_hori(local_block_t* local_block,) */

static void update_from_global_to_local_vert(local_block_t* local_block,
                                             grid_tile_t** local_grid,
                                             int x_start, int x_end,
                                             int y_start, int y_end)
{
    int x, y, z, block_moved;
    for (y = y_start; y < y_end; ++y) {
        for (x = x_start; x < x_end; ++x) {
            if (localvert_grid[y][x].type != EMPTY_TYPE) {
                local_grid[x][y].usage = localvert_grid[y][x].usage;

                for (z = 0; z < localvert_grid[y][x].type->capacity; ++z) {
                    if (localvert_grid[y][x].blocks[z] != local_grid[x][y].blocks[z]) {
                        //block has been moved
                        block_moved = localvert_grid[y][x].blocks[z];
                        local_grid[x][y].blocks[z] = block_moved;

                        if (block_moved == -1) {
                            continue;
                        }

                        local_block[block_moved].x = block[block_moved].x;
                        local_block[block_moved].y = block[block_moved].y;
                        local_block[block_moved].z = block[block_moved].z;
                    }
                }
            }
        }
    }
}  /* end of static void update_from_global_to_local_vert(local_block_t* local_block,) */

static void update_from_global_to_local_grid_only(grid_tile_t** local_grid,
                                                  const int x_start,
                                                  const int x_end,
                                                  const int y_start,
                                                  const int y_end)
{
    int x, y, z;
    for (x = x_start; x < x_end; ++x) {
        for (y = y_start; y < y_end; ++y) {
            if (grid[x][y].type != EMPTY_TYPE) {
                local_grid[x][y].usage = grid[x][y].usage;

                for (z = 0; z < grid[x][y].type->capacity; ++z) {
                    if (grid[x][y].blocks[z] != local_grid[x][y].blocks[z]) {
                        //block has been moved
                        local_grid[x][y].blocks[z] = grid[x][y].blocks[z];
                    }
                }
            }
        }
    }
} /* end of void update_from_global_to_local_grid_only(grid_tile_t** local_grid,) */

static double my_difftime2(struct timeval* start,
                           struct timeval* end)
{
    long usec = end->tv_usec - start->tv_usec;
    long sec = end->tv_sec - start->tv_sec;
    if (usec < 0) {
        sec--;
        usec += 1000000;
    }

    double ret = (double)sec;
    ret += (double)usec / 1000000;

    return ret;
}

static void print_grid(void)
{
    FILE* print_grid_ptr = fopen("print_grid.txt", "w");
    int col = 0;
    int row = 0;
    for (col = 0; col <= num_grid_columns + 1; ++col) {
        for (row = 0; row <= num_grid_rows + 1; ++row) {
            grid_tile_t grid_tile = grid[col][row];
            const char* grid_name = grid_tile.type->name;
            const int  grid_capacity = grid_tile.type->capacity;
            const int  grid_height = grid_tile.type->height;
            const int  grid_max_sblks = grid_tile.type->max_subblocks;
            fprintf(print_grid_ptr, "grid[%d][%d] is: %s, capacity: %d, height: %d, max_sblks: %d.\n",
                    col, row, grid_name, grid_capacity, grid_height, grid_max_sblks);
        }
    }
    fclose(print_grid_ptr);
} /* end of print_grid() */

static void tp_data_print_to_screen(int kthread_id,
                             double* timing_cost,
                             double* delay_cost,
                             double* crit_exponent,
                             double max_delay,
                             pthread_data_t* input_args,
                             placer_opts_t placer_opts,
                             double* cost,
                             double* bb_cost,
                             int* success_sum, double* sum_of_squares,
                             int* move_counter, int* tot_iter, double* success_rat,
                             double* av_cost, double* av_bb_cost, double* av_timing_cost, double* av_delay_cost,
                             double* std_dev,
                             double* range_limit, double* final_rlim, double* inverse_delta_rlim,
                             double* t, double* oldt, int* inner_iter_num,
                             double place_delay_value)
{
    if (kthread_id == 0) {
        *success_sum = *(input_args->success_sum);
        *sum_of_squares = *(input_args->sum_of_squares);
        *move_counter = *(input_args->move_lim);

        *tot_iter = *(input_args->tot_iter);
        *tot_iter += *move_counter;

        *success_rat = ((double) * success_sum) / *move_counter;

        *av_cost = *(input_args->av_cost);
        *av_bb_cost = *(input_args->av_bb_cost);
        *av_timing_cost = *(input_args->av_timing_cost);
        *av_delay_cost = *(input_args->av_delay_cost);

        *av_timing_cost = *av_timing_cost + 0.;
        *success_rat = *success_rat + 0.;

        if (success_sum == 0) {
            *av_cost = *cost;
            *av_bb_cost = *bb_cost;
            *av_timing_cost = *timing_cost;
            *av_delay_cost = *delay_cost;
        } else {
            *av_cost /= *success_sum;
            *av_bb_cost /= *success_sum;
            *av_timing_cost /= *success_sum;
            *av_delay_cost /= *success_sum;
        }

        *std_dev = get_std_dev(*success_sum, *sum_of_squares, *av_cost);

        update_rlim(range_limit, *success_rat);

        if (placer_opts.place_algorithm == NET_TIMING_DRIVEN_PLACE ||
                placer_opts.place_algorithm == PATH_TIMING_DRIVEN_PLACE) {
            *crit_exponent = (1 - (*range_limit - *final_rlim) *
                                (*inverse_delta_rlim)) *
                (placer_opts.td_place_exp_last - placer_opts.td_place_exp_first) +
                    placer_opts.td_place_exp_first;
        }


        if (*t != 0) {
            *oldt = *t;     // for finding and printing alpha.
            update_t_parallel(t,
                              inner_iter_num,
                              0,
                              *success_rat,
                              input_args->annealing_sched);
        }

#ifndef SPEC
        printf
        ("%11.5g  %10.6g %11.6g  %11.6g  %11.6g %11.6g %11.4g %9.4g %8.3g  %7.4g  %7.4g  %10d  ",
         *oldt, *av_cost, *av_bb_cost, *av_timing_cost, *av_delay_cost,
         place_delay_value, max_delay, *success_rat, *std_dev, *range_limit,
         *crit_exponent, *tot_iter);
#endif


#ifndef SPEC
        printf("%7.4g\n", *t / *oldt);
#endif

        char msg[BUFSIZE];
        sprintf(msg,
                "Cost: %g  BB Cost %g  TD Cost %g  Temperature: %g  max_delay: %g",
                *cost, *bb_cost, *timing_cost, *t, max_delay);
        update_screen(MINOR, msg, PLACEMENT, FALSE);

        if (*t != 0.0) {
            *(input_args->inner_iter_num) = *inner_iter_num;
            *(input_args->exit) = exit_crit(*t,
                                            *av_cost,
                                            input_args->annealing_sched);
            *(input_args->t) = *t;

            if (*(input_args->exit) == 0) {
                *(input_args->t) = *t;
            }

            else {
                *(input_args->t) = 0.;
            }

            *(input_args->av_cost) = *av_cost;
            *(input_args->av_bb_cost) = *av_bb_cost;
            *(input_args->av_timing_cost) = *av_timing_cost;
            *(input_args->av_delay_cost) = *av_delay_cost;
            *(input_args->std_dev) = *std_dev;
            *(input_args->range_limit ) = *range_limit;
            *(input_args->crit_exponent) = *crit_exponent;
            *(input_args->success_rat) = *success_rat;
            *(input_args->tot_iter) = *tot_iter;
        }

#ifdef VERBOSE
        dump_clbs();
#endif
    }
}

