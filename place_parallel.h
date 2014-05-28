#ifndef PLACE_H
#define PLACE_H

#include "vpr_types_parallel.h"
#include "mst.h"
#include "const.h"

/*The propose of start_finish_nets is to evenly divide up the work allocated
 * to each processor for:
 *  compute_net_slacks_parallel();
 *  compute_td_costs_parallel_with_update_crit();
 *  comp_bb_cost_parallel();
 *
 * allocating an equal # of nets per partition does not work since the inner loop
 * of each of the functions depend on other parameters , ie - number of edges
 * connected to each net or number of sinks.
 *
 * edges_in_this_partition and sinks_in_this_partition will eveutally balance
 * out for each thread as program progresses to provide an equal distribution
 * of work.  */

struct start_finish_nets {  /* FIXME, important data structs for timing_driven placement */
    int start_edge, finish_edge, start_sinks, finish_sinks;
    int edge_partition_size, sink_partition_size, counter_edge, counter_sink;
    unsigned long edges_in_this_partition;
    unsigned long sinks_in_this_partition;
} start_finish_nets[NUM_OF_THREADS] __attribute__ ((aligned(64)));

/*local block structure.
 * only x,y,z coordinates are needed since the other attributes do not change
 * ie - the nets connected to this block, or the name */
typedef struct s_local_block {
    /* (x,y) was the coordinate of the grid tile */
    short x;
    short y;
    /* z was the capacity of the grid tile */
    short z;
} local_block_t;

/*designed for 64b cachelines*/
typedef struct aligned_bar {
    int volatile arrived;
    int entry;
    int volatile proceed;
} __attribute__((aligned(64))) aligned_bar_t;

typedef struct aligned_neighbor_bar {
    int volatile arrived[4];
} __attribute__((aligned(64))) aligned_neighbor_bar_t;

/* structure to store local data to each thread */
typedef struct pthread_data {
    int thread_id;  /* Most important! which thread */
    /* x_{start|end} is the {starting|ending} column of this thread_region, *
     * y_{start|end} is the {starting|ending} row of this thread_region.    */
    int y_start, y_end, x_start, x_end;
    boolean fixed_pins;
    double* t;
    placer_opts_t placer_opts;
    annealing_sched_t annealing_sched;
    int* success_sum, *move_lim, *tot_iter;
    int* inner_iter_num;
    double* av_cost, *av_bb_cost, *av_timing_cost, *av_delay_cost, *sum_of_squares;
    double* timing_cost, *delay_cost;
    double* inverse_prev_bb_cost;
    double* inverse_prev_timing_cost;
    double* cost, *bb_cost, *success_rat, *std_dev;
    double* max_delay;
    int*    num_connections;
    double* place_delay_value;
    double** net_slack;
    double** net_delay;
    double* crit_exponent;
    int* exit;
    double* range_limit;
} __attribute__((aligned(64))) pthread_data_t;


/* Attention, although this struct' member variable was very similar with  *
 * pthread_data_t, but it restore the local copy of pthread_data_t in each *
 * placing thread. */
typedef struct s_thread_local_common_paras {
    placer_opts_t  local_placer_opts;
    boolean        local_fixed_pins;

    int            local_thread_id;

    int            local_x_start;
    int            local_x_end;
    int            local_y_start;
    int            local_y_end;

    int            local_region_x_boundary[3];
    int            local_region_y_boundary[3];
    /* both temperature and range_limit updating depended on success_ratio */
    double         local_temper;
    double         local_old_temper;

    double         local_range_limit;
    double         local_first_rlim;
    double         local_final_rlim;
    /* inverse_delta_rlim = 1 / (first_rlim - final_rlim) */
    double         local_inverse_delta_rlim;

    /* For placement costs */
    double         local_av_cost;
    double         local_av_bb_cost;
    double         local_av_timing_cost;
    double         local_av_delay_cost;

    double         local_inverse_prev_bb_cost;
    double         local_inverse_prev_timing_cost;

    double         local_total_cost;
    double         local_bb_cost;
    double         local_timing_cost;
    double         local_delay_cost;

    int            local_num_conns;
    /* place_delay_value = delay_cost / num_connections */
    double         local_place_delay_value;
    /* sum_of_squares = total_cost * total_cost */
    double         local_sum_of_squares;

    double**       local_net_slack;
    double**       local_net_delay;
    /* crit_exponent updating was depend on range_limit */
    double         local_crit_exponent;
    double         local_max_delay;

    int            local_success_sum;
    int            local_total_iter;
    /* success_ratio = success_sum / total_iter */
    double         local_success_ratio;
    double         local_std_dev;
} __attribute__((aligned(64))) thread_local_common_paras_t;

typedef struct s_thread_local_data_for_swap {
    grid_tile_t**   m_local_grid;
    local_block_t*  m_local_block;

    bbox_t*         m_bb_coord_new;
    bbox_t*         m_bb_edge_new;
    int*            m_nets_to_update;
    int*            m_net_block_moved;

    bbox_t*         m_local_bb_coord;
    bbox_t*         m_local_bb_edge;

    /* local copy of global variables */
    double*         m_local_net_cost;
    double*         m_local_temp_net_cost;

    double**        m_local_temp_point_to_point_timing_cost;
    double**        m_local_temp_point_to_point_delay_cost;
} __attribute__((aligned(64))) thread_local_data_for_swap_t;

/* local data for inside the try_swap loop
struct swap {
    int x_to, y_to, z_to, to_block;
    int i, k, inet, keep_switch, num_of_pins;
    int num_nets_affected, bb_index;
    double delta_c, bb_delta_c, timing_delta_c, delay_delta_c;
    int iNeighbor;
} __attribute__ ((aligned(64))); */

/* specialized mutex to prevent false sharing*/
typedef struct aligned_mutex {
    pthread_mutex_t mutex;
} __attribute__((aligned(64))) aligned_mutex_t;

/* global data declarations */
aligned_bar_t    thread_barriers[NUM_OF_THREADS];
aligned_mutex_t  global_data_access;
/* In VPR_ParaPlace, partial_result was used for bb_cost and timing_cost,
 * partial_results was used for delay_cost. */
double partial_timing_results[NUM_OF_THREADS];
double partial_delay_results[NUM_OF_THREADS];
double partial_bb_results[NUM_OF_THREADS];

/*calculate wall-time difference*/
typedef struct {
    int     secs;
    int     usecs;
} TIME_DIFF;

void try_place_use_multi_threads(placer_opts_t placer_opts,
                                 annealing_sched_t annealing_sched,
                                 chan_width_distr_t chan_width_dist,
                                 router_opts_t router_opts,
                                 detail_routing_arch_t det_routing_arch,
                                 segment_info_t* segment_inf,
                                 timing_info_t timing_inf,
                                 subblock_data_t* subblock_data_ptr,
                                 t_mst_edge** * mst,
                                 operation_types_t operation);

#endif

