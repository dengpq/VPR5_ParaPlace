#include <stdio.h>
#include <math.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "path_delay_parallel.h"
#include "path_delay2_parallel.h"
#include "net_delay.h"
#include "timing_place_lookup.h"
#include "timing_place.h"
#include "const.h"

double** timing_place_crit;  /*available externally */

static linked_vptr_t* timing_place_crit_chunk_list_head;
static linked_vptr_t* net_delay_chunk_list_head;


/******** prototypes ******************/
static double** alloc_crit(linked_vptr_t** chunk_list_head_ptr);

static void free_crit(linked_vptr_t** chunk_list_head_ptr);

/**************************************/

static double**
alloc_crit(linked_vptr_t** chunk_list_head_ptr)
{

    /* Allocates space for the timing_place_crit data structure *
     * [0..num_nets-1][1..num_pins-1].  I chunk the data to save space on large    *
     * problems.                                                                   */

    double* tmp_ptr;

    *chunk_list_head_ptr = NULL;
    int   chunk_bytes_avail = 0;
    char* chunk_next_avail_mem = NULL;

    double** local_crit = (double**)my_malloc(num_nets * sizeof(double*));

    int inet;
    for (inet = 0; inet < num_nets; ++inet) {
        tmp_ptr = (double*)my_chunk_malloc((net[inet].num_net_pins) *
                                          sizeof(double),
                                          chunk_list_head_ptr,
                                          &chunk_bytes_avail,
                                          &chunk_next_avail_mem);
        local_crit[inet] = tmp_ptr - 1; /* [1..num_sinks] */
    }

    return (local_crit);
}

/**************************************/
static void free_crit(linked_vptr_t** chunk_list_head_ptr)
{
    free_chunk_memory(*chunk_list_head_ptr);
    *chunk_list_head_ptr = NULL;
}

/**************************************/
void print_sink_delays(char* fname)
{
    int num_at_level, num_edges, inode, ilevel, i;
    FILE* fp = my_fopen(fname, "w");

    for (ilevel = num_tnode_levels - 1; ilevel >= 0; ilevel--) {
        num_at_level = tnodes_at_level[ilevel].nelem;

        for (i = 0; i < num_at_level; i++) {
            inode = tnodes_at_level[ilevel].list[i];
            num_edges = tnode[inode].num_edges;

            if (num_edges == 0) {
                /* sink */
                fprintf(fp, "%g\n", tnode[inode].arr_time);
            }
        }
    }

    fclose(fp);
}

/**************************************/
/*set criticality values, returns the maximum criticality found */
/*assumes that net_slack contains correct values, ie. assumes  *
 *that load_net_slack has been called*/
void load_criticalities(double** net_slack,
                        double max_delay,
                        double crit_exponent)
{
    int inet, ipin;
    double pin_crit;
    for (inet = 0; inet < num_nets; ++inet) {
        if (inet == OPEN || net[inet].is_global) {
            continue;
        }

        const int knum_net_pins = net[inet].num_net_pins;
        for (ipin = 1; ipin <= knum_net_pins; ++ipin) {
            /*clip the criticality to never go negative (could happen */
            /*for a constant generator since it's slack is huge) */
            pin_crit = max(1 - net_slack[inet][ipin] / max_delay, 0.);
            timing_place_crit[inet][ipin] =
                pow(pin_crit, crit_exponent);

        }
    }
}

/**************************************/
void alloc_lookups_and_criticalities(chan_width_distr_t chan_width_dist,
                                     router_opts_t router_opts,
                                     detail_routing_arch_t det_routing_arch,
                                     segment_info_t* segment_inf,
                                     timing_info_t timing_inf,
                                     subblock_data_t subblock_data,
                                     double** *net_delay,
                                     double** *net_slack)
{
    (*net_slack) = alloc_and_load_timing_graph(timing_inf, subblock_data);

    (*net_delay) = alloc_net_delay(&net_delay_chunk_list_head);

    compute_delay_lookup_tables(router_opts, det_routing_arch, segment_inf,
                                timing_inf, chan_width_dist, subblock_data);

    timing_place_crit = alloc_crit(&timing_place_crit_chunk_list_head);
}

/**************************************/
void free_lookups_and_criticalities(double** *net_delay,
                                    double** *net_slack)
{

    free(timing_place_crit);
    free_crit(&timing_place_crit_chunk_list_head);

    free_timing_graph(*net_slack);
    free_net_delay(*net_delay, &net_delay_chunk_list_head);

}

/**************************************/
