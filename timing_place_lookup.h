#ifndef TIMING_PLACE_LOOKUP_H
#define TIMING_PLACE_LOOKUP_H

#include "vpr_types_parallel.h"

#define IMPOSSIBLE -1       /*indicator of an array location that    */
/*should never be accessed */

void compute_delay_lookup_tables(router_opts_t router_opts,
                                 detail_routing_arch_t det_routing_arch,
                                 segment_info_t* segment_inf,
                                 timing_info_t timing_inf,
                                 chan_width_distr_t chan_width_dist,
                                 subblock_data_t subblock_data);

void free_place_lookup_structs(void);

extern double** delta_io_to_fb;
extern double** delta_fb_to_fb;
extern double** delta_fb_to_io;
extern double** delta_io_to_io;

#endif

