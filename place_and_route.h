#ifndef PLACE_AND_ROUTE_H
#define PLACE_AND_ROUTE_H

#include "vpr_types_parallel.h"

#define INFINITE -1
#define NOT_FOUND 0

#define WNEED 1
#define WL 2
#define PROC_TIME 3

typedef struct s_fmap_cell {
    int fs;         /* at this fs */
    int fc;         /* at this fc */
    int wneed;          /* need wneed to route */
    int wirelength;     /* corresponding wirelength of successful routing at wneed */
    int proc_time;
    struct s_fmap_cell* next;
} t_fmap_cell;

void place_and_route(operation_types_t operation,
                     placer_opts_t placer_opts,
                     char* place_file,
                     char* net_file,
                     char* arch_file,
                     char* route_file,
                     annealing_sched_t annealing_sched,
                     router_opts_t router_opts,
                     detail_routing_arch_t det_routing_arch,
                     segment_info_t* segment_inf,
                     timing_info_t timing_inf,
                     subblock_data_t* subblock_data_ptr,
                     chan_width_distr_t chan_width_dist);

void init_chan(int cfactor,
               chan_width_distr_t chan_width_dist);

#endif

