#ifndef STATS_H
#define STATS_H

#include "vpr_types_parallel.h"

void routing_stats(boolean full_stats,
                   router_types_t route_type,
                   int num_switch,
                   segment_info_t* segment_inf,
                   int num_segment,
                   double R_minW_nmos,
                   double R_minW_pmos,
                   directionality_t directionality,
                   boolean timing_analysis_enabled,
                   double** net_slack,
                   double** net_delay,
                   subblock_data_t subblock_data);

void print_wirelen_prob_dist(void);

void print_lambda(void);

#endif

