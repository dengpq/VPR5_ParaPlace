#ifndef TIMING_PLACE_H
#define TIMING_PLACE_H

#include "vpr_types_parallel.h"

void alloc_lookups_and_criticalities(chan_width_distr_t chan_width_dist,
                                     router_opts_t router_opts,
                                     detail_routing_arch_t det_routing_arch,
                                     segment_info_t* segment_inf,
                                     timing_info_t timing_inf,
                                     subblock_data_t subblock_data,
                                     double***  net_delay,
                                     double***  net_slack);

void load_criticalities_parallel(placer_opts_t placer_opts,
                                 double** net_slack,
                                 double max_delay,
                                 double crit_exponent,
                                 int start,
                                 int finish);

void load_criticalities(placer_opts_t placer_opts,
                        double** net_slack,
                        double max_delay,
                        double crit_exponent);

void free_lookups_and_criticalities(double***  net_delay,
                                    double***  net_slack);

/*void print_sink_delays(char *fname);*/
extern double** timing_place_crit;

#endif

