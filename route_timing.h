#ifndef ROUTE_TIMING_H
#define ROUTE_TIMING_H

#include "vpr_types_parallel.h"

boolean try_timing_driven_route(router_opts_t router_opts,
                                double** net_slack,
                                double** net_delay,
                                vector_t** fb_opins_used_locally);

boolean timing_driven_route_net(int inet,
                                double pres_fac,
                                double max_criticality,
                                double criticality_exp,
                                double astar_fac,
                                double bend_cost,
                                double* net_slack,
                                double* pin_criticality,
                                int* sink_order,
                                t_rt_node** rt_node_of_sink,
                                double T_crit,
                                double* net_delay);

void alloc_timing_driven_route_structs(double** pin_criticality_ptr,
                                       int** sink_order_ptr,
                                       t_rt_node** * rt_node_of_sink_ptr);

void free_timing_driven_route_structs(double* pin_criticality,
                                      int* sink_order,
                                      t_rt_node** rt_node_of_sink);

#endif

