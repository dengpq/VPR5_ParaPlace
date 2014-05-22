#ifndef PATH_DELAY_PARALLEL_H
#define PATH_DELAY_PARALLEL_H

#include "vpr_types_parallel.h"

double** alloc_and_load_timing_graph(timing_info_t timing_inf,
                                    subblock_data_t subblock_data);

t_linked_int* allocate_and_load_critical_path(void);

void  find_fanin_parallel(int thread_id);

unsigned long load_timing_graph_net_delays_parallel(double** net_delay,
                                                    int start,
                                                    int finish);

void load_timing_graph_net_delays(double** net_delay);

double calc_tnodes_arr_time_parallel(int start_node,
                                    int finish_node,
                                    int ilevel);

void  calc_tnodes_req_time_parallel(double T_cycle,
                                    int start_node,
                                    int finish_node,
                                    int ilevel);

double load_net_slack(double** net_slack,
                     double target_cycle_time);

unsigned long compute_net_slacks_parallel(double** net_slack,
                                          int start_net,
                                          int finish_net);

void free_timing_graph(double** net_slack);

void print_timing_graph(char* fname);

void print_net_slack(char* fname,
                     double** net_slack);

void print_critical_path(char* fname,
                         subblock_data_t subblock_data);

void get_tnode_block_and_output_net(int inode,
                                    int* iblk_ptr,
                                    int* inet_ptr);

void do_constant_net_delay_timing_analysis(timing_info_t timing_inf,
                                           subblock_data_t subblock_data,
                                           double constant_net_delay_value);

#endif

