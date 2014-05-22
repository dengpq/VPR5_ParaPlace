#ifndef RR_GRAPH_INDEX_DATA_H
#define RR_GRAPH_INDEX_DATA_H

#include "vpr_types_parallel.h"

void alloc_and_load_rr_indexed_data(segment_info_t* segment_inf,
                                    int num_segment,
                                    vector_t** * rr_node_indices,
                                    int nodes_per_chan,
                                    int wire_to_ipin_switch,
                                    router_base_cost_t base_cost_type);

#endif

