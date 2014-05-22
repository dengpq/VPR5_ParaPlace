#ifndef RR_GRAPH_SBOX_H
#define RR_GRAPH_SBOX_H

#include "vpr_types_parallel.h"

vector_t get_switch_box_tracks(IN int from_i,
                                    IN int from_j,
                                    IN int from_track,
                                    IN rr_type_t from_type,
                                    IN int to_i,
                                    IN int to_j,
                                    IN rr_type_t to_type,
                                    IN vector_t** *switch_block_conn);

void free_switch_block_conn(vector_t** *switch_block_conn,
                            int nodes_per_chan);

vector_t*** alloc_and_load_switch_block_conn(int nodes_per_chan,
                                                  switch_block_t
                                                  switch_block_type,
                                                  int Fs);

int get_simple_switch_block_track(side_types_t from_side,
                                  side_types_t to_side,
                                  int from_track,
                                  switch_block_t switch_block_type,
                                  int nodes_per_chan);
#endif

