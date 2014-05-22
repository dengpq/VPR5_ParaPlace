#ifndef RR_GRAPH2_H
#define RR_GRAPH2_H

#include "vpr_types_parallel.h"

/************** Global variables shared only by the rr_* modules. ************/
extern boolean* rr_edge_done;   /* [0..num_rr_nodes-1].  Used to keep track  *
                 * of whether or not a node has been put in  *
                 * an tedge list yet.                         */

/******************* Subroutines exported by rr_graph2.c *********************/
vector_t*** alloc_and_load_rr_node_indices(IN int nodes_per_chan,
                                                IN int num_grid_columns,
                                                IN int num_grid_rows,
                                                INOUT int* index,
                                                IN segment_details_t *
                                                seg_details);

void free_rr_node_indices(IN vector_t** * rr_node_indices);

int get_rr_node_index(int x,
                      int y,
                      rr_type_t rr_type,
                      int ptc,
                      vector_t** * rr_node_indices);

void free_seg_details(segment_details_t* seg_details,
                      int nodes_per_chan);

segment_details_t* alloc_and_load_seg_details(INOUT int* nodes_per_chan,
                                          IN int max_len,
                                          IN int num_seg_types,
                                          IN segment_info_t* segment_inf,
                                          IN boolean use_full_seg_groups,
                                          IN boolean is_global_graph,
                                          IN directionality_t
                                          directionality);

void dump_seg_details(segment_details_t* seg_details,
                      int nodes_per_chan,
                      char* fname);

int get_seg_start(IN segment_details_t* seg_details,
                  IN int itrack,
                  IN int chan_num,
                  IN int seg_num);

int get_seg_end(IN segment_details_t* seg_details,
                IN int itrack,
                IN int istart,
                IN int chan_num,
                IN int seg_max);

boolean is_cbox(IN int chan,
                IN int seg,
                IN int track,
                IN segment_details_t* seg_details,
                IN directionality_t directionality);

boolean is_sbox(IN int chan,
                IN int wire_seg,
                IN int sb_seg,
                IN int track,
                IN segment_details_t* seg_details,
                IN directionality_t directionality);

int get_bidir_opin_connections(IN int i,
                               IN int j,
                               IN int ipin,
                               IN struct s_linked_edge** edge_list,
                               IN int**** *opin_to_track_map,
                               IN int Fc,
                               IN boolean* rr_edge_done,
                               IN vector_t** * rr_node_indices,
                               IN segment_details_t* seg_details);

int get_unidir_opin_connections(IN int chan,
                                IN int seg,
                                IN int Fc,
                                IN rr_type_t chan_type,
                                IN segment_details_t* seg_details,
                                INOUT t_linked_edge** edge_list_ptr,
                                INOUT int** Fc_ofs,
                                INOUT boolean* rr_edge_done,
                                IN int max_len,
                                IN int nodes_per_chan,
                                IN vector_t** * rr_node_indices,
                                OUT boolean* Fc_clipped);

int get_track_to_ipins(int seg,
                       int chan,
                       int track,
                       t_linked_edge** edge_list_ptr,
                       vector_t** * rr_node_indices,
                       vector_t**** track_to_ipin_lookup,
                       segment_details_t* seg_details,
                       enum e_rr_type chan_type,
                       int chan_length,
                       int wire_to_ipin_switch,
                       directionality_t directionality);

int get_track_to_tracks(IN int from_chan,
                        IN int from_seg,
                        IN int from_track,
                        IN rr_type_t from_type,
                        IN int to_seg,
                        IN rr_type_t to_type,
                        IN int chan_len,
                        IN int nodes_per_chan,
                        IN int* opin_mux_size,
                        IN int Fs_per_side,
                        IN short**** *sblock_pattern,
                        INOUT struct s_linked_edge** edge_list,
                        IN segment_details_t* seg_details,
                        IN directionality_t directionality,
                        IN vector_t** * rr_node_indices,
                        INOUT boolean* rr_edge_done,
                        IN vector_t** *switch_block_conn);

short***** alloc_sblock_pattern_lookup(IN int num_grid_columns,
                                       IN int num_grid_rows,
                                       IN int nodes_per_chan);
void free_sblock_pattern_lookup(INOUT short**** *sblock_pattern);
void load_sblock_pattern_lookup(IN int i,
                                IN int j,
                                IN int nodes_per_chan,
                                IN segment_details_t* seg_details,
                                IN int Fs,
                                IN switch_block_t switch_block_type,
                                INOUT short**** *sblock_pattern);

#endif

