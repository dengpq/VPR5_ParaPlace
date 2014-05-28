#ifndef GLOBALS_H
#define GLOBALS_H

#include "vpr_types_parallel.h"

extern boolean WMF_DEBUG;
extern char* OutFilePrefix;

extern double grid_logic_tile_area;
extern double ipin_mux_trans_size;

extern int num_nets;
extern net_t* net;

extern int num_blocks;
extern block_t* block;

/* Physical FPGA architecture stuff */
extern int num_grid_columns, num_grid_rows;

/* chan_width_x is the x-directed channel; i.e. between rows */
extern int* chan_width_x, *chan_width_y;    /* numerical form */
extern grid_tile_t** grid;

/* [0..num_nets-1] of linked list start pointers.  Defines the routing.  */
extern trace_t** trace_head, ** trace_tail;

/* Structures to define the routing architecture of the FPGA.           */
extern int num_rr_nodes;
extern rr_node_t* rr_node;  /* [0..num_rr_nodes-1]          */
extern int num_rr_indexed_data;
extern rr_indexed_data_t* rr_indexed_data;  /* [0 .. num_rr_indexed_data-1] */
extern vector_t*** rr_node_indices;
extern int** net_rr_terminals;  /* [0..num_nets-1][0..num_pins-1] */
extern switch_info_t* switch_inf; /* [0..det_routing_arch.num_switch-1] */
extern int** rr_blk_source; /* [0..num_blocks-1][0..num_class-1] */

/* This identifies the block_type_ptr of an IO block */
extern block_type_ptr IO_TYPE;

/* This identifies the block_type_ptr of an Empty block */
extern block_type_ptr EMPTY_TYPE;

/* This identifies the block_type_ptr of the default logic block */
extern block_type_ptr CLB_TYPE;

/* Total number of type_descriptors */
extern int num_types;
extern struct s_type_descriptor* type_descriptors;

#endif

