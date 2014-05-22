#ifndef GLOBALS_DECLARE_H
#define GLOBALS_DECLARE_H

#include "vpr_types_parallel.h"

/* Netlist to be placed stuff. */
int num_nets;
int num_blocks;
net_t*   net;
block_t* block;
boolean* is_global;

/* Physical FPGA architecture stuff */
int num_grid_columns;
int num_grid_rows;

/* chan_width_x is the x-directed channel; i.e. between rows */
int*  chan_width_x;
int*  chan_width_y; /* numerical form */
grid_tile_t** grid;

/* [0..num_nets-1] of linked list start pointers.  Defines the routing.  */
trace_t** trace_head, ** trace_tail;


/* Structures to define the routing architecture of the FPGA.  */
int        num_rr_nodes;
rr_node_t* rr_node;     /* [0..num_rr_nodes-1]  */
vector_t***  rr_node_indices;
int        num_rr_indexed_data;
rr_indexed_data_t* rr_indexed_data; /* [0 .. num_rr_indexed_data-1] */
int**      net_rr_terminals;  /* [0..num_nets-1][0..num_pins-1] */
struct     s_switch_inf* switch_inf;    /* [0..det_routing_arch.num_switch-1] */
int**      rr_blk_source;  /* [0..num_blocks-1][0..num_class-1] */

#endif

