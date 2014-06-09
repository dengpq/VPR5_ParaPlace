#ifndef PATH_DELAY2_PARALLEL_H
#define PATH_DELAY2_PARALLEL_H

#include "vpr_types_parallel.h"

/*********** Types and defines used by all path_delay_parallel.h modules ****************/
typedef struct {
    int to_node;
    double Tdel;
} edge_t;

/* to_node: index of node at the sink end of this tedge.                      *
 * Tdel: Tdel to go to to_node along this tedge.                             */
typedef struct {
    edge_t* out_edges;
    int      num_out_edges;
    double   arr_time;
    double   req_time;

    edge_t* in_edges;
    int      num_parents;
} t_tnode;

/* out_edges: [0..num_edges - 1].  Array of the edges leaving this vertexes.    *
 * num_edges: Number of edges leaving this node.                             *
 * arr_time:  Arrival time of the last input signal to this node.               *
 * req_time:  Required arrival time of the last input signal to this node if    *
 *         the critical path is not to be lengthened.                        */


/* Info. below is only used to print out and display the critical path.  It  *
 * gives a mapping from each t_node to what circuit element it represents.   *
 * I put this info in a separate structure to maximize cache effectiveness,  *
 * since it's not used much.                                                 */
typedef enum {
    INPAD_SOURCE,
    INPAD_OPIN,
    OUTPAD_IPIN,
    OUTPAD_SINK,
    FB_IPIN,
    FB_OPIN,
    SUBBLK_IPIN,
    SUBBLK_OPIN,
    FF_SINK,
    FF_SOURCE,
    CONSTANT_GEN_SOURCE
} t_tnode_type;

typedef struct {
    t_tnode_type type;
    short ipin;
    short isubblk;
    int   iblk;
} t_tnode_descript;

/* type:  What is this vertexes? (Pad pin, clb pin, subblock pin, etc.)         *
 * ipin:  Number of the FB or subblock pin this vertexes represents, if        *
 *        applicable.                                                        *
 * isubblk: Number of the subblock this vertexes is part of, if applicable.     *
 * iblk:  Number of the block (FB or PAD) this vertexes is part of.            */

/*************** Variables shared only amongst path_delay_parallel.h modules ************/
extern t_tnode* vertexes;      /* [0..num_of_vertexs - 1] */
extern t_tnode_descript* tnode_descript;    /* [0..num_of_vertexs - 1] */
extern int num_of_vertexs;    /* Number of nodes in the timing graph */


/* [0..num_nets - 1].  Gives the index of the vertexes that drives each net. */
extern int* net_to_driver_tnode;


/* [0..num__tnode_levels - 1].  Count and list of tnodes at each level of    *
 * the timing graph, to make breadth-first searches easier.                  */
extern vector_t* tnodes_at_level;
extern int num_tnode_levels;    /* Number of levels in the timing graph. */


/***************** Subroutines exported by this module ***********************/
int alloc_and_load_timing_graph_levels(void);

void check_timing_graph(int num_const_gen,
                        int num_ff,
                        int num_sinks);

double print_critical_path_node(FILE* fp,
                                t_linked_int* critical_path_node,
                                subblock_data_t subblock_data);

#endif

