#include <assert.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "path_delay2_parallel.h"


/************* Variables (globals) shared by all path_delay_parallel.h modules **********/
t_tnode* vertexes;         /* [0..num_of_vertexs - 1] */
t_tnode_descript* tnode_descript;   /* [0..num_of_vertexs - 1] */
int num_of_vertexs;         /* Number of nodes (pins) in the timing graph */

/* [0..num_nets - 1].  Gives the index of the vertexes that drives each net. */
int* net_to_driver_tnode;


/* [0..num__tnode_levels - 1].  Count and list of tnodes at each level of    *
 * the timing graph, to make breadth-first searches easier.                  */
vector_t* tnodes_at_level;
int num_tnode_levels;       /* Number of levels in the timing graph. */


/******************* Subroutines local to this module ************************/
static int* alloc_and_load_tnode_fanin_and_check_edges(int* num_sinks_ptr);


/************************** Subroutine definitions ***************************/
static int* alloc_and_load_tnode_fanin_and_check_edges(int* num_sinks_ptr)
{
    /* Allocates an array and fills it with the number of in-edges (inputs) to   *
     * each vertexes.  While doing this it also checks that each tedge in the timing *
     * graph points to a valid vertexes. Also counts the number of sinks.           */
    int inode, iedge, to_node, num_out_edges;
    int* tnode_num_fanin;
    edge_t* tedge;
    tnode_num_fanin = (int*)my_calloc(num_of_vertexs, sizeof(int));
    int error = 0;
    int num_sinks = 0;

    for (inode = 0; inode < num_of_vertexs; inode++) {
        num_out_edges = vertexes[inode].num_out_edges;

        if (num_out_edges > 0) {
            tedge = vertexes[inode].out_edges;

            for (iedge = 0; iedge < num_out_edges; iedge++) {
                to_node = tedge[iedge].to_node;

                if (to_node < 0 || to_node >= num_of_vertexs) {
                    printf
                    ("Error in alloc_and_load_tnode_fanin_and_check_edges:\n"
                     "vertexes #%d tedge #%d goes to illegal node #%d.\n",
                     inode, iedge, to_node);
                    error++;
                }

                tnode_num_fanin[to_node]++;
            }
        } else if (num_out_edges == 0) {
            num_sinks++;
        } else {
            printf
            ("Error in alloc_and_load_tnode_fanin_and_check_edges: \n"
             "vertexes #%d has %d edges.\n", inode, num_out_edges);
            error++;
        }
    }

    if (error != 0) {
        printf("Found %d Errors in the timing graph.  Aborting.\n",
               error);
        exit(1);
    }

    *num_sinks_ptr = num_sinks;
    return (tnode_num_fanin);
}


int
alloc_and_load_timing_graph_levels(void)
{
    /* Does a breadth-first search through the timing graph in order to levelize *
     * it.  This allows subsequent breadth-first traversals to be faster. Also   *
     * returns the number of sinks in the graph (nodes with no fanout).          */
    t_linked_int* free_list_head, *nodes_at_level_head;
    int inode, num_at_level, iedge, to_node, num_out_edges, num_sinks,
        num_levels, i;
    edge_t* tedge;
    /* [0..num_of_vertexs-1]. # of in-edges to each vertexes that have not yet been    *
     * seen in this traversal.                                                  */
    int* tnode_fanin_left;
    tnode_fanin_left = alloc_and_load_tnode_fanin_and_check_edges(&num_sinks);
    free_list_head = NULL;
    nodes_at_level_head = NULL;
    /* Very conservative->max number of levels = num_of_vertexs.  Realloc later.  *
     * Temporarily need one extra level on the end because I look at the first  *
     * empty level.                                                             */
    tnodes_at_level = (vector_t*)my_malloc((num_of_vertexs + 1) *
                                                sizeof(vector_t));
    /* Scan through the timing graph, putting all the primary input nodes (no    *
     * fanin) into level 0 of the level structure.                               */
    num_at_level = 0;

    for (inode = 0; inode < num_of_vertexs; inode++) {
        if (tnode_fanin_left[inode] == 0) {
            num_at_level++;
            nodes_at_level_head =
                insert_in_int_list(nodes_at_level_head, inode,
                                   &free_list_head);
        }
    }

    alloc_ivector_and_copy_int_list(&nodes_at_level_head, num_at_level,
                                    &tnodes_at_level[0], &free_list_head);
    num_levels = 0;

    while (num_at_level != 0) {
        /* Until there's nothing in the queue. */
        num_levels++;
        num_at_level = 0;

        for (i = 0; i < tnodes_at_level[num_levels - 1].nelem; i++) {
            inode = tnodes_at_level[num_levels - 1].list[i];
            tedge = vertexes[inode].out_edges;
            num_out_edges = vertexes[inode].num_out_edges;

            for (iedge = 0; iedge < num_out_edges; iedge++) {
                to_node = tedge[iedge].to_node;
                tnode_fanin_left[to_node]--;

                if (tnode_fanin_left[to_node] == 0) {
                    num_at_level++;
                    nodes_at_level_head =
                        insert_in_int_list
                        (nodes_at_level_head, to_node,
                         &free_list_head);
                }
            }
        }

        alloc_ivector_and_copy_int_list(&nodes_at_level_head,
                                        num_at_level,
                                        &tnodes_at_level[num_levels],
                                        &free_list_head);
    }

    tnodes_at_level =
        (vector_t*)my_realloc(tnodes_at_level,
                                   num_levels * sizeof(vector_t));
    num_tnode_levels = num_levels;
    free(tnode_fanin_left);
    free_int_list(&free_list_head);
    return (num_sinks);
}


void
check_timing_graph(int num_const_gen,
                   int num_ff,
                   int num_sinks)
{
    /* Checks the timing graph to see that: (1) all the tnodes have been put    *
     * into some level of the timing graph; (2) the number of primary inputs    *
     * to the timing graph is equal to the number of input pads + the number of *
     * constant generators; and (3) the number of sinks (nodes with no fanout)  *
     * equals the number of output pads + the number of flip flops.             */
    int i, j, num_tnodes_check, ilevel, inet, inode, error, num_p_inputs,
        num_p_outputs;
    error = 0;
    num_tnodes_check = 0;
    num_p_inputs = 0;
    num_p_outputs = 0;

    /* Count I/O input and output pads */
    for (i = 0; i < num_blocks; i++) {
        if (blocks[i].block_type == IO_TYPE) {
            for (j = 0; j < IO_TYPE->num_pins; j++) {
                if (blocks[i].nets[j] != OPEN) {
                    if (IO_TYPE->
                            class_inf[IO_TYPE->pin_class[j]].
                            type == DRIVER) {
                        num_p_inputs++;
                    } else {
                        assert(IO_TYPE->
                               class_inf[IO_TYPE->
                                         pin_class[j]].
                               type == RECEIVER);
                        num_p_outputs++;
                    }
                }
            }
        }
    }

    for (ilevel = 0; ilevel < num_tnode_levels; ilevel++) {
        num_tnodes_check += tnodes_at_level[ilevel].nelem;
    }

    if (num_tnodes_check != num_of_vertexs) {
        printf
        ("Error in check_timing_graph: %d tnodes appear in the vertexes level "
         "structure.  Expected %d.\n", num_tnodes_check, num_of_vertexs);
        printf("Check the netlist for combinational cycles.\n");
        error++;
    }

    if (num_const_gen + num_p_inputs != tnodes_at_level[0].nelem) {
        printf
        ("Error in check_timing_graph: %d tnodes are sources (have no "
         "inputs -- expected %d.\n", tnodes_at_level[0].nelem,
         num_const_gen + num_p_inputs);
        error++;
    }

    if (num_sinks != num_p_outputs + num_ff) {
        printf
        ("Error in check_timing_graph: %d tnodes are sinks (have no "
         "outputs -- expected %d.\n", num_sinks,
         num_ff + num_p_outputs);
        error++;
    }

    for (inet = 0; inet < num_nets; inet++) {
        inode = net_to_driver_tnode[inet];

        if (inode < 0 || inode >= num_of_vertexs) {
            printf("Error in check_timing_graph:\n"
                   "\tdriver of net %d has a vertexes mapping of %d (out of range).\n",
                   inet, inode);
            error++;
        }
    }

    if (error != 0) {
        printf("Found %d Errors in the timing graph.  Aborting.\n",
               error);
        exit(1);
    }
}


double
print_critical_path_node(FILE* fp,
                         t_linked_int* critical_path_node,
                         subblock_data_t subblock_data)
{
    /* Prints one vertexes on the critical path out to fp. Returns the Tdel to    *
     * the next node.                                                           */
    int inode, iblk, ipin, inet, downstream_node;
    t_tnode_type type;
    static char* tnode_type_names[] = { "INPAD_SOURCE", "INPAD_OPIN",
                                        "OUTPAD_IPIN", "OUTPAD_SINK", "CLB_IPIN", "CLB_OPIN",
                                        "SUBBLK_IPIN", "SUBBLK_OPIN", "FF_SINK", "FF_SOURCE",
                                        "CONSTANT_GEN_SOURCE"
                                      };              /* NB: static for speed */
    t_linked_int* next_crit_node;
    double Tdel;
    inode = critical_path_node->data;
    type = tnode_descript[inode].type;
    iblk = tnode_descript[inode].iblk;
    ipin = tnode_descript[inode].ipin;
    fprintf(fp, "Node: %d  %s Block #%d (%s)\n", inode,
            tnode_type_names[type], iblk, blocks[iblk].name);

    if (type != INPAD_SOURCE && type != OUTPAD_SINK && type != FF_SINK &&
            type != FF_SOURCE && type != CONSTANT_GEN_SOURCE) {
        fprintf(fp, "Pin: %d ", ipin);
    }

    if (type == SUBBLK_IPIN || type == SUBBLK_OPIN || type == FF_SINK
            || type == FF_SOURCE || type == CONSTANT_GEN_SOURCE) {
        fprintf(fp, "Subblock #%d subpin #%d",
                tnode_descript[inode].isubblk,
                tnode_descript[inode].ipin);
    }

    if (type != INPAD_SOURCE && type != OUTPAD_SINK) {
        fprintf(fp, "\n");
    }

    fprintf(fp, "arr_time: %g  req_time: %g  ", vertexes[inode].arr_time,
            vertexes[inode].req_time);
    next_crit_node = critical_path_node->next;

    if (next_crit_node != NULL) {
        downstream_node = next_crit_node->data;
        Tdel = vertexes[downstream_node].arr_time - vertexes[inode].arr_time;
        fprintf(fp, "Tdel: %g\n", Tdel);
    } else {
        /* last node, no Tdel. */
        Tdel = 0.;
        fprintf(fp, "\n");
    }

    if (type == FB_OPIN) {
        inet = blocks[iblk].nets[ipin];
        fprintf(fp, "Net to next node: #%d (%s).  Pins on net: %d.\n",
                inet, net[inet].name, (net[inet].num_net_pins + 1));
    }

    if (type == INPAD_OPIN) {
        ipin =
            subblock_data.subblock_inf[iblk][tnode_descript[inode].
                                             isubblk].outputs[ipin];
        inet = blocks[iblk].nets[ipin];
        fprintf(fp, "Net to next node: #%d (%s).  Pins on net: %d.\n",
                inet, net[inet].name, (net[inet].num_net_pins + 1));
    }

    fprintf(fp, "\n");
    return (Tdel);
}
