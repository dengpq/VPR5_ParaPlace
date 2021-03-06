#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "rr_graph.h"
#include "check_rr_graph.h"


/********************** Local defines and types *****************************/

#define BUF_FLAG 1
#define PTRANS_FLAG 2
#define BUF_AND_PTRANS_FLAG 3


/*********************** Subroutines local to this module *******************/

static boolean rr_node_is_global_clb_ipin(int inode);

static void check_pass_transistors(int from_node);


/************************ Subroutine definitions ****************************/

void
check_rr_graph(IN t_graph_type graph_type,
               IN int num_types,
               IN block_type_ptr types,
               IN int num_grid_columns,
               IN int num_grid_rows,
               IN grid_tile_t** grid,
               IN int nodes_per_chan,
               IN int Fs,
               IN int num_seg_types,
               IN int num_switches,
               IN segment_info_t* segment_inf,
               IN int global_route_switch,
               IN int delayless_switch,
               IN int wire_to_ipin_switch,
               segment_details_t* seg_details,
               int* Fc_in,
               int* Fc_out,
               vector_t** * rr_node_indices,
               int**** *opin_to_track_map,
               int**** *ipin_to_track_map,
               vector_t**** track_to_ipin_lookup,
               vector_t** * switch_block_conn,
               boolean* perturb_ipins)
{
    int* num_edges_from_current_to_node;    /* [0..num_rr_nodes-1] */
    int* total_edges_to_node;   /* [0..num_rr_nodes-1] */
    char* switch_types_from_current_to_node;    /* [0..num_rr_nodes-1] */
    int inode, iedge, to_node, num_edges;
    short switch_type;
    rr_type_t rr_type, to_rr_type;
    router_types_t route_type;
    boolean is_fringe_warning_sent;
    route_type = DETAILED;

    if (graph_type == GRAPH_GLOBAL) {
        route_type = GLOBAL;
    }

    total_edges_to_node = (int*)my_calloc(num_rr_nodes, sizeof(int));
    num_edges_from_current_to_node = (int*)my_calloc(num_rr_nodes,
                                                     sizeof(int));
    switch_types_from_current_to_node = (char*)my_calloc(num_rr_nodes,
                                                         sizeof(char));

    for (inode = 0; inode < num_rr_nodes; inode++) {
        rr_type = rr_node[inode].type;
        num_edges = rr_node[inode].num_edges;
        check_node(inode, route_type);

        /* Check all the connectivity (edges, etc.) information.                    */

        for (iedge = 0; iedge < num_edges; iedge++) {
            to_node = rr_node[inode].edges[iedge];

            if (to_node < 0 || to_node >= num_rr_nodes) {
                printf
                ("Error in check_rr_graph:  node %d has an tedge %d.\n"
                 "Edge is out of range.\n", inode, to_node);
                exit(1);
            }

            num_edges_from_current_to_node[to_node]++;
            total_edges_to_node[to_node]++;
            switch_type = rr_node[inode].switches[iedge];

            if (switch_type < 0 || switch_type >= num_switches) {
                printf
                ("Error in check_rr_graph:  node %d has a switch type %d.\n"
                 "Switch type is out of range.\n", inode,
                 switch_type);
                exit(1);
            }

            if (switch_inf[switch_type].buffered)
                switch_types_from_current_to_node[to_node] |=
                    BUF_FLAG;
            else
                switch_types_from_current_to_node[to_node] |=
                    PTRANS_FLAG;
        }       /* End for all edges of node. */

        for (iedge = 0; iedge < num_edges; iedge++) {
            to_node = rr_node[inode].edges[iedge];

            if (num_edges_from_current_to_node[to_node] > 1) {
                to_rr_type = rr_node[to_node].type;

                if ((to_rr_type != CHANX && to_rr_type != CHANY) ||
                        (rr_type != CHANX && rr_type != CHANY)) {
                    printf
                    ("Error in check_rr_graph:  node %d connects to node %d "
                     "%d times.\n", inode, to_node,
                     num_edges_from_current_to_node
                     [to_node]);
                    exit(1);
                }
                /* Between two wire segments.  Two connections are legal only if  *
                 * one connection is a buffer and the other is a pass transistor. */
                else if (num_edges_from_current_to_node[to_node] !=
                         2
                         ||
                         switch_types_from_current_to_node[to_node]
                         != BUF_AND_PTRANS_FLAG) {
                    printf
                    ("Error in check_rr_graph:  node %d connects to node %d "
                     "%d times.\n", inode, to_node,
                     num_edges_from_current_to_node
                     [to_node]);
                    exit(1);
                }
            }

            num_edges_from_current_to_node[to_node] = 0;
            switch_types_from_current_to_node[to_node] = 0;
        }

        /* Slow test below.  Leave commented out most of the time. */
#ifdef DEBUG
        check_pass_transistors(inode);
#endif
    }           /* End for all rr_nodes */

    /* I built a list of how many edges went to everything in the code above -- *
     * now I check that everything is reachable.                                */
    is_fringe_warning_sent = FALSE;

    for (inode = 0; inode < num_rr_nodes; inode++) {
        rr_type = rr_node[inode].type;

        if (rr_type != SOURCE) {
            if (total_edges_to_node[inode] < 1 &&
                    !rr_node_is_global_clb_ipin(inode)) {
                boolean is_fringe;
                boolean is_wire;
                /* A global FB input pin will not have any edges, and neither will  *
                 * a SOURCE.  Anything else is an error.                             */
                is_fringe = ((rr_node[inode].xlow == 1) || (rr_node[inode].ylow == 1)
                             || (rr_node[inode].xhigh == num_grid_columns) || (rr_node[inode].yhigh == num_grid_rows));
                is_wire = (rr_node[inode].type == CHANX || rr_node[inode].type == CHANY);

                if (!is_fringe && !is_wire) {
                    printf("Error in check_rr_graph:  node %d has no fanin.\n", inode);
                    exit(1);
                } else if (!is_fringe_warning_sent) {
                    printf("WARNING: in check_rr_graph:  fringe node %d has no fanin.\n"
                           "This is possible on the fringe for low Fc_out, N, and certain Lengths\n"
                           , inode);
                    is_fringe_warning_sent = TRUE;
                }
            }
        } else {
            /* SOURCE.  No fanin for now; change if feedthroughs allowed. */
            if (total_edges_to_node[inode] != 0) {
                printf
                ("Error in check_rr_graph:  SOURCE node %d has a fanin\n"
                 "\tof %d, expected 0.\n", inode,
                 total_edges_to_node[inode]);
                exit(1);
            }
        }
    }

    free(num_edges_from_current_to_node);
    free(total_edges_to_node);
    free(switch_types_from_current_to_node);
}


static boolean rr_node_is_global_clb_ipin(int inode)
{
    if (rr_node[inode].type != IPIN) {
        return (FALSE);
    }
    /* Returns TRUE if inode refers to a global FB input pin node.   */
    block_type_ptr type =
        bin_grids[rr_node[inode].xlow][rr_node[inode].ylow].grid_type;


    int ipin = rr_node[inode].ptc_num;
    return (type->is_global_pin[ipin]);
}


void check_node(int inode,
                router_types_t route_type)
{
    /* This routine checks that the rr_node is inside the grid and has a valid  *
     * pin number, etc.                                                         */
    int xlow, ylow, xhigh, yhigh, ptc_num, capacity;
    rr_type_t rr_type;
    block_type_ptr type;
    int nodes_per_chan, tracks_per_node, num_edges, cost_index;
    double C, R;
    rr_type = rr_node[inode].type;
    xlow = rr_node[inode].xlow;
    xhigh = rr_node[inode].xhigh;
    ylow = rr_node[inode].ylow;
    yhigh = rr_node[inode].yhigh;
    ptc_num = rr_node[inode].ptc_num;
    capacity = rr_node[inode].capacity;
    type = NULL;

    if (xlow > xhigh || ylow > yhigh) {
        printf
        ("Error in check_node:  rr endpoints are (%d,%d) and (%d,%d).\n",
         xlow, ylow, xhigh, yhigh);
        exit(1);
    }

    if (xlow < 0 || xhigh > num_grid_columns + 1 || ylow < 0 || yhigh > num_grid_rows + 1) {
        printf
        ("Error in check_node:  rr endpoints, (%d,%d) and (%d,%d), \n"
         "are out of range.\n", xlow, ylow, xhigh, yhigh);
        exit(1);
    }

    if (ptc_num < 0) {
        printf("Error in check_node.  Inode %d (type %d) had a ptc_num\n"
               "of %d.\n", inode, rr_type, ptc_num);
        exit(1);
    }

    /* Check that the segment is within the array and such. */

    switch (rr_type) {
        case SOURCE:
        case SINK:
        case IPIN:
        case OPIN:
            /* This is used later as well */
            type = bin_grids[xlow][ylow].grid_type;

            if (type == NULL) {
                printf
                ("Error in check_node:  Node %d (type %d) is at an illegal\n"
                 " clb location (%d, %d).\n", inode, rr_type, xlow,
                 ylow);
                exit(1);
            }

            if (xlow != xhigh || ylow != (yhigh - type->height + 1)) {
                printf
                ("Error in check_node:  Node %d (type %d) has endpoints of\n"
                 "(%d,%d) and (%d,%d)\n", inode, rr_type, xlow, ylow,
                 xhigh, yhigh);
                exit(1);
            }

            break;

        case CHANX:
            if (xlow < 1 || xhigh > num_grid_columns || yhigh > num_grid_rows || yhigh != ylow) {
                printf("Error in check_node:  CHANX out of range.\n");
                printf("Endpoints: (%d,%d) and (%d,%d)\n", xlow, ylow,
                       xhigh, yhigh);
                exit(1);
            }

            if (route_type == GLOBAL && xlow != xhigh) {
                printf
                ("Error in check_node:  node %d spans multiple channel segments\n"
                 "which is not allowed with global routing.\n",
                 inode);
                exit(1);
            }

            break;

        case CHANY:
            if (xhigh > num_grid_columns || ylow < 1 || yhigh > num_grid_rows || xlow != xhigh) {
                printf("Error in check_node:  CHANY out of range.\n");
                printf("Endpoints: (%d,%d) and (%d,%d)\n", xlow, ylow,
                       xhigh, yhigh);
                exit(1);
            }

            if (route_type == GLOBAL && ylow != yhigh) {
                printf
                ("Error in check_node:  node %d spans multiple channel segments\n"
                 "which is not allowed with global routing.\n",
                 inode);
                exit(1);
            }

            break;

        default:
            printf("Error in check_node:  Unexpected segment type: %d\n",
                   rr_type);
            exit(1);
    }

    /* Check that it's capacities and such make sense. */

    switch (rr_type) {
        case SOURCE:
            if (ptc_num >= type->num_class
                    || type->class_inf[ptc_num].type != DRIVER) {
                printf
                ("Error in check_node.  Inode %d (type %d) had a ptc_num\n"
                 "of %d.\n", inode, rr_type, ptc_num);
                exit(1);
            }

            if (type->class_inf[ptc_num].num_pins != capacity) {
                printf
                ("Error in check_node.  Inode %d (type %d) had a capacity\n"
                 "of %d.\n", inode, rr_type, capacity);
                exit(1);
            }

            break;

        case SINK:
            if (ptc_num >= type->num_class
                    || type->class_inf[ptc_num].type != RECEIVER) {
                printf
                ("Error in check_node.  Inode %d (type %d) had a ptc_num\n"
                 "of %d.\n", inode, rr_type, ptc_num);
                exit(1);
            }

            if (type->class_inf[ptc_num].num_pins != capacity) {
                printf
                ("Error in check_node.  Inode %d (type %d) has a capacity\n"
                 "of %d.\n", inode, rr_type, capacity);
                exit(1);
            }

            break;

        case OPIN:
            if (ptc_num >= type->num_type_pins
                    || type->class_inf[type->pin_class[ptc_num]].type != DRIVER) {
                printf
                ("Error in check_node.  Inode %d (type %d) had a ptc_num\n"
                 "of %d.\n", inode, rr_type, ptc_num);
                exit(1);
            }

            if (capacity != 1) {
                printf
                ("Error in check_node:  Inode %d (type %d) has a capacity\n"
                 "of %d.\n", inode, rr_type, capacity);
                exit(1);
            }

            break;

        case IPIN:
            if (ptc_num >= type->num_type_pins
                    || type->class_inf[type->pin_class[ptc_num]].type != RECEIVER) {
                printf
                ("Error in check_node.  Inode %d (type %d) had a ptc_num\n"
                 "of %d.\n", inode, rr_type, ptc_num);
                exit(1);
            }

            if (capacity != 1) {
                printf
                ("Error in check_node:  Inode %d (type %d) has a capacity\n"
                 "of %d.\n", inode, rr_type, capacity);
                exit(1);
            }

            break;

        case CHANX:
            if (route_type == DETAILED) {
                nodes_per_chan = chan_width_x[ylow];
                tracks_per_node = 1;
            } else {
                nodes_per_chan = 1;
                tracks_per_node = chan_width_x[ylow];
            }

            if (ptc_num >= nodes_per_chan) {
                printf
                ("Error in check_node:  Inode %d (type %d) has a ptc_num\n"
                 "of %d.\n", inode, rr_type, ptc_num);
                exit(1);
            }

            if (capacity != tracks_per_node) {
                printf
                ("Error in check_node:  Inode %d (type %d) has a capacity\n"
                 "of %d.\n", inode, rr_type, capacity);
                exit(1);
            }

            break;

        case CHANY:
            if (route_type == DETAILED) {
                nodes_per_chan = chan_width_y[xlow];
                tracks_per_node = 1;
            } else {
                nodes_per_chan = 1;
                tracks_per_node = chan_width_y[xlow];
            }

            if (ptc_num >= nodes_per_chan) {
                printf
                ("Error in check_node:  Inode %d (type %d) has a ptc_num\n"
                 "of %d.\n", inode, rr_type, ptc_num);
                exit(1);
            }

            if (capacity != tracks_per_node) {
                printf
                ("Error in check_node:  Inode %d (type %d) has a capacity\n"
                 "of %d.\n", inode, rr_type, capacity);
                exit(1);
            }

            break;

        default:
            printf("Error in check_node:  Unexpected segment type: %d\n",
                   rr_type);
            exit(1);
    }

    /* Check that the number of (out) edges is reasonable. */
    num_edges = rr_node[inode].num_edges;

    if (rr_type != SINK) {
        if (num_edges <= 0) {
            printf("Error: in check_node: node %d has no edges.\n",
                   inode);
            exit(1);
        }
    } else {
        /* SINK -- remove this check if feedthroughs allowed */
        if (num_edges != 0) {
            printf("Error in check_node: node %d is a sink, but has "
                   "%d edges.\n", inode, num_edges);
            exit(1);
        }
    }

    /* Check that the capacitance, resistance and cost_index are reasonable. */
    C = rr_node[inode].C;
    R = rr_node[inode].R;

    if (rr_type == CHANX || rr_type == CHANY) {
        if (C < 0. || R < 0.) {
            printf
            ("Error in check_node: node %d of type %d has R = %g "
             "and C = %g.\n", inode, rr_type, R, C);
            exit(1);
        }
    } else {
        if (C != 0. || R != 0.) {
            printf
            ("Error in check_node: node %d of type %d has R = %g "
             "and C = %g.\n", inode, rr_type, R, C);
            exit(1);
        }
    }

    cost_index = rr_node[inode].cost_index;

    if (cost_index < 0 || cost_index >= num_rr_indexed_data) {
        printf("Error in check_node:  node %d cost index (%d) is out of "
               "range.\n", inode, cost_index);
        exit(1);
    }
}


static void
check_pass_transistors(int from_node)
{
    /* This routine checks that all pass transistors in the routing truly are  *
     * bidirectional.  It may be a slow check, so don't use it all the time.   */
    int from_edge, to_node, to_edge, from_num_edges, to_num_edges;
    rr_type_t from_rr_type, to_rr_type;
    short from_switch_type;
    boolean trans_matched;
    from_rr_type = rr_node[from_node].type;

    if (from_rr_type != CHANX && from_rr_type != CHANY) {
        return;
    }

    from_num_edges = rr_node[from_node].num_edges;

    for (from_edge = 0; from_edge < from_num_edges; from_edge++) {
        to_node = rr_node[from_node].edges[from_edge];
        to_rr_type = rr_node[to_node].type;

        if (to_rr_type != CHANX && to_rr_type != CHANY) {
            continue;
        }

        from_switch_type = rr_node[from_node].switches[from_edge];

        if (switch_inf[from_switch_type].buffered) {
            continue;
        }

        /* We know that we have a pass transitor from from_node to to_node.  Now *
         * check that there is a corresponding tedge from to_node back to         *
         * from_node.                                                            */
        to_num_edges = rr_node[to_node].num_edges;
        trans_matched = FALSE;

        for (to_edge = 0; to_edge < to_num_edges; to_edge++) {
            if (rr_node[to_node].edges[to_edge] == from_node &&
                    rr_node[to_node].switches[to_edge] == from_switch_type) {
                trans_matched = TRUE;
                break;
            }
        }

        if (trans_matched == FALSE) {
            printf
            ("Error in check_pass_transistors:  Connection from node %d to\n"
             "node %d uses a pass transistor (switch type %d), but there is\n"
             "no corresponding pass transistor tedge in the other direction.\n",
             from_node, to_node, from_switch_type);
            exit(1);
        }
    }           /* End for all from_node edges */
}
