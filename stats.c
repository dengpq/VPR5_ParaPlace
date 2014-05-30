#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "rr_graph_area.h"
#include "segment_stats.h"
#include "stats.h"
#include "net_delay.h"
#include "path_delay_parallel.h"


/********************** Subroutines local to this module *********************/

static void load_channel_occupancies(int** chanx_occ,
                                     int** chany_occ);

static void get_num_bends_and_length(int inet,
                                     int* bends,
                                     int* length,
                                     int* segments);

static void get_length_and_bends_stats(void);

static void get_channel_occupancy_stats(void);



/************************* Subroutine definitions ****************************/


void
routing_stats(boolean full_stats,
              router_types_t route_type,
              int num_switch,
              segment_info_t* segment_inf,
              int num_segment,
              double R_minW_nmos,
              double R_minW_pmos,
              directionality_t directionality,
              boolean timing_analysis_enabled,
              double** net_slack,
              double** net_delay,
              subblock_data_t subblock_data)
{
    /* Prints out various statistics about the current routing.  Both a routing *
     * and an rr_graph must exist when you call this routine.                   */
    double T_crit;
    get_length_and_bends_stats();
    get_channel_occupancy_stats();
    printf("Logic Area (in minimum width transistor areas):\n");
    printf("Total Logic Area: %g  Per 1x1 logic tile: %g\n",
           num_grid_columns * num_grid_rows * grid_logic_tile_area, grid_logic_tile_area);

    if (route_type == DETAILED) {
        count_routing_transistors(directionality, num_switch, segment_inf,
                                  R_minW_nmos, R_minW_pmos);
        get_segment_usage_stats(num_segment, segment_inf);

        if (timing_analysis_enabled) {
            load_net_delay_from_routing(net_delay);
#ifdef CREATE_ECHO_FILES
            print_net_delay(net_delay, "net_delay.echo");
#endif /* CREATE_ECHO_FILES */
            load_timing_graph_net_delays(net_delay);
            T_crit = load_net_slack(net_slack, 0);
#ifdef CREATE_ECHO_FILES
            print_timing_graph("timing_graph.echo");
            print_net_slack("net_slack.echo", net_slack);
            print_critical_path("critical_path.echo", subblock_data);
#endif /* CREATE_ECHO_FILES */
            printf("\n");
            printf("Critical Path: %g (s)\n", T_crit);
        }
    }

    if (full_stats == TRUE) {
        print_wirelen_prob_dist();
    }
}


static void
get_length_and_bends_stats(void)
{
    /* Figures out maximum, minimum and average number of bends and net length   *
     * in the routing.                                                           */
    int inet, bends, total_bends, max_bends;
    int length, total_length, max_length;
    int segments, total_segments, max_segments;
    double av_bends, av_length, av_segments;
    int num_global_nets;
    max_bends = 0;
    total_bends = 0;
    max_length = 0;
    total_length = 0;
    max_segments = 0;
    total_segments = 0;
    num_global_nets = 0;

    for (inet = 0; inet < num_nets; inet++) {
        if (net[inet].is_global == FALSE) {
            /* Globals don't count. */
            get_num_bends_and_length(inet, &bends, &length,
                                     &segments);
            total_bends += bends;
            max_bends = max(bends, max_bends);
            total_length += length;
            max_length = max(length, max_length);
            total_segments += segments;
            max_segments = max(segments, max_segments);
        } else {
            num_global_nets++;
        }
    }

    av_bends = (double)total_bends / (double)(num_nets - num_global_nets);
    printf
    ("\nAverage number of bends per net: %#g  Maximum # of bends: %d\n\n",
     av_bends, max_bends);
    av_length = (double)total_length / (double)(num_nets - num_global_nets);
    printf("\nThe number of routed nets (nonglobal): %d\n",
           num_nets - num_global_nets);
    printf("Wirelength results (all in units of 1 clb segments):\n");
    printf("\tTotal wirelength: %d   Average net length: %#g\n",
           total_length, av_length);
    printf("\tMaximum net length: %d\n\n", max_length);
    av_segments = (double)total_segments / (double)(num_nets - num_global_nets);
    printf("Wirelength results in terms of physical segments:\n");
    printf("\tTotal wiring segments used: %d   Av. wire segments per net: "
           "%#g\n", total_segments, av_segments);
    printf("\tMaximum segments used by a net: %d\n\n", max_segments);
}


static void
get_channel_occupancy_stats(void)
{
    /* Determines how many tracks are used in each channel.                    */
    int i, j, max_occ, total_x, total_y;
    double av_occ;
    int** chanx_occ;        /* [1..num_grid_columns][0..num_grid_rows] */
    int** chany_occ;        /* [0..num_grid_columns][1..num_grid_rows] */
    chanx_occ = (int**)alloc_matrix(1, num_grid_columns, 0, num_grid_rows, sizeof(int));
    chany_occ = (int**)alloc_matrix(0, num_grid_columns, 1, num_grid_rows, sizeof(int));
    load_channel_occupancies(chanx_occ, chany_occ);
    printf("\nX - Directed channels:\n\n");
    printf("j\tmax occ\tav_occ\t\tcapacity\n");
    total_x = 0;

    for (j = 0; j <= num_grid_rows; j++) {
        total_x += chan_width_x[j];
        av_occ = 0.;
        max_occ = -1;

        for (i = 1; i <= num_grid_columns; i++) {
            max_occ = max(chanx_occ[i][j], max_occ);
            av_occ += chanx_occ[i][j];
        }

        av_occ /= num_grid_columns;
        printf("%d\t%d\t%-#9g\t%d\n", j, max_occ, av_occ,
               chan_width_x[j]);
    }

    printf("\nY - Directed channels:\n\n");
    printf("i\tmax occ\tav_occ\t\tcapacity\n");
    total_y = 0;

    for (i = 0; i <= num_grid_columns; i++) {
        total_y += chan_width_y[i];
        av_occ = 0.;
        max_occ = -1;

        for (j = 1; j <= num_grid_rows; j++) {
            max_occ = max(chany_occ[i][j], max_occ);
            av_occ += chany_occ[i][j];
        }

        av_occ /= num_grid_rows;
        printf("%d\t%d\t%-#9g\t%d\n", i, max_occ, av_occ,
               chan_width_y[i]);
    }

    printf("\nTotal Tracks in X-direction: %d  in Y-direction: %d\n\n",
           total_x, total_y);
    free_matrix(chanx_occ, 1, num_grid_columns, 0, sizeof(int));
    free_matrix(chany_occ, 0, num_grid_columns, 1, sizeof(int));
}


static void
load_channel_occupancies(int** chanx_occ,
                         int** chany_occ)
{
    /* Loads the two arrays passed in with the total occupancy at each of the  *
     * channel segments in the FPGA.                                           */
    int i, j, inode, inet;
    trace_t* tptr;
    rr_type_t rr_type;

    /* First set the occupancy of everything to zero. */

    for (i = 1; i <= num_grid_columns; i++)
        for (j = 0; j <= num_grid_rows; j++) {
            chanx_occ[i][j] = 0;
        }

    for (i = 0; i <= num_grid_columns; i++)
        for (j = 1; j <= num_grid_rows; j++) {
            chany_occ[i][j] = 0;
        }

    /* Now go through each net and count the tracks and pins used everywhere */

    for (inet = 0; inet < num_nets; inet++) {
        if (net[inet].is_global) { /* Skip global nets. */
            continue;
        }

        tptr = trace_head[inet];

        while (tptr != NULL) {
            inode = tptr->index;
            rr_type = rr_node[inode].type;

            if (rr_type == SINK) {
                tptr = tptr->next;  /* Skip next segment. */

                if (tptr == NULL) {
                    break;
                }
            } else if (rr_type == CHANX) {
                j = rr_node[inode].ylow;

                for (i = rr_node[inode].xlow;
                        i <= rr_node[inode].xhigh; i++) {
                    chanx_occ[i][j]++;
                }
            } else if (rr_type == CHANY) {
                i = rr_node[inode].xlow;

                for (j = rr_node[inode].ylow;
                        j <= rr_node[inode].yhigh; j++) {
                    chany_occ[i][j]++;
                }
            }

            tptr = tptr->next;
        }
    }
}


static void
get_num_bends_and_length(int inet,
                         int* bends_ptr,
                         int* len_ptr,
                         int* segments_ptr)
{
    /* Counts and returns the number of bends, wirelength, and number of routing *
     * resource segments in net inet's routing.                                  */
    trace_t* tptr, *prevptr;
    int inode;
    rr_type_t curr_type, prev_type;
    int bends, length, segments;
    bends = 0;
    length = 0;
    segments = 0;
    prevptr = trace_head[inet]; /* Should always be SOURCE. */

    if (prevptr == NULL) {
        printf
        ("Error in get_num_bends_and_length:  net #%d has no traceback.\n",
         inet);
        exit(1);
    }

    inode = prevptr->index;
    prev_type = rr_node[inode].type;
    tptr = prevptr->next;

    while (tptr != NULL) {
        inode = tptr->index;
        curr_type = rr_node[inode].type;

        if (curr_type == SINK) {
            /* Starting a new segment */
            tptr = tptr->next;  /* Link to existing path - don't add to len. */

            if (tptr == NULL) {
                break;
            }

            curr_type = rr_node[tptr->index].type;
        } else if (curr_type == CHANX || curr_type == CHANY) {
            segments++;
            length += 1 + rr_node[inode].xhigh - rr_node[inode].xlow +
                      rr_node[inode].yhigh - rr_node[inode].ylow;

            if (curr_type != prev_type
                    && (prev_type == CHANX || prev_type == CHANY)) {
                bends++;
            }
        }

        prev_type = curr_type;
        tptr = tptr->next;
    }

    *bends_ptr = bends;
    *len_ptr = length;
    *segments_ptr = segments;
}


void print_wirelen_prob_dist(void)
{
    /* Prints out the probability distribution of the wirelength / number   *
     * input pins on a net -- i.e. simulates 2-point net length probability *
     * distribution.                                                        */
    double two_point_length;
    int inet, bends, length, segments, index;
    double av_length;
    int i, incr;
    int     prob_dist_size = num_grid_columns + num_grid_rows + 10;
    double* prob_dist = (double*)my_calloc(prob_dist_size, sizeof(double));
    double norm_fac = 0.0;

    for (inet = 0; inet < num_nets; ++inet) {
        if (net[inet].is_global == FALSE) {
            get_num_bends_and_length(inet, &bends, &length,
                                     &segments);
            /*  Assign probability to two integer lengths proportionately -- i.e.  *
             *  if two_point_length = 1.9, add 0.9 of the pins to prob_dist[2] and *
             *  only 0.1 to prob_dist[1].                                          */
            two_point_length =
                (double)length / (double)(net[inet].num_net_pins);
            index = (int)two_point_length;

            if (index >= prob_dist_size) {
                printf
                ("Warning: index (%d) to prob_dist exceeds its allocated size (%d)\n",
                 index, prob_dist_size);
                printf
                ("Realloc'ing to increase 2-pin wirelen prob distribution array\n");
                incr = index - prob_dist_size + 2;
                prob_dist_size += incr;
                prob_dist = (double*)my_realloc(prob_dist,
                                                prob_dist_size * sizeof(double));

                for (i = prob_dist_size - incr; i < prob_dist_size; i++) {
                    prob_dist[i] = 0.0;
                }
            }

            prob_dist[index] += (net[inet].num_net_pins) * (1 - two_point_length +
                                         index);
            index++;

            if (index >= prob_dist_size) {
                printf
                ("Warning: index (%d) to prob_dist exceeds its allocated size (%d)\n",
                 index, prob_dist_size);
                printf
                ("Realloc'ing to increase 2-pin wirelen prob distribution array\n");
                incr = index - prob_dist_size + 2;
                prob_dist_size += incr;
                prob_dist = (double*)my_realloc(prob_dist,
                                                prob_dist_size * sizeof(double));

                for (i = prob_dist_size - incr; i < prob_dist_size;
                        i++) {
                    prob_dist[i] = 0.0;
                }
            }

            prob_dist[index] += (net[inet].num_net_pins) * (1 - index +
                                                             two_point_length);
            norm_fac += net[inet].num_net_pins;
        }
    }

    /* Normalize so total probability is 1 and print out. */
    printf("\nProbability distribution of 2-pin net lengths:\n\n");
    printf("Length    p(Lenth)\n");
    av_length = 0;

    for (index = 0; index < prob_dist_size; index++) {
        prob_dist[index] /= norm_fac;
        printf("%6d  %10.6f\n", index, prob_dist[index]);
        av_length += prob_dist[index] * index;
    }

    printf("\nThe number of 2-pin nets is ;%g;\n", norm_fac);
    printf("\nExpected value of 2-pin net length (R) is ;%g;\n", av_length);
    printf("\nTotal wire length is ;%g;\n", norm_fac * av_length);
    free(prob_dist);
}


void
print_lambda(void)
{
    /* Finds the average number of input pins used per fb.  Does not    *
     * count inputs which are hooked to global nets (i.e. the clock     *
     * when it is marked global).                                       */
    int block_num, ipin;
    int num_inputs_used = 0;
    int iclass, inet;
    double lambda;
    block_type_ptr type;

    for (block_num = 0; block_num < num_blocks; block_num++) {
        type = block[block_num].type;
        assert(type != NULL);

        if (type != IO_TYPE) {
            for (ipin = 0; ipin < type->num_pins; ipin++) {
                iclass = type->pin_class[ipin];

                if (type->class_inf[iclass].type == RECEIVER) {
                    inet = block[block_num].nets[ipin];

                    if (inet != OPEN)   /* Pin is connected? */
                        if (net[inet].is_global == FALSE) { /* Not a global clock */
                            num_inputs_used++;
                        }
                }
            }
        }
    }

    lambda = (double)num_inputs_used / (double)num_blocks;
    printf("Average lambda (input pins used per fb) is: %g\n", lambda);
}
