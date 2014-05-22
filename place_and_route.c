#include <stdio.h>
#include <sys/types.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "vpr_utils.h"
#include "globals.h"
#include "place_and_route.h"
#include "mst.h"
#include "place_parallel.h"
#include "read_place.h"
#include "route_export.h"
#include "draw.h"
#include "stats.h"
#include "check_route.h"
#include "rr_graph.h"
#include "path_delay_parallel.h"
#include "net_delay.h"
#include "timing_place.h"
#include <assert.h>

/******************* Subroutines local to this module ************************/
static int binary_search_place_and_route(placer_opts_t placer_opts,
                                         char* place_file,
                                         char* net_file,
                                         char* arch_file,
                                         char* route_file,
                                         boolean full_stats,
                                         boolean verify_binary_search,
                                         annealing_sched_t
                                         annealing_sched,
                                         router_opts_t router_opts,
                                         detail_routing_arch_t
                                         det_routing_arch,
                                         segment_info_t* segment_inf,
                                         timing_info_t timing_inf,
                                         subblock_data_t* subblock_data_ptr,
                                         chan_width_distr_t chan_width_dist,
                                         t_mst_edge** mst);

static double comp_width(channel_t* chan,
                        double x,
                        double separation);


void post_place_sync(IN int num_blocks,
                     INOUT const block_t block_list[],
                     INOUT subblock_data_t* subblock_data_ptr);

void free_subblock_data(subblock_data_t* subblock_data_ptr);


/************************* Subroutine Definitions ****************************/

void
place_and_route(operation_types_t operation,
                placer_opts_t placer_opts,
                char* place_file,
                char* net_file,
                char* arch_file,
                char* route_file,
                annealing_sched_t annealing_sched,
                router_opts_t router_opts,
                detail_routing_arch_t det_routing_arch,
                segment_info_t* segment_inf,
                timing_info_t timing_inf,
                subblock_data_t* subblock_data_ptr,
                chan_width_distr_t chan_width_dist)
{

    /* This routine controls the overall placement and routing of a circuit. */

    char msg[BUFSIZE];
    int width_fac, inet, i;
    boolean success, Fc_clipped;
    double** net_delay, **net_slack;
    linked_vptr_t* net_delay_chunk_list_head;
    vector_t** fb_opins_used_locally; /* [0..num_blocks-1][0..num_class-1] */
    t_mst_edge** mst = NULL;    /* Make sure mst is never undefined */
    int max_pins_per_fb;

    Fc_clipped = FALSE;

    max_pins_per_fb = 0;

    for (i = 0; i < num_types; i++) {
        if (type_descriptors[i].num_pins > max_pins_per_fb) {
            max_pins_per_fb = type_descriptors[i].num_pins;
        }
    }

    if (placer_opts.place_freq == PLACE_NEVER) {
        /* Read the placement from a file */
        read_place(place_file, net_file, arch_file, num_grid_columns, num_grid_rows, num_blocks,
                   block);
        sync_grid_to_blocks(num_blocks, block, num_grid_columns, num_grid_rows, grid);
    } else {
        assert((PLACE_ONCE == placer_opts.place_freq) ||
               (PLACE_ALWAYS == placer_opts.place_freq));

        try_place_use_multi_threads(placer_opts, annealing_sched, chan_width_dist,
                  router_opts, det_routing_arch, segment_inf,
                  timing_inf, subblock_data_ptr, &mst, operation);
        print_place(place_file, net_file, arch_file);
    }

    post_place_sync(num_blocks, block, subblock_data_ptr);


    fflush(stdout);

    if (operation == PLACE_ONLY) {
        return;
    }

    width_fac = router_opts.fixed_channel_width;

    /* If channel width not fixed, use binary search to find min W */
    if (NO_FIXED_CHANNEL_WIDTH == width_fac) {
        binary_search_place_and_route(placer_opts, place_file,
                                      net_file, arch_file, route_file,
                                      router_opts.full_stats, router_opts.verify_binary_search,
                                      annealing_sched, router_opts,
                                      det_routing_arch, segment_inf,
                                      timing_inf, subblock_data_ptr,
                                      chan_width_dist, mst);
    } else {
        if (det_routing_arch.directionality == UNI_DIRECTIONAL) {
            if (width_fac % 2 != 0) {
                printf
                ("Error: place_and_route.c: given odd chan width (%d) for udsd architecture\n",
                 width_fac);
                exit(1);
            }
        }

        /* Other constraints can be left to rr_graph to check since this is one pass routing */


        /* Allocate the major routing structures. */

        fb_opins_used_locally = alloc_route_structs(*subblock_data_ptr);

        if (timing_inf.timing_analysis_enabled) {
            net_slack =
                alloc_and_load_timing_graph(timing_inf,
                                            *subblock_data_ptr);
            net_delay = alloc_net_delay(&net_delay_chunk_list_head);
        } else {
            net_delay = NULL;   /* Defensive coding. */
            net_slack = NULL;
        }

        success =
            try_route(width_fac, router_opts, det_routing_arch,
                      segment_inf, timing_inf, net_slack, net_delay,
                      chan_width_dist, fb_opins_used_locally, mst,
                      &Fc_clipped);

        if (Fc_clipped) {
            printf
            ("Warning: Fc_output was too high and was clipped to full (maximum) connectivity.\n");
        }

        if (success == FALSE) {
            printf
            ("Circuit is unrouteable with a channel width factor of %d\n\n",
             width_fac);
            sprintf(msg,
                    "Routing failed with a channel width factor of %d.  ILLEGAL routing shown.",
                    width_fac);
        }

        else {
            check_route(router_opts.route_type,
                        det_routing_arch.num_switch,
                        fb_opins_used_locally);
            get_serial_num();

            printf
            ("Circuit successfully routed with a channel width factor of %d.\n\n",
             width_fac);

            routing_stats(router_opts.full_stats, router_opts.route_type,
                          det_routing_arch.num_switch, segment_inf,
                          det_routing_arch.num_segment,
                          det_routing_arch.R_minW_nmos,
                          det_routing_arch.R_minW_pmos,
                          det_routing_arch.directionality,
                          timing_inf.timing_analysis_enabled,
                          net_slack, net_delay, *subblock_data_ptr);

            print_route(route_file);

#ifdef CREATE_ECHO_FILES
            /*print_sink_delays("routing_sink_delays.echo"); */
#endif /* CREATE_ECHO_FILES */

            sprintf(msg,
                    "Routing succeeded with a channel width factor of %d.\n\n",
                    width_fac);
        }

        free_subblock_data(subblock_data_ptr);

        init_draw_coords(max_pins_per_fb);
        update_screen(MAJOR, msg, ROUTING,
                      timing_inf.timing_analysis_enabled);

        if (timing_inf.timing_analysis_enabled) {
            assert(net_slack);
            free_timing_graph(net_slack);

            assert(net_delay);
            free_net_delay(net_delay, &net_delay_chunk_list_head);
        }

        free_route_structs(fb_opins_used_locally);
        fflush(stdout);
    }

    /*WMF: cleaning up memory usage */
    if (mst) {
        for (inet = 0; inet < num_nets; inet++) {
            assert(mst[inet]);
            free(mst[inet]);
        }

        free(mst);
        mst = NULL;
    }

}

/* This routine performs a binary search to find the minimum number of      *
 * tracks per channel required to successfully route a circuit, and returns *
 * that minimum width_fac.                                                  */
static int binary_search_place_and_route(placer_opts_t placer_opts,
                                         char* place_file, char* net_file,
                                         char* arch_file, char* route_file,
                                         boolean full_stats,
                                         boolean verify_binary_search,
                                         annealing_sched_t annealing_sched,
                                         router_opts_t router_opts,
                                         detail_routing_arch_t det_routing_arch,
                                         segment_info_t* segment_inf,
                                         timing_info_t timing_inf,
                                         subblock_data_t* subblock_data_ptr,
                                         chan_width_distr_t chan_width_dist,
                                         t_mst_edge** mst)
{
    boolean success, prev_success, prev2_success, Fc_clipped = FALSE;
    char msg[BUFSIZE];
    linked_vptr_t* net_delay_chunk_list_head;

    /* [0..num_blocks-1][0..num_class-1] */
    vector_t** saved_clb_opins_used_locally;

    int warnings;

    /* Allocate the major routing structures. */
    t_graph_type graph_type;
    if (router_opts.route_type == GLOBAL) {
        graph_type = GRAPH_GLOBAL;
    } else {
        graph_type = (det_routing_arch.directionality ==
                      BI_DIRECTIONAL ? GRAPH_BIDIR : GRAPH_UNIDIR);
    }

    int max_pins_per_fb = 0;
    int i = -1;
    for (i = 0; i < num_types; ++i) {
        max_pins_per_fb = max(max_pins_per_fb,
                              type_descriptors[i].num_pins);
    }

    vector_t** fb_opins_used_locally = alloc_route_structs(*subblock_data_ptr);
    trace_t** best_routing = alloc_saved_routing(fb_opins_used_locally,
                                                        &saved_clb_opins_used_locally);
    
    double** net_delay = NULL;
    double** net_slack = NULL;
    if (timing_inf.timing_analysis_enabled) {
        net_slack = alloc_and_load_timing_graph(timing_inf,
                                                *subblock_data_ptr);
        net_delay = alloc_net_delay(&net_delay_chunk_list_head);
    } else {
        net_delay = NULL;   /* Defensive coding. */
        net_slack = NULL;
    }

    /* UDSD by AY Start */
    int udsd_multiplier = 0;
    if (det_routing_arch.directionality == BI_DIRECTIONAL) {
        udsd_multiplier = 1;
    } else {
        udsd_multiplier = 2;
    }

    /* UDSD by AY End */
    int current = 0;
    int low = 0;
    if (router_opts.fixed_channel_width != NO_FIXED_CHANNEL_WIDTH) {
        current = router_opts.fixed_channel_width + 5 * udsd_multiplier;
        low = router_opts.fixed_channel_width - 1 * udsd_multiplier;
    } else {
        current = max_pins_per_fb + max_pins_per_fb % 2;    /* Binary search part */
        low = -1;
    }

    /* Constraints must be checked to not break rr_graph generator */
    if (det_routing_arch.directionality == UNI_DIRECTIONAL) {
        if (current % 2 != 0) {
            printf("Error: place_and_route.c: tried odd chan width (%d) for udsd architecture\n",
                   current);
            exit(1);
        }
    } else {
        if (det_routing_arch.Fs % 3) {
            printf("Fs must be three in bidirectional mode\n");
            exit(1);
        }
    }

    int high = -1;
    int final = -1;
    int try_w_limit = 0;

    int attempt_count = 0;
    while (final == -1) {
        printf("low, high, current %d %d %d\n", low, high, current);
        fflush(stdout);

        /* Check if the channel width is huge to avoid overflow.  Assume the *
         * circuit is unroutable with the current router options if we're    *
         * going to overflow.                                                */
        if (router_opts.fixed_channel_width != NO_FIXED_CHANNEL_WIDTH) {
            if (current > router_opts.fixed_channel_width * 4) {
                printf("This circuit appears to be unroutable with the current "
                       "router options. Last failed at %d\n", low);
                printf("Aborting routing procedure.\n");
                exit(1);
            }
        } else {
            if (current > 305) {
                printf("This circuit appears to be unroutable with the current "
                       "router options. Last failed at %d\n", low);
                printf("Aborting routing procedure.\n");
                exit(1);
            }
        }

        if ((current * 3) < det_routing_arch.Fs) {
            printf("width factor is now below specified Fs. Stop search.\n");
            final = high;
            break;
        }

        if (placer_opts.place_freq == PLACE_ALWAYS) {
            placer_opts.place_chan_width = current;
            try_place_use_multi_threads(placer_opts, annealing_sched, chan_width_dist,
                      router_opts, det_routing_arch, segment_inf,
                      timing_inf, subblock_data_ptr, &mst, 0);
        }

        success =
            try_route(current, router_opts, det_routing_arch, segment_inf,
                      timing_inf, net_slack, net_delay, chan_width_dist,
                      fb_opins_used_locally, mst, &Fc_clipped);
        attempt_count++;
        fflush(stdout);
#if 1

        if (success && (Fc_clipped == FALSE)) {
#else

        if (success &&
              (Fc_clipped == FALSE || det_routing_arch.Fc_type == FRACTIONAL)) {
#endif
            high = current;
            /* If Fc_output is too high, set to full connectivity but warn the user */
            if (Fc_clipped) {
                printf("Warning: Fc_output was too high and was clipped to full (maximum) connectivity.\n");
            }
            /* If we're re-placing constantly, save placement in case it is best. */
#if 0

            if (placer_opts.place_freq == PLACE_ALWAYS) {
                print_place(place_file, net_file, arch_file);
            }

#endif

            /* Save routing in case it is best. */
            save_routing(best_routing, fb_opins_used_locally,
                         saved_clb_opins_used_locally);
            if ((high - low) <= 1 * udsd_multiplier) {
                final = high;
            }

            if (low != -1) {
                current = (high + low) / 2;
            } else {
                current = high / 2; /* haven't found lower bound yet */
            }
        } else {
            /* last route not successful */
            if (success && Fc_clipped) {
                printf
                ("Routing rejected, Fc_output was too high.\n");
                success = FALSE;
            }

            low = current;

            if (high != -1) {
                if ((high - low) <= 1 * udsd_multiplier) {
                    final = high;
                }
                current = (high + low) / 2;
            } else {
                if (router_opts.fixed_channel_width !=
                        NO_FIXED_CHANNEL_WIDTH) {
                    /* FOR Wneed = f(Fs) search */
                    if (low < router_opts.fixed_channel_width + 30) {
                        current = low + 5 * udsd_multiplier;
                    } else {
                        printf("Aborting: Wneed = f(Fs) search found exceedingly large Wneed (at least %d)\n",
                         low);
                        exit(1);
                    }
                } else {
                    current = low * 2;  /* Haven't found upper bound yet */
                }
            }
        }

        current = current + current % udsd_multiplier;
    }
    /* The binary search above occassionally does not find the minimum    *
     * routeable channel width.  Sometimes a circuit that will not route  *
     * in 19 channels will route in 18, due to router flukiness.  If      *
     * verify_binary_search is set, the code below will ensure that FPGAs *
     * with channel widths of final-2 and final-3 wil not route           *
     * successfully.  If one does route successfully, the router keeps    *
     * trying smaller channel widths until two in a row (e.g. 8 and 9)    *
     * fail.                                                              */
    if (verify_binary_search) {
        printf("\nVerifying that binary search found min. channel width ...\n");
        prev_success = TRUE;    /* Actually final - 1 failed, but this makes router */
        /* try final-2 and final-3 even if both fail: safer */
        prev2_success = TRUE;

        current = final - 2;

        while (prev2_success || prev_success) {
            if ((router_opts.fixed_channel_width !=
                    NO_FIXED_CHANNEL_WIDTH)
                    && (current < router_opts.fixed_channel_width)) {
                break;
            }

            fflush(stdout);

            if (current < 1) {
                break;
            }

            if (placer_opts.place_freq == PLACE_ALWAYS) {
                placer_opts.place_chan_width = current;
                try_place_use_multi_threads(placer_opts, annealing_sched,
                          chan_width_dist, router_opts,
                          det_routing_arch, segment_inf,
                          timing_inf, subblock_data_ptr, &mst, 0);
            }

            success =
                try_route(current, router_opts, det_routing_arch,
                          segment_inf, timing_inf, net_slack,
                          net_delay, chan_width_dist,
                          fb_opins_used_locally, mst, &Fc_clipped);

            if (success && Fc_clipped == FALSE) {
                final = current;
                save_routing(best_routing, fb_opins_used_locally,
                             saved_clb_opins_used_locally);

                if (placer_opts.place_freq == PLACE_ALWAYS) {
                    print_place(place_file, net_file,
                                arch_file);
                }
            }

            prev2_success = prev_success;
            prev_success = success;
            current--;

            if (det_routing_arch.directionality == UNI_DIRECTIONAL) {
                current--; /* width must be even */
            }
        }
    }

    /* End binary search verification. */
    /* Restore the best placement (if necessary), the best routing, and  *
     * * the best channel widths for final drawing and statistics output.  */
    init_chan(final, chan_width_dist);
#if 0

    if (placer_opts.place_freq == PLACE_ALWAYS) {
        printf("Reading best placement back in.\n");
        placer_opts.place_chan_width = final;
        read_place(place_file, net_file, arch_file, placer_opts,
                   router_opts, chan_width_dist, det_routing_arch,
                   segment_inf, timing_inf, subblock_data_ptr, &mst);
    }

#endif
    free_rr_graph();

    build_rr_graph(graph_type,
                   num_types, type_descriptors, num_grid_columns, num_grid_rows, grid,
                   chan_width_x[0], NULL,
                   det_routing_arch.switch_block_type, det_routing_arch.Fs,
                   det_routing_arch.num_segment, det_routing_arch.num_switch, segment_inf,
                   det_routing_arch.global_route_switch,
                   det_routing_arch.delayless_switch, timing_inf,
                   det_routing_arch.wire_to_ipin_switch,
                   router_opts.base_cost_type, &warnings);

    restore_routing(best_routing, fb_opins_used_locally,
                    saved_clb_opins_used_locally);
    check_route(router_opts.route_type, det_routing_arch.num_switch,
                fb_opins_used_locally);
    get_serial_num();

    if (Fc_clipped) {
        printf
        ("Warning: Best routing Fc_output too high, clipped to full (maximum) connectivity.\n");
    }

    printf("Best routing used a channel width factor of %d.\n\n", final);

    routing_stats(full_stats, router_opts.route_type,
                  det_routing_arch.num_switch, segment_inf,
                  det_routing_arch.num_segment,
                  det_routing_arch.R_minW_nmos,
                  det_routing_arch.R_minW_pmos,
                  det_routing_arch.directionality,
                  timing_inf.timing_analysis_enabled, net_slack, net_delay,
                  *subblock_data_ptr);

    print_route(route_file);

#ifdef CREATE_ECHO_FILES
    /* print_sink_delays("routing_sink_delays.echo"); */
#endif /* CREATE_ECHO_FILES */

    init_draw_coords(max_pins_per_fb);
    sprintf(msg, "Routing succeeded with a channel width factor of %d.",
            final);
    update_screen(MAJOR, msg, ROUTING, timing_inf.timing_analysis_enabled);

    free_subblock_data(subblock_data_ptr);

    if (timing_inf.timing_analysis_enabled) {
        free_timing_graph(net_slack);
        free_net_delay(net_delay, &net_delay_chunk_list_head);
    }

    free_route_structs(fb_opins_used_locally);
    free_saved_routing(best_routing, saved_clb_opins_used_locally);
    fflush(stdout);

    return (final);
}


/* Assigns widths to channels (in tracks).  Minimum one track          *
 * per channel.  io channels are io_rat * maximum in interior          *
 * tracks wide.  The channel distributions read from the architecture  *
 * file are scaled by cfactor.                                         */
void init_chan(int cfactor,
               chan_width_distr_t chan_width_dist)
{
    double  chan_width_io = chan_width_dist.chan_width_io;
    channel_t chan_x_dist = chan_width_dist.chan_x_dist;
    channel_t chan_y_dist = chan_width_dist.chan_y_dist;

    /* io channel widths */
    int nio = (int)floor(cfactor * chan_width_io + 0.5);
    if (nio == 0) {
        nio = 1;    /* No zero width channels */
    }

    chan_width_x[0] = chan_width_x[num_grid_rows] = nio;
    chan_width_y[0] = chan_width_y[num_grid_columns] = nio;
    double separation = 0.0;
    double x = 0.0;
    int i = 0;
    if (num_grid_rows > 1) {
        separation = 1.0 / (num_grid_rows - 2.0); /* Norm. distance between two channels. */
        x = 0.0; /* This avoids div by zero if num_grid_rows = 2. */
        chan_width_x[1] = (int)floor(cfactor * comp_width(&chan_x_dist, x,
                                                          separation) + 0.5);
        /* No zero width channels */
        chan_width_x[1] = max(chan_width_x[1], 1);

        for (i = 1; i < num_grid_rows - 1; ++i) {
            x = (double)i / ((double)(num_grid_rows - 2.));
            chan_width_x[i + 1] =
                (int)floor(cfactor *
                           comp_width(&chan_x_dist, x,
                                      separation) + 0.5);
            chan_width_x[i + 1] = max(chan_width_x[i + 1], 1);
        }
    }

    if (num_grid_columns > 1) {
        separation = 1. / (num_grid_columns - 2.);    /* Norm. distance between two channels. */
        x = 0.;     /* Avoids div by zero if num_grid_columns = 2. */
        chan_width_y[1] = (int)floor(cfactor * comp_width(&chan_y_dist, x,
                                                          separation) + 0.5);
        chan_width_y[1] = max(chan_width_y[1], 1);

        for (i = 1; i < num_grid_columns - 1; i++) {
            x = (double)i / ((double)(num_grid_columns - 2.));
            chan_width_y[i + 1] = (int)floor(cfactor * comp_width(&chan_y_dist,
                                                                  x, separation) + 0.5);
            chan_width_y[i + 1] = max(chan_width_y[i + 1], 1);
        }
    }

#ifdef VERBOSE
    printf("\nchan_width_x:\n");

    for (i = 0; i <= num_grid_rows; ++i) {
        printf("%d  ", chan_width_x[i]);
    }
    printf("\n\nchan_width_y:\n");

    for (i = 0; i <= num_grid_columns; ++i) {
        printf("%d  ", chan_width_y[i]);
    }
    printf("\n\n");
#endif
}


static double
comp_width(channel_t* chan,
           double x,
           double separation)
{

    /* Return the relative channel density.  *chan points to a channel   *
     * functional description data structure, and x is the distance      *
     * (between 0 and 1) we are across the chip.  separation is the      *
     * distance between two channels, in the 0 to 1 coordinate system.   */

    double val;

    switch (chan->type) {

        case UNIFORM:
            val = chan->peak;
            break;

        case GAUSSIAN:
            val = (x - chan->xpeak) * (x - chan->xpeak) / (2 * chan->width *
                                                           chan->width);
            val = chan->peak * exp(-val);
            val += chan->dc;
            break;

        case PULSE:
            val = (double)fabs((double)(x - chan->xpeak));

            if (val > chan->width / 2.) {
                val = 0;
            } else {
                val = chan->peak;
            }

            val += chan->dc;
            break;

        case DELTA:
            val = x - chan->xpeak;

            if (val > -separation / 2. && val <= separation / 2.) {
                val = chan->peak;
            } else {
                val = 0.;
            }

            val += chan->dc;
            break;

        default:
            printf("Error in comp_width:  Unknown channel type %d.\n",
                   chan->type);
            exit(1);
            break;
    }

    return (val);
}

/* After placement, logical pins for subblocks, blocks, and nets must be updated to correspond with physical pins of type */
/* This function should only be called once */
void
post_place_sync(IN int num_blocks,
                INOUT const block_t block_list[],
                INOUT subblock_data_t* subblock_data_ptr)
{
    int iblk, j, k, inet;
    block_type_ptr type;
    int max_num_block_pins;
    subblock_t** subblock_inf;

    subblock_inf = subblock_data_ptr->subblock_inf;


    /* Go through each block */
    for (iblk = 0; iblk < num_blocks; ++iblk) {
        type = block[iblk].type;
        assert(type->num_pins % type->capacity == 0);
        max_num_block_pins = type->num_pins / type->capacity;

        /* Logical location and physical location is offset by z * max_num_block_pins */
        /* Sync blocks and nets */
        for (j = 0; j < max_num_block_pins; j++) {
            inet = block[iblk].nets[j];

            if (inet != OPEN && block[iblk].z > 0) {
                assert(block[iblk].
                       nets[j +
                            block[iblk].z * max_num_block_pins] ==
                       OPEN);
                block[iblk].nets[j +
                                 block[iblk].z *
                                 max_num_block_pins] =
                                     block[iblk].nets[j];
                block[iblk].nets[j] = OPEN;

                for (k = 0; k < type->num_pins; k++) {
                    if (net[inet].node_block[k] == iblk) {
                        assert(net[inet].
                               node_block_pin[k] == j);
                        net[inet].node_block_pin[k] =
                            j +
                            block[iblk].z *
                            max_num_block_pins;
                        break;
                    }
                }
            }
        }

        /* Shift subblock pins to point to synced block pins */
        for (j = 0; j < subblock_data_ptr->num_subblocks_per_block[iblk];
                j++) {
            for (k = 0; k < type->max_subblock_inputs; k++) {
                if (subblock_inf[iblk][j].inputs[k] != OPEN) {
                    subblock_inf[iblk][j].inputs[k] =
                        subblock_inf[iblk][j].inputs[k] +
                        block[iblk].z * max_num_block_pins;
                }
            }

            for (k = 0; k < type->max_subblock_outputs; k++) {
                if (subblock_inf[iblk][j].outputs[k] != OPEN) {
                    subblock_inf[iblk][j].outputs[k] =
                        subblock_inf[iblk][j].outputs[k] +
                        block[iblk].z * max_num_block_pins;
                }
            }

            if (subblock_inf[iblk][j].clock != OPEN) {
                subblock_inf[iblk][j].clock =
                    subblock_inf[iblk][j].clock +
                    block[iblk].z * max_num_block_pins;
            }
        }
    }
}

void
free_subblock_data(subblock_data_t* subblock_data_ptr)
{
    int i, j;

    /* Frees all the subblock data structures.                                 */
    for (i = 0; i < num_blocks; i++) {
        if (subblock_data_ptr->subblock_inf[i] != NULL) {
            for (j = 0; j < subblock_data_ptr->num_subblocks_per_block[i]; j++) {
                free(subblock_data_ptr->subblock_inf[i][j].name);
                free(subblock_data_ptr->subblock_inf[i][j].outputs);
                free(subblock_data_ptr->subblock_inf[i][j].inputs);
            }

            free(subblock_data_ptr->subblock_inf[i]);

        }
    }

    free(subblock_data_ptr->subblock_inf);
    free(subblock_data_ptr->num_subblocks_per_block);

    /* Mark as freed. */
    subblock_data_ptr->num_subblocks_per_block = NULL;
    subblock_data_ptr->subblock_inf = NULL;
}
