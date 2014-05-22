#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "OptionTokens.h"
#include "ReadOptions.h"
#include "xml_arch.h"
#include "SetupVPR.h"


void
CheckSetup(IN operation_types_t Operation,
           IN placer_opts_t PlacerOpts,
           IN annealing_sched_t AnnealSched,
           IN router_opts_t RouterOpts,
           IN detail_routing_arch_t RoutingArch,
           IN segment_info_t* Segments,
           IN timing_info_t Timing,
           IN subblock_data_t Subblocks,
           IN chan_width_distr_t Chans)
{
    int i;
    int Tmp;

    if ((NONLINEAR_CONG == PlacerOpts.place_cost_type) &&
            (Operation != PLACE_AND_ROUTE) &&
            (PLACE_ALWAYS == PlacerOpts.place_freq)) {
        printf(ERRTAG "Replacing using the nonlinear congestion option "
               "for each channel width makes sense only for full "
               "place and route.\n");
        exit(1);
    }

    if ((NONLINEAR_CONG == PlacerOpts.place_cost_type) &&
            (BOUNDING_BOX_PLACE != PlacerOpts.place_algorithm)) {
        /* Note that this may work together, but I have not tested it */
        printf(ERRTAG
               "Cannot use non-linear placement only supported with "
               "bounding box placement\n");
        exit(1);
    }

    if ((GLOBAL == RouterOpts.route_type) &&
            (TIMING_DRIVEN == RouterOpts.router_algorithm)) {
        printf(ERRTAG "The global router does not support timing-drvien "
               "routing.\n");
        exit(1);
    }

    if ((GLOBAL == RouterOpts.route_type) &&
            (BOUNDING_BOX_PLACE != PlacerOpts.place_algorithm)) {
        /* Works, but very weird.  Can't optimize timing well, since you're
         * not doing proper architecture Tdel modelling. */
        printf(WARNTAG
               "Using global routing with timing-driven placement. "
               "This is allowed, but strange, and circuit speed will suffer.\n");
    }

    if ((FALSE == Timing.timing_analysis_enabled) &&
            ((PlacerOpts.place_algorithm == NET_TIMING_DRIVEN_PLACE) ||
             (PlacerOpts.place_algorithm == PATH_TIMING_DRIVEN_PLACE))) {
        /* May work, not tested */
        printf(ERRTAG
               "Timing analysis must be enabled for timing-driven placement\n");
        exit(1);
    }

    if ((ROUTE_ONLY == Operation) && (USER == PlacerOpts.pad_loc_type)) {
        printf(ERRTAG "You cannot specify both a full placement file and "
               "a pad location file.\n");
        exit(1);
    }

    if ((ROUTE_ONLY == Operation) || (PLACE_AND_ROUTE == Operation)) {
        if ((TIMING_DRIVEN == RouterOpts.router_algorithm) &&
                (FALSE == Timing.timing_analysis_enabled)) {
            printf(ERRTAG
                   "Cannot perform timing-driven routing when timing "
                   "analysis is disabled.\n");
            exit(1);
        }

        if ((FALSE == Timing.timing_analysis_enabled) &&
                (DEMAND_ONLY != RouterOpts.base_cost_type)) {
            printf(ERRTAG
                   "base_cost_type must be demand_only when timing "
                   "analysis is disabled.\n");
            exit(1);
        }
    }

    if ((TIMING_ANALYSIS_ONLY == Operation) &&
            (FALSE == Timing.timing_analysis_enabled)) {
        printf(ERRTAG
               "-timing_analyze_only_with_net_delay option requires "
               "that timing analysis not be disabled.\n");
        exit(1);
    }

    if ((NONLINEAR_CONG == PlacerOpts.place_cost_type) &&
            ((PlacerOpts.num_regions > num_grid_columns) || (PlacerOpts.num_regions > num_grid_rows))) {
        printf(ERRTAG "Cannot use more regions than clbs in "
               "placement cost function.\n");
        exit(1);
    }

    if (DETAILED == RouterOpts.route_type) {
        if ((Chans.chan_x_dist.type != UNIFORM) ||
                (Chans.chan_y_dist.type != UNIFORM) ||
                (Chans.chan_x_dist.peak != Chans.chan_y_dist.peak) ||
                (Chans.chan_x_dist.peak != Chans.chan_width_io)) {
            printf(ERRTAG "Detailed routing currently only supported "
                   "on FPGAs with all channels of equal width.\n");
            exit(1);
        }
    }

    for (i = 0; i < RoutingArch.num_segment; ++i) {
        Tmp = Segments[i].opin_switch;

        if (FALSE == switch_inf[Tmp].buffered) {
            printf(ERRTAG "opin_switch (#%d) of segment type #%d "
                   "is not buffered.\n", Tmp, i);
            exit(1);
        }
    }

    if (UNI_DIRECTIONAL == RoutingArch.directionality) {
        if ((RouterOpts.fixed_channel_width != NO_FIXED_CHANNEL_WIDTH) &&
                (RouterOpts.fixed_channel_width % 2 > 0)) {
            printf(ERRTAG
                   "Routing channel width must be even for unidirectional\n");
            exit(1);
        }

        if ((PlacerOpts.place_chan_width != NO_FIXED_CHANNEL_WIDTH) &&
                (PlacerOpts.place_chan_width % 2 > 0)) {
            printf(ERRTAG
                   "Place channel width must be even for unidirectional\n");
            exit(1);
        }
    }
}
