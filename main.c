#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "vpr_types_parallel.h"
#include "graphics.h"
#include "read_netlist.h"
#include "print_netlist.h"
#include "draw.h"
#include "place_and_route.h"
#include "stats.h"
#include "path_delay_parallel.h"
#include "OptionTokens.h"
#include "ReadOptions.h"
#include "xml_arch.h"
#include "SetupVPR.h"
#include "rr_graph.h"
#include <sched.h>
#include "const.h"

#include <time.h>
#include <sys/time.h> //gettimeofday
typedef struct {
    int     secs;
    int     usecs;
} TIME_DIFF;

TIME_DIFF* my_difftime (struct timeval* start, struct timeval* end)
{
    TIME_DIFF* diff = (TIME_DIFF*) malloc ( sizeof (TIME_DIFF) );

    if (start->tv_sec == end->tv_sec) {
        diff->secs = 0;
        diff->usecs = end->tv_usec - start->tv_usec;
    } else {
        diff->usecs = 1000000 - start->tv_usec;
        diff->secs = end->tv_sec - (start->tv_sec + 1);
        diff->usecs += end->tv_usec;

        if (diff->usecs >= 1000000) {
            diff->usecs -= 1000000;
            diff->secs += 1;
        }
    }

    return diff;
}

/******** Global variables ********/
int Fs_seed = -1;
boolean WMF_DEBUG = FALSE;

int W_seed = -1;
int binary_search = -1;
char* OutFilePrefix = NULL;

double grid_logic_tile_area = 0;
double ipin_mux_trans_size = 0;

/******** Netlist to be mapped stuff ********/

int num_nets = 0;
net_t* net = NULL;

int num_blocks = 0;
block_t* blocks = NULL;


/* This identifies the block_type_ptr of an IO block */
int num_types = 0;
struct s_type_descriptor* type_descriptors = NULL;

block_type_ptr IO_TYPE = NULL;
block_type_ptr EMPTY_TYPE = NULL;
block_type_ptr CLB_TYPE = NULL;


/******** Physical architecture stuff ********/
int num_grid_columns = 0;
int num_grid_rows = 0;

/* TRUE if this is a global clb pin -- an input pin to which the netlist can *
 * connect global signals, but which does not connect into the normal        *
 * routing via muxes etc.  Marking pins like this (only clocks in my work)   *
 * stops them from screwing up the input switch pattern in the rr_graph      *
 * generator and from creating extra switches that the area model would      *
 * count.                                                                    */

int* chan_width_x = NULL;   /* [0..num_grid_rows] */
int* chan_width_y = NULL;   /* [0..num_grid_columns] */

/* [0..(num_grid_columns+1)][0..(num_grid_rows+1)] Physical block list */
grid_tile_t** clb_grids = NULL;


/******** Structures defining the routing ********/

/* Linked list start pointers.  Define the routing. */
trace_t** trace_head = NULL; /* [0..(num_nets-1)] */
trace_t** trace_tail = NULL; /* [0..(num_nets-1)] */


/******** Structures defining the FPGA routing architecture ********/

int num_rr_nodes = 0;
rr_node_t* rr_node = NULL;  /* [0..(num_rr_nodes-1)] */
vector_t*** rr_node_indices = NULL;

int num_rr_indexed_data = 0;
rr_indexed_data_t* rr_indexed_data = NULL;  /* [0..(num_rr_indexed_data-1)] */

/* Gives the rr_node indices of net terminals. */

int** net_rr_terminals = NULL;  /* [0..num_nets-1][0..num_pins-1] */

/* Gives information about all the switch types                      *
 * (part of routing architecture, but loaded in read_arch.c          */

switch_info_t* switch_inf = NULL; /* [0..(det_routing_arch.num_switch-1)] */

/* Stores the SOURCE and SINK nodes of all CLBs (not valid for pads).     */

int** rr_blk_source = NULL; /* [0..(num_blocks-1)][0..(num_class-1)] */


/********************** Subroutines local to this module ********************/
static void PrintUsage();
static void PrintTitle();
static void freeArch(t_arch* Arch);



/************************* Subroutine definitions ***************************/
int main(int argc,
         char** argv)
{
    t_options Options;
    t_arch Arch = {0};

    operation_types_t Operation;
    placer_opts_t PlacerOpts;
    annealing_sched_t AnnealSched;
    router_opts_t RouterOpts;
    detail_routing_arch_t RoutingArch;
    segment_info_t* Segments;
    timing_info_t Timing;
    subblock_data_t Subblocks;
    boolean ShowGraphics;
    boolean TimingEnabled;
    int GraphPause;

    TIME_DIFF* diff;
    struct timeval t_start, t_finish;
    clock_t start, finish;


    start = clock();
    gettimeofday(&t_start, NULL);

    /* Print title message */
    PrintTitle();

    /* Print usage message if no args */
    if (argc < 2) {
        PrintUsage();
        exit(1);
    }

    /* Read in available inputs  */
    ReadOptions(argc, argv, &Options);

    /* Determine whether timing is on or off */
    TimingEnabled = IsTimingEnabled(Options);

    /* Use inputs to configure VPR */
    SetupVPR(Options, TimingEnabled, &Arch, &Operation, &PlacerOpts,
             &AnnealSched, &RouterOpts, &RoutingArch, &Segments,
             &Timing, &Subblocks, &ShowGraphics, &GraphPause);

    /* Check inputs are reasonable */
    CheckOptions(Options, TimingEnabled);
    CheckArch(Arch, TimingEnabled);

    /* Verify settings don't conflict or otherwise not make sense */
    CheckSetup(Operation, PlacerOpts, AnnealSched, RouterOpts,
               RoutingArch, Segments, Timing, Subblocks, Arch.Chans);

    /* Output the current settings to console. */
    ShowSetup(Options, Arch, TimingEnabled, Operation, PlacerOpts,
              AnnealSched, RouterOpts, RoutingArch, Segments, Timing,
              Subblocks);

    if (Operation == TIMING_ANALYSIS_ONLY) {
        do_constant_net_delay_timing_analysis(
            Timing, Subblocks, Options.constant_net_delay);
        return 0;
    }

    /* Startup X graphics */
    set_graphics_state(ShowGraphics, GraphPause, RouterOpts.route_type);

    if (ShowGraphics) {
        init_graphics("VPR:  Versatile Place and Route for FPGAs");
        alloc_draw_structs();
    }

    /* Do the actual operation */
    place_and_route(Operation, PlacerOpts, Options.PlaceFile,
                    Options.NetFile, Options.ArchFile, Options.RouteFile,
                    AnnealSched, RouterOpts, RoutingArch,
                    Segments, Timing, &Subblocks, Arch.Chans);

    /* Close down X Display */
    if (ShowGraphics) {
        close_graphics();
    }

    /* free data structures */
    free(Options.PlaceFile);
    free(Options.NetFile);
    free(Options.ArchFile);
    free(Options.RouteFile);

    freeArch(&Arch);

    finish = clock();
    gettimeofday(&t_finish, NULL);

    diff = my_difftime(&t_start, &t_finish);
    printf("entire prog: %d.%d, cpu: %f\n", diff->secs, diff->usecs, (double) (finish - start) / CLOCKS_PER_SEC);
    /* Return 0 to single success to scripts */
    return 0;
}



/* Outputs usage message */
static void PrintUsage()
{
    puts("Usage:  vpr circuit.net fpga.arch placed.out routed.out [Options ...]");
    puts("");
    puts("General Options:  [-nodisp] [-auto <int>] [-route_only]");
    puts("\t[-place_only] [-timing_analyze_only_with_net_delay <double>]");
    puts("\t[-fast] [-full_stats] [-timing_analysis on | off] [-outfile_prefix <string>]");
    puts("");
    puts("Placer Options:");
    puts("\t[-place_algorithm bounding_box | net_timing_driven | path_timing_driven]");
    puts("\t[-init_t <double>] [-exit_t <double>]");
    puts("\t[-alpha_t <double>] [-inner_num <double>] [-seed <int>]");
    puts("\t[-place_cost_exp <double>] [-place_cost_type linear | nonlinear]");
    puts("\t[-place_chan_width <int>] [-num_regions <int>] ");
    puts("\t[-fix_pins random | <file.pads>]");
    puts("\t[-enable_timing_computations on | off]");
    puts("\t[-block_dist <int>]");
    puts("");
    puts("Placement Options Valid Only for Timing-Driven Placement:");
    puts("\t[-timing_tradeoff <double>]");
    puts("\t[-recompute_crit_iter <int>]");
    puts("\t[-inner_loop_recompute_divider <int>]");
    puts("\t[-td_place_exp_first <double>]");
    puts("\t[-td_place_exp_last <double>]");
    puts("");
    puts("Router Options:  [-max_router_iterations <int>] [-bb_factor <int>]");
    puts("\t[-initial_pres_fac <double>] [-pres_fac_mult <double>]");
    puts("\t[-acc_fac <double>] [-first_iter_pres_fac <double>]");
    puts("\t[-bend_cost <double>] [-route_type global | detailed]");
    puts("\t[-verify_binary_search] [-route_chan_width <int>]");
    puts("\t[-router_algorithm breadth_first | timing_driven]");
    puts("\t[-base_cost_type intrinsic_delay | delay_normalized | demand_only]");
    puts("");
    puts("Routing options valid only for timing-driven routing:");
    puts("\t[-astar_fac <double>] [-max_criticality <double>]");
    puts("\t[-criticality_exp <double>]");
    puts("");
}



static void
PrintTitle()
{
    puts("");
    puts("VPR FPGA Placement and Routing.");
    puts("Version: Version 5.0.2");
    puts("Compiled: " __DATE__ ".");
    puts("Original VPR by V. Betz.");
    puts("Timing-driven placement enhancements by A. Marquardt.");
    puts("Single-drivers enhancements by Andy Ye with additions by.");
    puts("Mark Fang, Jason Luu, Ted Campbell");
    puts("Heterogeneous stucture support by Jason Luu and Ted Campbell.");
    puts("This code is licensed only for non-commercial use.");
    puts("");
}

static void freeArch(t_arch* Arch)
{
    int i;

    for (i = 0; i < Arch->num_switches; i++) {
        if (Arch->Switches->name != NULL) {
            free(Arch->Switches[i].name);
        }
    }

    free(Arch->Switches);

    for (i = 0; i < Arch->num_segments; i++) {
        if (Arch->Segments->cb != NULL) {
            free(Arch->Segments[i].cb);
        }

        if (Arch->Segments->sb != NULL) {
            free(Arch->Segments[i].sb);
        }
    }

    free(Arch->Segments);
}
