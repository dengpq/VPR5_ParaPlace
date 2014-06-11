#include <stdio.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "route_common.h"
#include "place_and_route.h"
#include "route_tree_timing.h"
#include "route_timing.h"
#include "timing_place_lookup.h"
#include "rr_graph.h"
#include "mst.h"
#include "route_export.h"
#include <assert.h>


/*this file contains routines that generate the array containing*/
/*the delays between blocks, this is used in the timing driven  */
/*placement routines */

/*To compute Tdel between blocks we place temporary blocks at */
/*different locations in the FPGA and route nets  between      */
/*the blocks.  From this procedure we generate a lookup table  */
/*which tells us the Tdel between different locations in      */
/*the FPGA */

/*Note: these routines assume that there is a uniform and even */
/*distribution of the different wire segments. If this is not  */
/*the case, then this lookup table will be off */

/*Note: This code removes all heterogeneous types and creates an
artificial 1x1 tile.  A good lookup for heterogeniety
requires more research */

#define NET_COUNT 1     /*we only use one net in these routines,   */
/*it is repeatedly routed and ripped up    */
/*to compute delays between different      */
/*locations, this value should not change  */
#define NET_USED 0      /*we use net at location zero of the net    */
/*structure                                 */
#define NET_USED_SOURCE_BLOCK 0 /*net.blocks[0] is source blocks */
#define NET_USED_SINK_BLOCK 1   /*net.blocks[1] is sink blocks */
#define SOURCE_BLOCK 0      /*blocks[0] is source */
#define SINK_BLOCK 1        /*blocks[1] is sink */

#define BLOCK_COUNT 2       /*use 2 blocks to compute Tdel between  */
/*the various FPGA locations             */
/*do not change this number unless you   */
/*really know what you are doing, it is  */
/*assumed that the net only connects to  */
/*two blocks */

#define NUM_TYPES_USED 3    /* number of types used in look up */

#define DEBUG_TIMING_PLACE_LOOKUP   /*initialize arrays to known state */

#define DUMPFILE "lookup_dump.echo"
/* #define PRINT_ARRAYS *//*only used during debugging, calls routine to  */
/*print out the various lookup arrays           */

/***variables that are exported to other modules***/

/*the delta arrays are used to contain the best case routing Tdel */
/*between different locations on the FPGA. */
double** delta_inpad_to_clb;
double** delta_clb_to_clb;
double** delta_clb_to_outpad;
double** delta_inpad_to_outpad;


/*** Other Global Arrays ******/
/* I could have allocated these as local variables, and passed them all */
/* around, but was too lazy, since this is a small file, it should not  */
/* be a big problem */

static double** net_delay;
static double** net_slack;
static double* pin_criticality;
static int* sink_order;
static t_rt_node** rt_node_of_sink;
static block_type_ptr IO_TYPE_BACKUP;
static block_type_ptr EMPTY_TYPE_BACKUP;
static block_type_ptr FILL_TYPE_BACKUP;
static type_descriptor_t dummy_type_descriptors[NUM_TYPES_USED];
static type_descriptor_t* type_descriptors_backup;
static grid_tile_t** grid_backup;
static int num_types_backup;

static vector_t** clb_opins_used_locally;

#ifdef PRINT_ARRAYS
static FILE* lookup_dump;   /* If debugging mode is on, print out to
                 * the file defined in DUMPFILE */
#endif /* PRINT_ARRAYS */

/*** Function Prototypes *****/

static void alloc_net(void);

static void alloc_block(void);

static void load_simplified_device(void);
static void restore_original_device(void);

static void alloc_and_assign_internal_structures(net_t** original_net,
                                                 block_t
                                                 ** original_block,
                                                 int* original_num_nets,
                                                 int* original_num_blocks);

static void free_and_reset_internal_structures(net_t* original_net,
                                               block_t* original_block,
                                               int original_num_nets,
                                               int original_num_blocks);

static void setup_chan_width(router_opts_t router_opts,
                             chan_width_distr_t chan_width_dist);

static void alloc_routing_structs(router_opts_t router_opts,
                                  detail_routing_arch_t det_routing_arch,
                                  segment_info_t* segment_inf,
                                  timing_info_t timing_inf,
                                  subblock_data_t subblock_data);

static void free_routing_structs(router_opts_t router_opts,
                                 detail_routing_arch_t det_routing_arch,
                                 segment_info_t* segment_inf,
                                 timing_info_t timing_inf);

static void assign_locations(block_type_ptr source_type,
                             int source_x_loc,
                             int source_y_loc,
                             int source_z_loc,
                             block_type_ptr sink_type,
                             int sink_x_loc,
                             int sink_y_loc,
                             int sink_z_loc);

static double assign_blocks_and_route_net(block_type_ptr source_type,
                                         int source_x_loc,
                                         int source_y_loc,
                                         block_type_ptr sink_type,
                                         int sink_x_loc,
                                         int sink_y_loc,
                                         router_opts_t router_opts);

static void alloc_delta_arrays(void);

static void free_delta_arrays(void);

static void generic_compute_matrix(double** *matrix_ptr,
                                   block_type_ptr source_type,
                                   block_type_ptr sink_type,
                                   int source_x,
                                   int source_y,
                                   int start_x,
                                   int end_x,
                                   int start_y,
                                   int end_y,
                                   router_opts_t router_opts,
                                   detail_routing_arch_t det_routing_arch,
                                   segment_info_t* segment_inf,
                                   timing_info_t timing_inf);

static void compute_delta_fb_to_fb(router_opts_t router_opts,
                                   int longest_length);

static void compute_delta_io_to_fb(router_opts_t router_opts,
                                   detail_routing_arch_t det_routing_arch,
                                   segment_info_t* segment_inf,
                                   timing_info_t timing_inf);

static void compute_delta_fb_to_io(router_opts_t router_opts);

static void compute_delta_io_to_io(router_opts_t router_opts);

static void compute_delta_arrays(router_opts_t router_opts,
                                 detail_routing_arch_t det_routing_arch,
                                 segment_info_t* segment_inf,
                                 timing_info_t timing_inf,
                                 int longest_length);

static int get_first_pin(pin_types_t pintype,
                         block_type_ptr type);

static int get_longest_segment_length(detail_routing_arch_t
                                      det_routing_arch,
                                      segment_info_t* segment_inf);
static void reset_placement(void);

#ifdef PRINT_ARRAYS
static void print_array(double** array_to_print,
                        int x1,
                        int x2,
                        int y1,
                        int y2);
#endif
/**************************************/
static int
get_first_pin(pin_types_t pintype,
              block_type_ptr type)
{

    /*this code assumes logical equivilance between all driving pins */
    /*global pins are not hooked up to the temporary net */

    int i, currpin;

    currpin = 0;

    for (i = 0; i < type->num_class; i++) {
        if (type->class_inf[i].type == pintype
                && !type->is_global_pin[currpin]) {
            return (type->class_inf[i].pinlist[0]);
        } else {
            currpin += type->class_inf[i].num_pins;
        }
    }

    assert(0);
    exit(0);            /*should never hit this line */
}

/**************************************/
static int
get_longest_segment_length(detail_routing_arch_t det_routing_arch,
                           segment_info_t* segment_inf)
{

    int i, length;

    length = 0;

    for (i = 0; i < det_routing_arch.num_segment; i++) {
        if (segment_inf[i].length > length) {
            length = segment_inf[i].length;
        }
    }

    return (length);
}

/**************************************/
static void alloc_net(void)
{
    int i, len;
    net = (net_t*)my_malloc(num_nets * sizeof(net_t));

    for (i = 0; i < NET_COUNT; i++) {
        /* FIXME: We *really* shouldn't be allocating write-once copies */
        len = strlen("TEMP_NET");
        net[i].name = (char*)my_malloc((len + 1) * sizeof(char));
        net[i].is_global = FALSE;
        strcpy(net[NET_USED].name, "TEMP_NET");

        net[i].num_net_pins = (BLOCK_COUNT - 1);
        net[i].node_blocks = (int*)my_malloc(BLOCK_COUNT * sizeof(int));
        net[i].node_blocks[NET_USED_SOURCE_BLOCK] = NET_USED_SOURCE_BLOCK; /*driving blocks */
        net[i].node_blocks[NET_USED_SINK_BLOCK] = NET_USED_SINK_BLOCK;   /*target blocks */

        net[i].node_block_pins = (int*)my_malloc(BLOCK_COUNT * sizeof(int));
        /*the values for this are allocated in assign_blocks_and_route_net */
    }
}

/**************************************/
static void
alloc_block(void)
{

    /*allocates blocks structure, and assigns values to known parameters */
    /*type and x,y fields are left undefined at this stage since they  */
    /*are not known until we start moving blocks through the clb array */

    int ix_b, ix_p, len, i;
    int max_pins;

    max_pins = 0;

    for (i = 0; i < NUM_TYPES_USED; i++) {
        max_pins = max(max_pins, type_descriptors[i].num_pins);
    }

    blocks = (block_t*)my_malloc(num_blocks * sizeof(block_t));

    for (ix_b = 0; ix_b < BLOCK_COUNT; ix_b++) {
        len = strlen("TEMP_BLOCK");
        blocks[ix_b].name = (char*)my_malloc((len + 1) * sizeof(char));
        strcpy(blocks[ix_b].name, "TEMP_BLOCK");

        blocks[ix_b].nets = (int*)my_malloc(max_pins * sizeof(int));
        blocks[ix_b].nets[0] = 0;

        for (ix_p = 1; ix_p < max_pins; ix_p++) {
            blocks[ix_b].nets[ix_p] = OPEN;
        }
    }
}

/**************************************/
static void load_simplified_device(void)
{
    /* Backup original globals */
    EMPTY_TYPE_BACKUP = EMPTY_TYPE;
    IO_TYPE_BACKUP = IO_TYPE;
    FILL_TYPE_BACKUP = CLB_TYPE;
    type_descriptors_backup = type_descriptors;
    num_types_backup = num_types;
    num_types = NUM_TYPES_USED;

    /* Fill in homogeneous core type info */
    dummy_type_descriptors[0] = *EMPTY_TYPE;
    dummy_type_descriptors[0].index = 0;
    dummy_type_descriptors[1] = *IO_TYPE;
    dummy_type_descriptors[1].index = 1;
    dummy_type_descriptors[2] = *CLB_TYPE;
    dummy_type_descriptors[2].index = 2;
    type_descriptors = dummy_type_descriptors;
    EMPTY_TYPE = &dummy_type_descriptors[0];
    IO_TYPE = &dummy_type_descriptors[1];
    CLB_TYPE = &dummy_type_descriptors[2];

    /* Fill in homogeneous core bin_grids info */
    grid_backup = bin_grids;
    bin_grids = (grid_tile_t**)alloc_matrix(0, num_grid_columns + 1,
                                            0, num_grid_rows + 1,
                                            sizeof(grid_tile_t));
    int i, j;
    for (i = 0; i < num_grid_columns + 2; i++) {
        for (j = 0; j < num_grid_rows + 2; j++) {
            if ((i == 0 && j == 0) ||
                    (i == num_grid_columns + 1 && j == 0) ||
                    (i == 0 && j == num_grid_rows + 1) || (i == num_grid_columns + 1
                                                         && j == num_grid_rows + 1)) {
                bin_grids[i][j].grid_type = EMPTY_TYPE;
            } else if (i == 0 || i == num_grid_columns + 1 || j == 0 || j == num_grid_rows + 1) {
                bin_grids[i][j].grid_type = IO_TYPE;
            } else {
                bin_grids[i][j].grid_type = CLB_TYPE;
            }

            bin_grids[i][j].in_blocks =
                (int*)my_malloc(bin_grids[i][j].grid_type->capacity * sizeof(int));
            bin_grids[i][j].m_offset = 0;
        }
    }
}

static void restore_original_device(void)
{
    /* restore previous globals */
    IO_TYPE = IO_TYPE_BACKUP;
    EMPTY_TYPE = EMPTY_TYPE_BACKUP;
    CLB_TYPE = FILL_TYPE_BACKUP;
    type_descriptors = type_descriptors_backup;
    num_types = num_types_backup;

    /* free allocatd data */
    int i, j;
    for (i = 0; i < num_grid_columns + 2; ++i) {
        for (j = 0; j < num_grid_rows + 2; ++j) {
            free(bin_grids[i][j].in_blocks);
        }
    }

    free_matrix(bin_grids, 0, num_grid_columns + 1, 0, sizeof(grid_tile_t));
    bin_grids = grid_backup;
}

/**************************************/
static void reset_placement(void)
{
    int i, j, k;
    for (i = 0; i <= num_grid_columns + 1; i++) {
        for (j = 0; j <= num_grid_rows + 1; j++) {
            bin_grids[i][j].m_usage = 0;

            for (k = 0; k < bin_grids[i][j].grid_type->capacity; k++) {
                bin_grids[i][j].in_blocks[k] = EMPTY;
            }
        }
    }
}

/**************************************/
static void
alloc_and_assign_internal_structures(net_t** original_net,
                                     block_t** original_block,
                                     int* original_num_nets,
                                     int* original_num_blocks)
{
    /*allocate new data structures to hold net, and blocks info */
    *original_net = net;
    *original_num_nets = num_nets;
    num_nets = NET_COUNT;
    alloc_net();

    *original_block = blocks;
    *original_num_blocks = num_blocks;
    num_blocks = BLOCK_COUNT;
    alloc_block();


    /* [0..num_nets-1][1..num_pins-1] */
    net_delay =
        (double**)alloc_matrix(0, NET_COUNT - 1, 1, BLOCK_COUNT - 1,
                              sizeof(double));
    net_slack =
        (double**)alloc_matrix(0, NET_COUNT - 1, 1, BLOCK_COUNT - 1,
                              sizeof(double));

    reset_placement();
}

/**************************************/
static void
free_and_reset_internal_structures(net_t* original_net,
                                   block_t* original_block,
                                   int original_num_nets,
                                   int original_num_blocks)
{
    /*reset gloabal data structures to the state that they were in before these */
    /*lookup computation routines were called */

    int i;


    /*there should be only one net to free, but this is safer */
    for (i = 0; i < NET_COUNT; i++) {
        free(net[i].name);
        free(net[i].node_blocks);
        free(net[i].node_block_pins);
    }

    free(net);
    net = original_net;

    for (i = 0; i < BLOCK_COUNT; i++) {
        free(blocks[i].name);
        free(blocks[i].nets);
    }

    free(blocks);
    blocks = original_block;

    num_nets = original_num_nets;
    num_blocks = original_num_blocks;

    free_matrix(net_delay, 0, NET_COUNT - 1, 1, sizeof(double));
    free_matrix(net_slack, 0, NET_COUNT - 1, 1, sizeof(double));

}

/**************************************/
static void
setup_chan_width(router_opts_t router_opts,
                 chan_width_distr_t chan_width_dist)
{
    /*we give plenty of tracks, this increases routability for the */
    /*lookup table generation */
    int width_fac, i, max_pins_per_fb;

    max_pins_per_fb = 0;

    for (i = 0; i < num_types; i++) {
        max_pins_per_fb =
            max(max_pins_per_fb, type_descriptors[i].num_pins);
    }

    if (router_opts.fixed_channel_width == NO_FIXED_CHANNEL_WIDTH) {
        width_fac = 4 * max_pins_per_fb;    /*this is 2x the value that binary search starts */
    }
    /*this should be enough to allow most pins to   */
    /*connect to tracks in the architecture */
    else {
        width_fac = router_opts.fixed_channel_width;
    }

    init_chan(width_fac, chan_width_dist);
}

/**************************************/
static void alloc_routing_structs(router_opts_t router_opts,
                                  detail_routing_arch_t det_routing_arch,
                                  segment_info_t* segment_inf,
                                  timing_info_t timing_inf,
                                  subblock_data_t subblock_data)
{

    int bb_factor;
    int warnings;
    t_graph_type graph_type;

    /*calls routines that set up routing resource graph and associated structures */


    /*must set up dummy blocks for the first pass through to setup locally used opins */
    /* Only one blocks per tile */
    assign_locations(CLB_TYPE, 1, 1, 0,
                     CLB_TYPE, num_grid_columns, num_grid_rows, 0);

    clb_opins_used_locally = alloc_route_structs(subblock_data);

    free_rr_graph();

    if (router_opts.route_type == GLOBAL) {
        graph_type = GRAPH_GLOBAL;
    } else {
        graph_type = (det_routing_arch.directionality ==
                      BI_DIRECTIONAL ? GRAPH_BIDIR : GRAPH_UNIDIR);
    }

    build_rr_graph(graph_type, num_types, dummy_type_descriptors, num_grid_columns,
                   //             num_grid_rows, bin_grids, chan_width_x[0], NULL,
                   num_grid_rows, bin_grids,    16          , NULL,
                   det_routing_arch.switch_block_type,
                   det_routing_arch.Fs, det_routing_arch.num_segment, det_routing_arch.num_switch,
                   segment_inf, det_routing_arch.global_route_switch,
                   det_routing_arch.delayless_switch, timing_inf,
                   det_routing_arch.wire_to_ipin_switch,
                   router_opts.base_cost_type, &warnings);

    alloc_and_load_rr_node_route_structs();

    alloc_timing_driven_route_structs(&pin_criticality, &sink_order,
                                      &rt_node_of_sink);


    bb_factor = num_grid_columns + num_grid_rows;   /*set it to a huge value */
    init_route_structs(bb_factor);
}

/**************************************/
static void
free_routing_structs(router_opts_t router_opts,
                     detail_routing_arch_t det_routing_arch,
                     segment_info_t* segment_inf,
                     timing_info_t timing_inf)
{
    free_rr_graph();

    free_rr_node_route_structs();
    free_route_structs(clb_opins_used_locally);
    free_trace_structs();

    free_timing_driven_route_structs(pin_criticality, sink_order,
                                     rt_node_of_sink);
}

/**************************************/
static void assign_locations(block_type_ptr source_type,
                             int source_x_loc,
                             int source_y_loc,
                             int source_z_loc,
                             block_type_ptr sink_type,
                             int sink_x_loc,
                             int sink_y_loc,
                             int sink_z_loc)
{
    /*all routing occurs between blocks 0 (source) and blocks 1 (sink) */
    blocks[SOURCE_BLOCK].block_type = source_type;
    blocks[SOURCE_BLOCK].x = source_x_loc;
    blocks[SOURCE_BLOCK].y = source_y_loc;
    blocks[SOURCE_BLOCK].z = source_z_loc;

    blocks[SINK_BLOCK].block_type = sink_type;
    blocks[SINK_BLOCK].x = sink_x_loc;
    blocks[SINK_BLOCK].y = sink_y_loc;
    blocks[SINK_BLOCK].z = sink_z_loc;

    bin_grids[source_x_loc][source_y_loc].in_blocks[source_z_loc] = SOURCE_BLOCK;
    bin_grids[sink_x_loc][sink_y_loc].in_blocks[sink_z_loc] = SINK_BLOCK;

    net[NET_USED].node_block_pins[NET_USED_SOURCE_BLOCK] =
        get_first_pin(DRIVER, blocks[SOURCE_BLOCK].block_type);
    net[NET_USED].node_block_pins[NET_USED_SINK_BLOCK] =
        get_first_pin(RECEIVER, blocks[SINK_BLOCK].block_type);

    bin_grids[source_x_loc][source_y_loc].m_usage += 1;
    bin_grids[sink_x_loc][sink_y_loc].m_usage += 1;

}

/**************************************/
static double assign_blocks_and_route_net(block_type_ptr source_type,
                                          int source_x_loc,
                                          int source_y_loc,
                                          block_type_ptr sink_type,
                                          int sink_x_loc,
                                          int sink_y_loc,
                                          router_opts_t router_opts)
{
    /*places blocks at the specified locations, and routes a net between them */
    /*returns the Tdel of this net */
    /* Only one blocks per tile */
    int source_z_loc = 0;
    int sink_z_loc = 0;

    assign_locations(source_type,
                     source_x_loc,
                     source_y_loc,
                     source_z_loc,
                     sink_type,
                     sink_x_loc,
                     sink_y_loc,
                     sink_z_loc);

    load_net_rr_terminals(rr_node_indices);

    double T_crit = 1.0;
    double pres_fac = 0;  /* ignore congestion */

    const int knum_net_pins = net[NET_USED].num_net_pins;
    int ipin;
    for (ipin = 1; ipin <= knum_net_pins; ++ipin) {
        net_slack[NET_USED][ipin] = 0;
    }

    boolean is_routeable = timing_driven_route_net(NET_USED, pres_fac,
                                           router_opts.max_criticality,
                                           router_opts.criticality_exp,
                                           router_opts.astar_fac,
                                           router_opts.bend_cost,
                                           net_slack[NET_USED],
                                           pin_criticality, sink_order,
                                           rt_node_of_sink, T_crit,
                                           net_delay[NET_USED]);

    double net_delay_value = net_delay[NET_USED][NET_USED_SINK_BLOCK];

    bin_grids[source_x_loc][source_y_loc].m_usage = 0;
    bin_grids[source_x_loc][source_y_loc].in_blocks[source_z_loc] = EMPTY;
    bin_grids[sink_x_loc][sink_y_loc].m_usage = 0;
    bin_grids[sink_x_loc][sink_y_loc].in_blocks[sink_z_loc] = EMPTY;

    return (net_delay_value);
}

/**************************************/
static void alloc_delta_arrays(void)
{
    delta_clb_to_clb = (double**)alloc_matrix(0, num_grid_columns - 1,
                                              0, num_grid_rows - 1,
                                              sizeof(double));
    delta_inpad_to_clb = (double**)alloc_matrix(0, num_grid_columns,
                                                0, num_grid_rows,
                                                sizeof(double));
    delta_clb_to_outpad = (double**)alloc_matrix(0, num_grid_columns,
                                                 0, num_grid_rows,
                                                 sizeof(double));
    delta_inpad_to_outpad = (double**)alloc_matrix(0, num_grid_columns + 1,
                                                   0, num_grid_rows + 1,
                                                   sizeof(double));


    /*initialize all of the array locations to -1 */
    int id_x, id_y;
    for (id_x = 0; id_x <= num_grid_columns; id_x++) {
        for (id_y = 0; id_y <= num_grid_rows; id_y++) {
            delta_inpad_to_clb[id_x][id_y] = IMPOSSIBLE;
        }
    }

    for (id_x = 0; id_x <= num_grid_columns - 1; id_x++) {
        for (id_y = 0; id_y <= num_grid_rows - 1; id_y++) {
            delta_clb_to_clb[id_x][id_y] = IMPOSSIBLE;
        }
    }

    for (id_x = 0; id_x <= num_grid_columns; id_x++) {
        for (id_y = 0; id_y <= num_grid_rows; id_y++) {
            delta_clb_to_outpad[id_x][id_y] = IMPOSSIBLE;
        }
    }

    for (id_x = 0; id_x <= num_grid_columns + 1; id_x++) {
        for (id_y = 0; id_y <= num_grid_rows + 1; id_y++) {
            delta_inpad_to_outpad[id_x][id_y] = IMPOSSIBLE;
        }
    }
}

/**************************************/
static void free_delta_arrays(void)
{
    free_matrix(delta_inpad_to_clb, 0, num_grid_columns, 0, sizeof(double));
    free_matrix(delta_clb_to_clb, 0, num_grid_columns - 1, 0, sizeof(double));
    free_matrix(delta_clb_to_outpad, 0, num_grid_columns, 0, sizeof(double));
    free_matrix(delta_inpad_to_outpad, 0, num_grid_columns + 1, 0, sizeof(double));
}

/**************************************/
static void
generic_compute_matrix(double** *matrix_ptr,
                       block_type_ptr source_type,
                       block_type_ptr sink_type,
                       int source_x,
                       int source_y,
                       int start_x,
                       int end_x,
                       int start_y,
                       int end_y,
                       router_opts_t router_opts,
                       detail_routing_arch_t det_routing_arch,
                       segment_info_t* segment_inf,
                       timing_info_t timing_inf)
{

    int delta_x, delta_y;
    int sink_x, sink_y;

    for (sink_x = start_x; sink_x <= end_x; sink_x++) {
        for (sink_y = start_y; sink_y <= end_y; sink_y++) {
            delta_x = abs(sink_x - source_x);
            delta_y = abs(sink_y - source_y);

            if (delta_x == 0 && delta_y == 0) {
                continue;    /*do not compute distance from a blocks to itself     */
            }

            /*if a value is desired, pre-assign it somewhere else */

            (*matrix_ptr)[delta_x][delta_y] =
                assign_blocks_and_route_net(source_type, source_x,
                                            source_y, sink_type,
                                            sink_x, sink_y,
                                            router_opts);
        }
    }
}

/**************************************/
/*this routine must compute Tdel values in a slightly different way than the */
/*other compute routines. We cannot use a location close to the tedge as the  */
/*source location for the majority of the Tdel computations  because this   */
/*would give gradually increasing Tdel values. To avoid this from happening */
/*a clb that is at least longest_length away from an tedge should be chosen   */
/*as a source , if longest_length is more than 0.5 of the total size then    */
/*choose a FB at the center as the source FB */
static void compute_delta_fb_to_fb(router_opts_t router_opts,
                                   int longest_length)
{
    block_type_ptr source_type = CLB_TYPE;
    block_type_ptr sink_type = CLB_TYPE;

    int start_x;
    if (longest_length < 0.5 * (num_grid_columns)) {
        start_x = longest_length;
    } else {
        start_x = (int)(0.5 * num_grid_columns);
    }

    int end_x = num_grid_columns;
    int source_x = start_x;

    int start_y = 0;
    if (longest_length < 0.5 * (num_grid_rows)) {
        start_y = longest_length;
    } else {
        start_y = (int)(0.5 * num_grid_rows);
    }

    int end_y = num_grid_rows;
    int source_y = start_y;

    /*don't put the sink all the way to the corner, until it is necessary */
    int sink_x, sink_y;
    int delta_x, delta_y;
    for (sink_x = start_x; sink_x <= end_x - 1; sink_x++) {
        for (sink_y = start_y; sink_y <= end_y - 1; sink_y++) {
            delta_x = abs(sink_x - source_x);
            delta_y = abs(sink_y - source_y);

            if (delta_x == 0 && delta_y == 0) {
                delta_clb_to_clb[delta_x][delta_y] = 0.0;
                continue;
            }

            delta_clb_to_clb[delta_x][delta_y] =
                assign_blocks_and_route_net(source_type, source_x,
                                            source_y, sink_type,
                                            sink_x, sink_y,
                                            router_opts);
        }

    }


    sink_x = end_x - 1;
    sink_y = end_y - 1;

    for (source_x = start_x - 1; source_x >= 1; source_x--) {
        for (source_y = start_y; source_y <= end_y - 1; source_y++) {
            delta_x = abs(sink_x - source_x);
            delta_y = abs(sink_y - source_y);

            delta_clb_to_clb[delta_x][delta_y] =
                assign_blocks_and_route_net(source_type, source_x,
                                            source_y, sink_type,
                                            sink_x, sink_y,
                                            router_opts);
        }
    }

    for (source_x = 1; source_x <= end_x - 1; source_x++) {
        for (source_y = 1; source_y < start_y; source_y++) {
            delta_x = abs(sink_x - source_x);
            delta_y = abs(sink_y - source_y);

            delta_clb_to_clb[delta_x][delta_y] =
                assign_blocks_and_route_net(source_type, source_x,
                                            source_y, sink_type,
                                            sink_x, sink_y,
                                            router_opts);
        }
    }


    /*now move sink into the top right corner */
    sink_x = end_x;
    sink_y = end_y;
    source_x = 1;

    for (source_y = 1; source_y <= end_y; source_y++) {
        delta_x = abs(sink_x - source_x);
        delta_y = abs(sink_y - source_y);

        delta_clb_to_clb[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);

    }

    sink_x = end_x;
    sink_y = end_y;
    source_y = 1;

    for (source_x = 1; source_x <= end_x; source_x++) {
        delta_x = abs(sink_x - source_x);
        delta_y = abs(sink_y - source_y);

        delta_clb_to_clb[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);
    }
}

/**************************************/
static void
compute_delta_io_to_fb(router_opts_t router_opts,
                       detail_routing_arch_t det_routing_arch,
                       segment_info_t* segment_inf,
                       timing_info_t timing_inf)
{
    int source_x, source_y;
    int start_x, start_y, end_x, end_y;
    block_type_ptr source_type, sink_type;

    source_type = IO_TYPE;
    sink_type = CLB_TYPE;

    delta_inpad_to_clb[0][0] = IMPOSSIBLE;
    delta_inpad_to_clb[num_grid_columns][num_grid_rows] = IMPOSSIBLE;

    source_x = 0;
    source_y = 1;

    start_x = 1;
    end_x = num_grid_columns;
    start_y = 1;
    end_y = num_grid_rows;
    generic_compute_matrix(&delta_inpad_to_clb, source_type, sink_type,
                           source_x, source_y, start_x, end_x, start_y,
                           end_y, router_opts, det_routing_arch,
                           segment_inf, timing_inf);

    source_x = 1;
    source_y = 0;

    start_x = 1;
    end_x = 1;
    start_y = 1;
    end_y = num_grid_rows;
    generic_compute_matrix(&delta_inpad_to_clb, source_type, sink_type,
                           source_x, source_y, start_x, end_x, start_y,
                           end_y, router_opts, det_routing_arch,
                           segment_inf, timing_inf);

    start_x = 1;
    end_x = num_grid_columns;
    start_y = num_grid_rows;
    end_y = num_grid_rows;
    generic_compute_matrix(&delta_inpad_to_clb, source_type, sink_type,
                           source_x, source_y, start_x, end_x, start_y,
                           end_y, router_opts, det_routing_arch,
                           segment_inf, timing_inf);
}

/**************************************/
static void compute_delta_fb_to_io(router_opts_t router_opts)
{
    block_type_ptr source_type = CLB_TYPE;
    block_type_ptr sink_type = IO_TYPE;

    delta_clb_to_outpad[0][0] = IMPOSSIBLE;
    delta_clb_to_outpad[num_grid_columns][num_grid_rows] = IMPOSSIBLE;

    int sink_x = 0;
    int sink_y = 1;

    int delta_x, delta_y;
    int source_x, source_y;
    for (source_x = 1; source_x <= num_grid_columns; source_x++) {
        for (source_y = 1; source_y <= num_grid_rows; source_y++) {
            delta_x = abs(source_x - sink_x);
            delta_y = abs(source_y - sink_y);

            delta_clb_to_outpad[delta_x][delta_y] =
                assign_blocks_and_route_net(source_type, source_x,
                                            source_y, sink_type,
                                            sink_x, sink_y,
                                            router_opts);
        }
    }

    sink_x = 1;
    sink_y = 0;
    source_x = 1;
    delta_x = abs(source_x - sink_x);

    for (source_y = 1; source_y <= num_grid_rows; source_y++) {
        delta_y = abs(source_y - sink_y);
        delta_clb_to_outpad[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);
    }

    sink_x = 1;
    sink_y = 0;
    source_y = num_grid_rows;
    delta_y = abs(source_y - sink_y);

    for (source_x = 2; source_x <= num_grid_columns; source_x++) {
        delta_x = abs(source_x - sink_x);
        delta_clb_to_outpad[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);
    }
}

/**************************************/
static void compute_delta_io_to_io(router_opts_t router_opts)
{
    block_type_ptr source_type = IO_TYPE;
    block_type_ptr sink_type = IO_TYPE;

    delta_inpad_to_outpad[0][0] = 0;   /*Tdel to itself is 0 (this can happen) */
    delta_inpad_to_outpad[num_grid_columns + 1][num_grid_rows + 1] = IMPOSSIBLE;
    delta_inpad_to_outpad[0][num_grid_rows] = IMPOSSIBLE;
    delta_inpad_to_outpad[num_grid_columns][0] = IMPOSSIBLE;
    delta_inpad_to_outpad[num_grid_columns][num_grid_rows + 1] = IMPOSSIBLE;
    delta_inpad_to_outpad[num_grid_columns + 1][num_grid_rows] = IMPOSSIBLE;


    int source_x = 0;
    int source_y = 1;
    int sink_x = 0;
    int delta_x = abs(sink_x - source_x);

    int sink_y;
    int delta_y;
    for (sink_y = 2; sink_y <= num_grid_rows; sink_y++) {
        delta_y = abs(sink_y - source_y);
        delta_inpad_to_outpad[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);

    }

    source_x = 0;
    source_y = 1;
    sink_x = num_grid_columns + 1;
    delta_x = abs(sink_x - source_x);
    for (sink_y = 1; sink_y <= num_grid_rows; sink_y++) {
        delta_y = abs(sink_y - source_y);
        delta_inpad_to_outpad[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);

    }


    source_x = 1;
    source_y = 0;
    sink_y = 0;
    delta_y = abs(sink_y - source_y);

    for (sink_x = 2; sink_x <= num_grid_columns; sink_x++) {
        delta_x = abs(sink_x - source_x);
        delta_inpad_to_outpad[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);

    }

    source_x = 1;
    source_y = 0;
    sink_y = num_grid_rows + 1;
    delta_y = abs(sink_y - source_y);

    for (sink_x = 1; sink_x <= num_grid_columns; sink_x++) {
        delta_x = abs(sink_x - source_x);
        delta_inpad_to_outpad[delta_x][delta_y] =
            assign_blocks_and_route_net(source_type, source_x, source_y,
                                        sink_type, sink_x, sink_y,
                                        router_opts);

    }

    source_x = 0;
    sink_y = num_grid_rows + 1;

    for (source_y = 1; source_y <= num_grid_rows; source_y++) {
        for (sink_x = 1; sink_x <= num_grid_columns; sink_x++) {
            delta_y = abs(source_y - sink_y);
            delta_x = abs(source_x - sink_x);
            delta_inpad_to_outpad[delta_x][delta_y] =
                assign_blocks_and_route_net(source_type, source_x,
                                            source_y, sink_type,
                                            sink_x, sink_y,
                                            router_opts);

        }
    }
}

/**************************************/
#ifdef PRINT_ARRAYS
static void
print_array(double** array_to_print,
            int x1,
            int x2,
            int y1,
            int y2)
{


    int idx_x, idx_y;

    fprintf(lookup_dump, "\nPrinting Array \n\n");

    for (idx_y = y2; idx_y >= y1; idx_y--) {
        for (idx_x = x1; idx_x <= x2; idx_x++) {
            fprintf(lookup_dump, " %9.2e",
                    array_to_print[idx_x][idx_y]);
        }

        fprintf(lookup_dump, "\n");
    }

    fprintf(lookup_dump, "\n\n");
}
#endif
/**************************************/
static void
compute_delta_arrays(router_opts_t router_opts,
                     detail_routing_arch_t det_routing_arch,
                     segment_info_t* segment_inf,
                     timing_info_t timing_inf,
                     int longest_length)
{

    printf
    ("Computing delta_inpad_to_outpad lookup matrix, may take a few seconds, please wait...\n");
    compute_delta_io_to_io(router_opts);
    printf
    ("Computing delta_inpad_to_clb> lookup matrix, may take a few seconds, please wait...\n");
    compute_delta_io_to_fb(router_opts,
                           det_routing_arch,
                           segment_inf,
                           timing_inf);
    printf
    ("Computing delta_clb_to_outpad> lookup matrix, may take a few seconds, please wait...\n");
    compute_delta_fb_to_io(router_opts);
    printf
    ("Computing delta_clb_to_clb lookup matrix, may take a few seconds, please wait...\n");
    compute_delta_fb_to_fb(router_opts,
                           longest_length);

#ifdef PRINT_ARRAYS
    lookup_dump = my_fopen(DUMPFILE, "w");
    fprintf(lookup_dump, "\n\nprinting delta_clb_to_clb\n");
    print_array(delta_clb_to_clb, 0, num_grid_columns - 1, 0, num_grid_rows - 1);
    fprintf(lookup_dump, "\n\nprinting delta_inpad_to_clb>\n");
    print_array(delta_inpad_to_clb>, 0, num_grid_columns, 0, num_grid_rows);
    fprintf(lookup_dump, "\n\nprinting delta_clb_to_outpad>\n");
    print_array(delta_clb_to_outpad>, 0, num_grid_columns, 0, num_grid_rows);
    fprintf(lookup_dump, "\n\nprinting delta_inpad_to_outpad\n");
    print_array(delta_inpad_to_outpad, 0, num_grid_columns + 1, 0, num_grid_rows + 1);
    fclose(lookup_dump);
#endif

}

/******* Globally Accessable Functions **********/

/**************************************/
void compute_delay_lookup_tables(router_opts_t router_opts,
                                 detail_routing_arch_t det_routing_arch,
                                 segment_info_t* segment_inf,
                                 timing_info_t timing_inf,
                                 chan_width_distr_t chan_width_dist,
                                 subblock_data_t subblock_data)
{
    static net_t* original_net;  /*this will be used as a pointer to remember what */

    /*the "real" nets in the circuit are. This is    */
    /*required because we are using the net structure */
    /*in these routines to find delays between blocks */
    static block_t* original_block;  /*same def as original_nets, but for blocks  */

    static int original_num_nets;
    static int original_num_blocks;
    static int longest_length;

    load_simplified_device();

    alloc_and_assign_internal_structures(&original_net,
                                         &original_block,
                                         &original_num_nets,
                                         &original_num_blocks);
    setup_chan_width(router_opts, chan_width_dist);

    alloc_routing_structs(router_opts, det_routing_arch, segment_inf,
                          timing_inf, subblock_data);

    longest_length =
        get_longest_segment_length(det_routing_arch, segment_inf);


    /*now setup and compute the actual arrays */
    alloc_delta_arrays();
    compute_delta_arrays(router_opts, det_routing_arch, segment_inf,
                         timing_inf, longest_length);

    /*free all data structures that are no longer needed */
    free_routing_structs(router_opts, det_routing_arch, segment_inf,
                         timing_inf);

    restore_original_device();

    free_and_reset_internal_structures(original_net, original_block,
                                       original_num_nets,
                                       original_num_blocks);
}

/**************************************/
void
free_place_lookup_structs(void)
{

    free_delta_arrays();

}
