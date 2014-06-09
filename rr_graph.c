#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "rr_graph_util.h"
#include "rr_graph.h"
#include "rr_graph2.h"
#include "rr_graph_sbox.h"
#include "check_rr_graph.h"
#include "rr_graph_timing_params.h"
#include "rr_graph_indexed_data.h"
#include "vpr_utils.h"

/* UDSD Modifications by WMF Begin */
/* WMF: I'm using this macro to reverse the order for wire muxes, so opins
 * connect to the highest track first. This way, the leftover imbalance
 * is in reverse to the leftover imbalance due to Fs generation
 * May 15, I revamped the imbalance solution, so this feature is not useful
 * anymore. */
#define REVERSE_OPIN_ORDER

/* #define ENABLE_DUMP */
#define MUX_SIZE_DIST_DISPLAY

/* mux size statistic data structures */

typedef struct s_mux {
    int size;
    struct s_mux* next;
}
t_mux;

typedef struct s_mux_size_distribution {
    int mux_count;
    int max_index;
    int* distr;
    struct s_mux_size_distribution* next;
}
t_mux_size_distribution;

/* UDSD Modifications by WMF End */


/******************* Variables local to this module. ***********************/

static linked_vptr_t* rr_mem_chunk_list_head = NULL;

/* Used to free "chunked" memory.  If NULL, no rr_graph exists right now.  */

static int chunk_bytes_avail = 0;
static char* chunk_next_avail_mem = NULL;

/* Status of current chunk being dished out by calls to my_chunk_malloc.   */

/********************* Subroutines local to this module. *******************/
static int**** alloc_and_load_pin_to_track_map(IN pin_types_t pin_type,
                                               IN int nodes_per_chan,
                                               IN int Fc,
                                               IN block_type_ptr Type,
                                               IN boolean
                                               perturb_switch_pattern,
                                               IN directionality_t
                                               directionality);

static vector_t*** alloc_and_load_track_to_pin_lookup(IN int
                                                           **** pin_to_track_map,
                                                           IN int Fc,
                                                           IN int height,
                                                           IN int num_pins,
                                                           IN int
                                                           nodes_per_chan);

static void build_bidir_rr_opins(IN int i,
                                 IN int j,
                                 INOUT rr_node_t* rr_node,
                                 IN vector_t** * rr_node_indices,
                                 IN int**** *opin_to_track_map,
                                 IN int* Fc_out,
                                 IN boolean* rr_edge_done,
                                 IN segment_details_t* seg_details,
                                 IN grid_tile_t** grid);

static void build_unidir_rr_opins(IN int i,
                                  IN int j,
                                  IN grid_tile_t** grid,
                                  IN int* Fc_out,
                                  IN int nodes_per_chan,
                                  IN segment_details_t* seg_details,
                                  INOUT int** Fc_xofs,
                                  INOUT int** Fc_yofs,
                                  INOUT rr_node_t* rr_node,
                                  INOUT boolean* rr_edge_done,
                                  OUT boolean* Fc_clipped,
                                  IN vector_t** * rr_node_indices);

static void alloc_and_load_rr_graph(IN int num_nodes,
                                    IN rr_node_t* rr_node,
                                    IN int num_seg_types,
                                    IN segment_details_t* seg_details,
                                    IN boolean* rr_edge_done,
                                    IN vector_t**** track_to_ipin_lookup,
                                    IN int**** *opin_to_track_map,
                                    IN vector_t** *switch_block_conn,
                                    IN grid_tile_t** grid,
                                    IN int num_grid_columns,
                                    IN int num_grid_rows,
                                    IN int Fs,
                                    IN short**** *sblock_pattern,
                                    IN int* Fc_out,
                                    IN int** Fc_xofs,
                                    IN int** Fc_yofs,
                                    IN vector_t** * rr_node_indices,
                                    IN int nodes_per_chan,
                                    IN switch_block_t sb_type,
                                    IN int delayless_switch,
                                    IN directionality_t directionality,
                                    IN int wire_to_ipin_switch,
                                    OUT boolean* Fc_clipped);

static void load_uniform_switch_pattern(IN block_type_ptr type,
                                        INOUT int**** tracks_connected_to_pin,
                                        IN int num_phys_pins,
                                        IN int* pin_num_ordering,
                                        IN int* side_ordering,
                                        IN int* offset_ordering,
                                        IN int nodes_per_chan,
                                        IN int Fc,
                                        IN directionality_t
                                        directionality);

static void load_perturbed_switch_pattern(IN block_type_ptr type,
                                          INOUT int
                                          **** tracks_connected_to_pin,
                                          IN int num_phys_pins,
                                          IN int* pin_num_ordering,
                                          IN int* side_ordering,
                                          IN int* offset_ordering,
                                          IN int nodes_per_chan,
                                          IN int Fc,
                                          IN directionality_t
                                          directionality);

static void check_all_tracks_reach_pins(block_type_ptr type,
                                        int**** tracks_connected_to_pin,
                                        int nodes_per_chan,
                                        int Fc,
                                        pin_types_t ipin_or_opin);

static boolean* alloc_and_load_perturb_ipins(IN int nodes_per_chan,
                                             IN int num_types,
                                             IN int* Fc_in,
                                             IN int* Fc_out,
                                             IN directionality_t
                                             directionality);


static void print_rr_node(FILE* fp,
                          rr_node_t* rr_node,
                          int inode);

static void build_rr_sinks_sources(IN int i,
                                   IN int j,
                                   IN rr_node_t* rr_node,
                                   IN vector_t** * rr_node_indices,
                                   IN int delayless_switch,
                                   IN grid_tile_t** grid);

static void build_rr_xchan(IN int i,
                           IN int j,
                           IN vector_t**** track_to_ipin_lookup,
                           IN vector_t** *switch_block_conn,
                           IN int cost_index_offset,
                           IN int nodes_per_chan,
                           IN int* opin_mux_size,
                           IN short**** *sblock_pattern,
                           IN int Fs_per_side,
                           IN segment_details_t* seg_details,
                           IN vector_t** * rr_node_indices,
                           IN boolean* rr_edge_done,
                           INOUT rr_node_t* rr_node,
                           IN int wire_to_ipin_switch,
                           IN directionality_t directionality);

static void build_rr_ychan(IN int i,
                           IN int j,
                           IN vector_t**** track_to_ipin_lookup,
                           IN vector_t** *switch_block_conn,
                           IN int cost_index_offset,
                           IN int nodes_per_chan,
                           IN int* opin_mux_size,
                           IN short**** *sblock_pattern,
                           IN int Fs_per_side,
                           IN segment_details_t* seg_details,
                           IN vector_t** * rr_node_indices,
                           IN boolean* rr_edge_done,
                           INOUT rr_node_t* rr_node,
                           IN int wire_to_ipin_switch,
                           IN directionality_t directionality);


static void rr_graph_externals(
    timing_info_t timing_inf,
    segment_info_t* segment_inf,
    int num_seg_types,
    int nodes_per_chan,
    int wire_to_ipin_switch,
    router_base_cost_t base_cost_type);

void alloc_and_load_edges_and_switches(IN rr_node_t* rr_node,
                                       IN int inode,
                                       IN int num_edges,
                                       IN boolean* rr_edge_done,
                                       IN t_linked_edge* edge_list_head);

static void alloc_net_rr_terminals(void);

static void alloc_and_load_rr_clb_source(vector_t** * rr_node_indices);

void watch_edges(int inode,
                 t_linked_edge* edge_list_head);

static void free_type_pin_to_track_map(int**** * ipin_to_track_map, block_type_ptr types, int* fc_in);

static void free_type_track_to_ipin_map(vector_t**** track_to_pin_map, block_type_ptr types, int nodes_per_chan);

static segment_details_t* alloc_and_load_global_route_seg_details(IN int
                                                              nodes_per_chan,
                                                              IN int
                                                              global_route_switch);

/* UDSD Modifications by WMF End */

static int* alloc_and_load_actual_fc(IN int num_types,
                                     IN block_type_ptr types,
                                     IN int nodes_per_chan,
                                     IN boolean is_Fc_out,
                                     IN directionality_t directionality,
                                     OUT boolean* Fc_clipped);

/******************* Subroutine definitions *******************************/


void
build_rr_graph(IN t_graph_type graph_type,
               IN int num_types,
               IN block_type_ptr types,
               IN int num_grid_columns,
               IN int num_grid_rows,
               IN grid_tile_t** grid,
               IN int chan_width,
               IN struct s_chan_width_dist* chan_capacity_inf,
               IN switch_block_t sb_type,
               IN int Fs,
               IN int num_seg_types,
               IN int num_switches,
               IN segment_info_t* segment_inf,
               IN int global_route_switch,
               IN int delayless_switch,
               IN timing_info_t timing_inf,
               IN int wire_to_ipin_switch,
               IN router_base_cost_t base_cost_type,
               OUT int* Warnings)
{
    /* Temp structures used to build graph */
    int nodes_per_chan, i, j;
    segment_details_t* seg_details = NULL;
    int* Fc_in = NULL;
    int* Fc_out = NULL;
    int**** *opin_to_track_map = NULL; /* [0..num_types-1][0..num_pins-1][0..height][0..3][0..Fc-1] */
    int**** *ipin_to_track_map = NULL; /* [0..num_types-1][0..num_pins-1][0..height][0..3][0..Fc-1] */
    vector_t**** track_to_ipin_lookup = NULL; /* [0..num_types-1][0..nodes_per_chan-1][0..height][0..3] */
    vector_t** *switch_block_conn = NULL;
    short**** *unidir_sb_pattern = NULL;
    boolean* rr_edge_done = NULL;
    boolean is_global_graph;
    boolean Fc_clipped;
    boolean use_full_seg_groups;
    boolean* perturb_ipins = NULL;
    directionality_t directionality;
    int** Fc_xofs = NULL;   /* [0..num_grid_rows-1][0..num_grid_columns-1] */
    int** Fc_yofs = NULL;   /* [0..num_grid_columns-1][0..num_grid_rows-1] */
    rr_node_indices = NULL;
    rr_node = NULL;
    num_rr_nodes = 0;
    /* Reset warning flag */
    *Warnings = RR_GRAPH_NO_WARN;
    /* Decode the graph_type */
    is_global_graph = FALSE;

    if (GRAPH_GLOBAL == graph_type) {
        is_global_graph = TRUE;
    }

    use_full_seg_groups = FALSE;

    if (GRAPH_UNIDIR_TILEABLE == graph_type) {
        use_full_seg_groups = TRUE;
    }

    directionality = UNI_DIRECTIONAL;

    if (GRAPH_BIDIR == graph_type) {
        directionality = BI_DIRECTIONAL;
    }

    if (is_global_graph) {
        directionality = BI_DIRECTIONAL;
    }

    /* Global routing uses a single longwire track */
    nodes_per_chan = (is_global_graph ? 1 : chan_width);
    assert(nodes_per_chan > 0);

    /* START SEG_DETAILS */
    if (is_global_graph) {
        /* Sets up a single unit length segment type for global routing. */
        seg_details =
            alloc_and_load_global_route_seg_details(nodes_per_chan,
                                                    global_route_switch);
    } else {
        /* Setup segments including distrubuting tracks and staggering.
         * If use_full_seg_groups is specified, nodes_per_chan may be
         * changed. Warning should be singled to caller if this happens. */
        seg_details = alloc_and_load_seg_details(&nodes_per_chan,
                                                 max(num_grid_columns, num_grid_rows),
                                                 num_seg_types,
                                                 segment_inf,
                                                 use_full_seg_groups,
                                                 is_global_graph,
                                                 directionality);

        if ((is_global_graph ? 1 : chan_width) != nodes_per_chan) {
            *Warnings |= RR_GRAPH_WARN_CHAN_WIDTH_CHANGED;
        }

#ifdef CREATE_ECHO_FILES
        dump_seg_details(seg_details, nodes_per_chan, "seg_details.txt");
#endif /* CREATE_ECHO_FILES */
    }

    /* END SEG_DETAILS */

    /* START FC */
    /* Determine the actual value of Fc */
    if (is_global_graph) {
        Fc_in = (int*)my_malloc(sizeof(int) * num_types);
        Fc_out = (int*)my_malloc(sizeof(int) * num_types);

        for (i = 0; i < num_types; ++i) {
            Fc_in[i] = 1;
            Fc_out[i] = 1;
        }
    } else {
        Fc_clipped = FALSE;
        Fc_in =
            alloc_and_load_actual_fc(num_types, types, nodes_per_chan,
                                     FALSE, directionality, &Fc_clipped);

        if (Fc_clipped) {
            *Warnings |= RR_GRAPH_WARN_FC_CLIPPED;
        }

        Fc_clipped = FALSE;
        Fc_out =
            alloc_and_load_actual_fc(num_types, types, nodes_per_chan,
                                     TRUE, directionality, &Fc_clipped);

        if (Fc_clipped) {
            *Warnings |= RR_GRAPH_WARN_FC_CLIPPED;
        }

#ifdef VERBOSE

        for (i = 1; i < num_types; ++i) {
            /* Skip "<EMPTY>" */
            if (type_descriptors[i].is_Fc_out_full_flex) {
                printf
                ("Fc Actual Values: Type = %s, Fc_out = full, Fc_in = %d.\n",
                 type_descriptors[i].name, Fc_input[i]);
            } else {
                printf
                ("Fc Actual Values: Type = %s, Fc_out = %d, Fc_in = %d.\n",
                 type_descriptors[i].name, Fc_output[i],
                 Fc_input[i]);
            }
        }

#endif /* VERBOSE */
    }

    perturb_ipins =
        alloc_and_load_perturb_ipins(nodes_per_chan, num_types, Fc_in,
                                     Fc_out, directionality);
    /* END FC */
    /* Alloc node lookups, count nodes, alloc rr nodes */
    num_rr_nodes = 0;
    rr_node_indices =
        alloc_and_load_rr_node_indices(nodes_per_chan, num_grid_columns, num_grid_rows, &num_rr_nodes,
                                       seg_details);
    rr_node = (rr_node_t*) my_malloc(sizeof(rr_node_t) * num_rr_nodes);
    memset(rr_node, 0, sizeof(rr_node_t) * num_rr_nodes);
    rr_edge_done = (boolean*) my_malloc(sizeof(boolean) * num_rr_nodes);
    memset(rr_edge_done, 0, sizeof(boolean) * num_rr_nodes);

    /* These are data structures used by the the unidir opin mapping. */
    if (UNI_DIRECTIONAL == directionality) {
        Fc_xofs = (int**)alloc_matrix(0, num_grid_rows, 0, num_grid_columns, sizeof(int));
        Fc_yofs = (int**)alloc_matrix(0, num_grid_columns, 0, num_grid_rows, sizeof(int));

        for (i = 0; i <= num_grid_columns; ++i) {
            for (j = 0; j <= num_grid_rows; ++j) {
                Fc_xofs[j][i] = 0;
                Fc_yofs[i][j] = 0;
            }
        }
    }

    /* START SB LOOKUP */
    /* Alloc and load the switch block lookup */
    if (is_global_graph) {
        assert(nodes_per_chan == 1);
        switch_block_conn =
            alloc_and_load_switch_block_conn(1, SUBSET, 3);
    } else if (BI_DIRECTIONAL == directionality) {
        switch_block_conn =
            alloc_and_load_switch_block_conn(nodes_per_chan, sb_type, Fs);
    } else {
        assert(UNI_DIRECTIONAL == directionality);
        unidir_sb_pattern =
            alloc_sblock_pattern_lookup(num_grid_columns, num_grid_rows, nodes_per_chan);

        for (i = 0; i <= num_grid_columns; i++) {
            for (j = 0; j <= num_grid_rows; j++) {
                load_sblock_pattern_lookup(i, j, nodes_per_chan,
                                           seg_details, Fs,
                                           sb_type,
                                           unidir_sb_pattern);
            }
        }
    }

    /* END SB LOOKUP */
    /* START IPIN MAP */
    /* Create ipin map lookups */
    ipin_to_track_map = (int*****)my_malloc(sizeof(int****) * num_types);
    track_to_ipin_lookup =
        (vector_t****)my_malloc(sizeof(vector_t***) * num_types);

    for (i = 0; i < num_types; ++i) {
        ipin_to_track_map[i] = alloc_and_load_pin_to_track_map(RECEIVER,
                                                               nodes_per_chan,
                                                               Fc_in[i],
                                                               &types[i],
                                                               perturb_ipins
                                                               [i],
                                                               directionality);
        track_to_ipin_lookup[i] =
            alloc_and_load_track_to_pin_lookup(ipin_to_track_map[i],
                                               Fc_in[i], types[i].height,
                                               types[i].num_pins,
                                               nodes_per_chan);
    }

    /* END IPIN MAP */

    /* START OPIN MAP */
    /* Create opin map lookups */
    if (BI_DIRECTIONAL == directionality) {
        opin_to_track_map =
            (int*****)my_malloc(sizeof(int****) * num_types);

        for (i = 0; i < num_types; ++i) {
            opin_to_track_map[i] =
                alloc_and_load_pin_to_track_map(DRIVER,
                                                nodes_per_chan,
                                                Fc_out[i], &types[i],
                                                FALSE,
                                                directionality);
        }
    }

    /* END OPIN MAP */

    /* UDSD Modifications by WMF begin */
    /* I'm adding 2 new fields to rr_node_t, and I want them initialized to 0. */
    for (i = 0; i < num_rr_nodes; i++) {
        rr_node[i].num_wire_drivers = 0;
        rr_node[i].num_opin_drivers = 0;
    }

    alloc_and_load_rr_graph(num_rr_nodes, rr_node, num_seg_types,
                            seg_details, rr_edge_done,
                            track_to_ipin_lookup, opin_to_track_map,
                            switch_block_conn, grid, num_grid_columns,
                            num_grid_rows, Fs, unidir_sb_pattern, Fc_out, Fc_xofs,
                            Fc_yofs, rr_node_indices, nodes_per_chan,
                            sb_type, delayless_switch, directionality,
                            wire_to_ipin_switch, &Fc_clipped);
#if 0
#ifdef MUX_SIZE_DIST_DISPLAY

    if (UNI_DIRECTIONAL == directionality) {
        view_mux_size_distribution(rr_node_indices, nodes_per_chan,
                                   seg_details, seg_details);
    }

#endif
#endif

    /* Update rr_nodes capacities if global routing */
    if (graph_type == GRAPH_GLOBAL) {
        for (i = 0; i < num_rr_nodes; i++) {
            if (rr_node[i].type == CHANX || rr_node[i].type == CHANY) {
                rr_node[i].capacity = chan_width;
            }
        }
    }

    rr_graph_externals(timing_inf, segment_inf, num_seg_types,
                       nodes_per_chan, wire_to_ipin_switch, base_cost_type);
#ifdef CREATE_ECHO_FILES
    dump_rr_graph("rr_graph.echo");
#endif /* CREATE_ECHO_FILES */
    check_rr_graph(graph_type, num_types, types, num_grid_columns, num_grid_rows,
                   grid, nodes_per_chan, Fs, num_seg_types, num_switches, segment_inf,
                   global_route_switch, delayless_switch,
                   wire_to_ipin_switch, seg_details, Fc_in, Fc_out,
                   rr_node_indices, opin_to_track_map, ipin_to_track_map,
                   track_to_ipin_lookup, switch_block_conn, perturb_ipins);

    /* Free all temp structs */
    if (seg_details) {
        free_seg_details(seg_details, nodes_per_chan);
        seg_details = NULL;
    }

    if (Fc_in) {
        free(Fc_in);
        Fc_in = NULL;
    }

    if (Fc_out) {
        free(Fc_out);
        Fc_out = NULL;
    }

    if (perturb_ipins) {
        free(perturb_ipins);
        perturb_ipins = NULL;
    }

    if (switch_block_conn) {
        free_switch_block_conn(switch_block_conn, nodes_per_chan);
        switch_block_conn = NULL;
    }

    if (rr_edge_done) {
        free(rr_edge_done);
        rr_edge_done = NULL;
    }

    if (Fc_xofs) {
        free_matrix(Fc_xofs, 0, num_grid_rows, 0, sizeof(int));
        Fc_xofs = NULL;
    }

    if (Fc_yofs) {
        free_matrix(Fc_yofs, 0, num_grid_columns, 0, sizeof(int));
        Fc_yofs = NULL;
    }

    if (unidir_sb_pattern) {
        free_sblock_pattern_lookup(unidir_sb_pattern);
        unidir_sb_pattern = NULL;
    }

    if (opin_to_track_map) {
        for (i = 0; i < num_types; ++i) {
            free_matrix4(opin_to_track_map[i], 0, types[i].num_pins - 1,
                         0, types[i].height - 1, 0, 3, 0, sizeof(int));
        }

        free(opin_to_track_map);
    }

    free_type_pin_to_track_map(ipin_to_track_map, types, Fc_in);
    free_type_track_to_ipin_map(track_to_ipin_lookup, types, nodes_per_chan);
}


static void
rr_graph_externals(timing_info_t timing_inf,
                   segment_info_t* segment_inf,
                   int num_seg_types,
                   int nodes_per_chan,
                   int wire_to_ipin_switch,
                   router_base_cost_t base_cost_type)
{
    add_rr_graph_C_from_switches(timing_inf.C_ipin_cblock);
    alloc_and_load_rr_indexed_data(segment_inf, num_seg_types,
                                   rr_node_indices, nodes_per_chan,
                                   wire_to_ipin_switch, base_cost_type);
    alloc_net_rr_terminals();
    load_net_rr_terminals(rr_node_indices);
    alloc_and_load_rr_clb_source(rr_node_indices);
}


static boolean*
alloc_and_load_perturb_ipins(IN int nodes_per_chan,
                             IN int num_types,
                             IN int* Fc_in,
                             IN int* Fc_out,
                             IN directionality_t directionality)
{
    int i;
    double Fc_ratio;
    boolean* result = NULL;
    result = (boolean*) my_malloc(num_types * sizeof(boolean));

    if (BI_DIRECTIONAL == directionality) {
        for (i = 0; i < num_types; ++i) {
            result[i] = FALSE;

            if (Fc_in[i] > Fc_out[i]) {
                Fc_ratio = (double)Fc_in[i] / (double)Fc_out[i];
            } else {
                Fc_ratio = (double)Fc_out[i] / (double)Fc_in[i];
            }

            if ((Fc_in[i] <= nodes_per_chan - 2) &&
                    (fabs(Fc_ratio - nint(Fc_ratio)) <
                     (0.5 / (double)nodes_per_chan))) {
                result[i] = TRUE;
            }
        }
    } else {
        /* Unidirectional routing uses mux balancing patterns and
         * thus shouldn't need perturbation. */
        assert(UNI_DIRECTIONAL == directionality);

        for (i = 0; i < num_types; ++i) {
            result[i] = FALSE;
        }
    }

    return result;
}


static segment_details_t*
alloc_and_load_global_route_seg_details(IN int nodes_per_chan,
                                        IN int global_route_switch)
{
    segment_details_t* result = NULL;
    assert(nodes_per_chan == 1);
    result = (segment_details_t*) my_malloc(sizeof(segment_details_t));
    result->index = 0;
    result->length = 1;
    result->wire_switch = global_route_switch;
    result->opin_switch = global_route_switch;
    result->longline = FALSE;
    result->direction = BI_DIRECTION;
    result->Cmetal = 0.0;
    result->Rmetal = 0.0;
    result->start = 1;
    result->drivers = MULTI_BUFFERED;
    result->cb = (boolean*) my_malloc(sizeof(boolean) * 1);
    result->cb[0] = TRUE;
    result->sb = (boolean*) my_malloc(sizeof(boolean) * 2);
    result->sb[0] = TRUE;
    result->sb[1] = TRUE;
    result->group_size = 1;
    result->group_start = 0;
    return result;
}




/* Calculates the actual Fc values for the given nodes_per_chan value */
static int*
alloc_and_load_actual_fc(IN int num_types,
                         IN block_type_ptr types,
                         IN int nodes_per_chan,
                         IN boolean is_Fc_out,
                         IN directionality_t directionality,
                         OUT boolean* Fc_clipped)
{
    int i;
    int* Result = NULL;
    int fac, num_sets;
    double Fc;
    *Fc_clipped = FALSE;
    /* Unidir tracks formed in pairs, otherwise no effect. */
    fac = 1;

    if (UNI_DIRECTIONAL == directionality) {
        fac = 2;
    }

    assert((nodes_per_chan % fac) == 0);
    num_sets = nodes_per_chan / fac;
    Result = (int*)my_malloc(sizeof(int) * num_types);

    for (i = 0; i < num_types; ++i) {
        Fc = (is_Fc_out ? type_descriptors[i].
              Fc_out : type_descriptors[i].Fc_in);

        if (type_descriptors[i].is_Fc_frac) {
            Result[i] = fac * nint(num_sets * Fc);
        } else {
            Result[i] = Fc;
        }

        if (is_Fc_out && type_descriptors[i].is_Fc_out_full_flex) {
            Result[i] = nodes_per_chan;
        }

        Result[i] = max(Result[i], fac);

        if (Result[i] > nodes_per_chan) {
            *Fc_clipped = TRUE;
            Result[i] = nodes_per_chan;
        }

        assert(Result[i] % fac == 0);
    }

    return Result;
}

/* frees the track to ipin mapping for each physical grid type */
static void
free_type_track_to_ipin_map(vector_t**** track_to_pin_map, block_type_ptr types, int nodes_per_chan)
{
    int i, itrack, ioff, iside;

    for (i = 0; i < num_types; i++) {
        if (track_to_pin_map[i] != NULL) {
            for (itrack = 0; itrack < nodes_per_chan; itrack++) {
                for (ioff = 0; ioff < types[i].height; ioff++) {
                    for (iside = 0; iside < 4; iside++) {
                        if (track_to_pin_map[i][itrack][ioff][iside].list != NULL) {
                            free(track_to_pin_map[i][itrack][ioff][iside].list);
                        }
                    }
                }
            }

            free_matrix3(track_to_pin_map[i], 0, nodes_per_chan - 1, 0,
                         types[i].height - 1, 0,
                         sizeof(vector_t));
        }
    }

    free(track_to_pin_map);
}

/* frees the ipin to track mapping for each physical grid type */
static void
free_type_pin_to_track_map(int**** * ipin_to_track_map, block_type_ptr types, int* fc_in)
{
    int i;

    for (i = 0; i < num_types; i++) {
        free_matrix4(ipin_to_track_map[i], 0, types[i].num_pins - 1,
                     0, types[i].height - 1, 0, 3, 0, sizeof(int));
    }

    free(ipin_to_track_map);
}

/* Does the actual work of allocating the rr_graph and filling all the *
 * appropriate values.  Everything up to this was just a prelude!      */
static void
alloc_and_load_rr_graph(IN int num_nodes,
                        IN rr_node_t* rr_node,
                        IN int num_seg_types,
                        IN segment_details_t* seg_details,
                        IN boolean* rr_edge_done,
                        IN vector_t**** track_to_ipin_lookup,
                        IN int**** *opin_to_track_map,
                        IN vector_t** *switch_block_conn,
                        IN grid_tile_t** grid,
                        IN int num_grid_columns,
                        IN int num_grid_rows,
                        IN int Fs,
                        IN short**** *sblock_pattern,
                        IN int* Fc_out,
                        IN int** Fc_xofs,
                        IN int** Fc_yofs,
                        IN vector_t** * rr_node_indices,
                        IN int nodes_per_chan,
                        IN switch_block_t sb_type,
                        IN int delayless_switch,
                        IN directionality_t directionality,
                        IN int wire_to_ipin_switch,
                        OUT boolean* Fc_clipped)
{
    int i, j;
    boolean clipped;
    int* opin_mux_size = NULL;
    /* If Fc gets clipped, this will be flagged to true */
    *Fc_clipped = FALSE;

    /* Connection SINKS and SOURCES to their pins. */
    for (i = 0; i <= (num_grid_columns + 1); i++) {
        for (j = 0; j <= (num_grid_rows + 1); j++) {
            build_rr_sinks_sources(i, j, rr_node, rr_node_indices,
                                   delayless_switch, grid);
        }
    }

    /* Build opins */
    for (i = 0; i <= (num_grid_columns + 1); ++i) {
        for (j = 0; j <= (num_grid_rows + 1); ++j) {
            if (BI_DIRECTIONAL == directionality) {
                build_bidir_rr_opins(i, j, rr_node,
                                     rr_node_indices,
                                     opin_to_track_map, Fc_out,
                                     rr_edge_done, seg_details,
                                     grid);
            } else {
                assert(UNI_DIRECTIONAL == directionality);
                build_unidir_rr_opins(i, j, grid, Fc_out,
                                      nodes_per_chan, seg_details,
                                      Fc_xofs, Fc_yofs, rr_node,
                                      rr_edge_done, &clipped,
                                      rr_node_indices);

                if (clipped) {
                    *Fc_clipped = TRUE;
                }
            }
        }
    }

    /* We make a copy of the current fanin values for the nodes to
     * know the number of OPINs driving each mux presently */
    opin_mux_size = (int*)my_malloc(sizeof(int) * num_nodes);

    for (i = 0; i < num_nodes; ++i) {
        opin_mux_size[i] = rr_node[i].fan_in;
    }

    /* Build channels */
    assert(Fs % 3 == 0);

    for (i = 0; i <= num_grid_columns; i++) {
        for (j = 0; j <= num_grid_rows; j++) {
            if (i > 0) {
                build_rr_xchan(i, j, track_to_ipin_lookup,
                               switch_block_conn,
                               CHANX_COST_INDEX_START,
                               nodes_per_chan, opin_mux_size,
                               sblock_pattern, Fs / 3,
                               seg_details, rr_node_indices,
                               rr_edge_done, rr_node,
                               wire_to_ipin_switch,
                               directionality);
            }

            if (j > 0) {
                build_rr_ychan(i, j, track_to_ipin_lookup,
                               switch_block_conn,
                               CHANX_COST_INDEX_START +
                               num_seg_types, nodes_per_chan,
                               opin_mux_size, sblock_pattern,
                               Fs / 3, seg_details,
                               rr_node_indices, rr_edge_done,
                               rr_node, wire_to_ipin_switch,
                               directionality);
            }
        }
    }

    free(opin_mux_size);
}



static void
build_bidir_rr_opins(IN int i,
                     IN int j,
                     INOUT rr_node_t* rr_node,
                     IN vector_t** * rr_node_indices,
                     IN int**** *opin_to_track_map,
                     IN int* Fc_out,
                     IN boolean* rr_edge_done,
                     IN segment_details_t* seg_details,
                     IN grid_tile_t** grid)
{
    int ipin, inode, num_edges, Fc, ofs;
    block_type_ptr type;
    struct s_linked_edge* edge_list, *next;

    /* OPIN edges need to be done at once so let the offset 0
     * block do the work. */
    if (grid[i][j].m_offset > 0) {
        return;
    }

    type = grid[i][j].grid_type;
    Fc = Fc_out[type->index];

    for (ipin = 0; ipin < type->num_pins; ++ipin) {
        /* We only are working with opins so skip non-drivers */
        if (type->class_inf[type->pin_class[ipin]].type != DRIVER) {
            continue;
        }

        num_edges = 0;
        edge_list = NULL;

        for (ofs = 0; ofs < type->height; ++ofs) {
            num_edges +=
                get_bidir_opin_connections(i, j + ofs, ipin,
                                           &edge_list,
                                           opin_to_track_map, Fc,
                                           rr_edge_done,
                                           rr_node_indices,
                                           seg_details);
        }

        inode = get_rr_node_index(i, j, OPIN, ipin, rr_node_indices);
        alloc_and_load_edges_and_switches(rr_node, inode, num_edges,
                                          rr_edge_done, edge_list);

        while (edge_list != NULL) {
            next = edge_list->next;
            free(edge_list);
            edge_list = next;
        }
    }
}



void
free_rr_graph(void)
{
    int i;

    /* Frees all the routing graph data structures, if they have been       *
     * allocated.  I use rr_mem_chunk_list_head as a flag to indicate       *
     * whether or not the graph has been allocated -- if it is not NULL,    *
     * a routing graph exists and can be freed.  Hence, you can call this   *
     * routine even if you're not sure of whether a rr_graph exists or not. */

    if (rr_mem_chunk_list_head == NULL) { /* Nothing to free. */
        return;
    }

    free_chunk_memory(rr_mem_chunk_list_head);  /* Frees ALL "chunked" data */
    rr_mem_chunk_list_head = NULL;  /* No chunks allocated now. */
    chunk_bytes_avail = 0;  /* 0 bytes left in current "chunk". */
    chunk_next_avail_mem = NULL;    /* No current chunk.                */
    /* Before adding any more free calls here, be sure the data is NOT chunk *
     * allocated, as ALL the chunk allocated data is already free!           */
    free(net_rr_terminals);

    for (i = 0; i < num_rr_nodes; i++) {
        if (rr_node[i].edges != NULL) {
            free(rr_node[i].edges);
        }

        if (rr_node[i].switches != NULL) {
            free(rr_node[i].switches);
        }
    }

    assert(rr_node_indices);
    free_rr_node_indices(rr_node_indices);
    free(rr_node);
    free(rr_indexed_data);

    for (i = 0; i < num_blocks; i++) {
        free(rr_blk_source[i]);
    }

    free(rr_blk_source);
    rr_blk_source = NULL;
    net_rr_terminals = NULL;
    rr_node = NULL;
    rr_node_indices = NULL;
    rr_indexed_data = NULL;
}


static void
alloc_net_rr_terminals(void)
{
    int inet;
    net_rr_terminals = (int**)my_malloc(num_nets * sizeof(int*));

    for (inet = 0; inet < num_nets; inet++) {
        net_rr_terminals[inet] =
            (int*)my_chunk_malloc((net[inet].num_net_pins + 1) *
                                  sizeof(int), &rr_mem_chunk_list_head,
                                  &chunk_bytes_avail,
                                  &chunk_next_avail_mem);
    }
}


void
load_net_rr_terminals(vector_t** * rr_node_indices)
{
    /* Allocates and loads the net_rr_terminals data structure.  For each net   *
     * it stores the rr_node index of the SOURCE of the net and all the SINKs   *
     * of the net.  [0..num_nets-1][0..num_pins-1].  Entry [inet][pnum] stores  *
     * the rr index corresponding to the SOURCE (opin) or SINK (ipin) of pnum.  */
    int inet, ipin, inode, iblk, i, j, node_block_pin, iclass;
    block_type_ptr type;

    for (inet = 0; inet < num_nets; ++inet) {
        const int knum_net_pins = net[inet].num_net_pins;

        for (ipin = 0; ipin <= knum_net_pins; ++ipin) {
            iblk = net[inet].node_blocks[ipin];
            i = blocks[iblk].x;
            j = blocks[iblk].y;
            type = blocks[iblk].block_type;
            /* In the routing graph, each (x, y) location has unique pins on it
             * so when there is capacity, blocks are packed and their pin numbers
             * are offset to get their actual rr_node */
            node_block_pin = net[inet].node_block_pins[ipin];
            iclass = type->pin_class[node_block_pin];
            inode = get_rr_node_index(i, j, (ipin == 0 ? SOURCE : SINK),    /* First pin is driver */
                                      iclass, rr_node_indices);
            net_rr_terminals[inet][ipin] = inode;
        }
    }
}

static void
alloc_and_load_rr_clb_source(vector_t** * rr_node_indices)
{
    /* Saves the rr_node corresponding to each SOURCE and SINK in each FB      *
     * in the FPGA.  Currently only the SOURCE rr_node values are used, and     *
     * they are used only to reserve pins for locally used OPINs in the router. *
     * [0..num_blocks-1][0..num_class-1].  The values for blocks that are pads  *
     * are NOT valid.                                                           */
    int iblk, i, j, iclass, inode;
    int class_low, class_high;
    rr_type_t rr_type;
    block_type_ptr type;
    rr_blk_source = (int**)my_malloc(num_blocks * sizeof(int*));

    for (iblk = 0; iblk < num_blocks; iblk++) {
        type = blocks[iblk].block_type;
        get_class_range_for_block(iblk, &class_low, &class_high);
        rr_blk_source[iblk] =
            (int*)my_malloc(type->num_class * sizeof(int));

        for (iclass = 0; iclass < type->num_class; iclass++) {
            if (iclass >= class_low && iclass <= class_high) {
                i = blocks[iblk].x;
                j = blocks[iblk].y;

                if (type->class_inf[iclass].type == DRIVER) {
                    rr_type = SOURCE;
                } else {
                    rr_type = SINK;
                }

                inode =
                    get_rr_node_index(i, j, rr_type, iclass,
                                      rr_node_indices);
                rr_blk_source[iblk][iclass] = inode;
            } else {
                rr_blk_source[iblk][iclass] = OPEN;
            }
        }
    }
}

static void
build_rr_sinks_sources(IN int i,
                       IN int j,
                       IN rr_node_t* rr_node,
                       IN vector_t** * rr_node_indices,
                       IN int delayless_switch,
                       IN grid_tile_t** grid)
{
    /* Loads IPIN, SINK, SOURCE, and OPIN.
     * Loads IPIN to SINK edges, and SOURCE to OPIN edges */
    int ipin, iclass, inode, pin_num, to_node, num_edges;
    int num_class, num_pins;
    pin_class_t* class_inf;
    int* pin_class;

    /* Since we share nodes within a large block, only
     * start tile can initialize sinks, sources, and pins */
    if (grid[i][j].m_offset > 0) {
        return;
    }

    block_type_ptr type = grid[i][j].grid_type;
    num_class = type->num_class;
    class_inf = type->class_inf;
    num_pins = type->num_pins;
    pin_class = type->pin_class;

    /* SINKS and SOURCE to OPIN edges */
    for (iclass = 0; iclass < num_class; iclass++) {
        if (class_inf[iclass].type == DRIVER) {
            /* SOURCE */
            inode =
                get_rr_node_index(i, j, SOURCE, iclass,
                                  rr_node_indices);
            num_edges = class_inf[iclass].num_pins;
            rr_node[inode].num_edges = num_edges;
            rr_node[inode].edges =
                (int*)my_malloc(num_edges * sizeof(int));
            rr_node[inode].switches =
                (short*)my_malloc(num_edges * sizeof(short));

            for (ipin = 0; ipin < class_inf[iclass].num_pins; ipin++) {
                pin_num = class_inf[iclass].pinlist[ipin];
                to_node =
                    get_rr_node_index(i, j, OPIN, pin_num,
                                      rr_node_indices);
                rr_node[inode].edges[ipin] = to_node;
                rr_node[inode].switches[ipin] = delayless_switch;
                ++rr_node[to_node].fan_in;
            }

            rr_node[inode].cost_index = SOURCE_COST_INDEX;
            rr_node[inode].type = SOURCE;
        } else {
            /* SINK */
            assert(class_inf[iclass].type == RECEIVER);
            inode =
                get_rr_node_index(i, j, SINK, iclass,
                                  rr_node_indices);
            /* NOTE:  To allow route throughs through clbs, change the lines below to  *
             * make an tedge from the input SINK to the output SOURCE.  Do for just the *
             * special case of INPUTS = class 0 and OUTPUTS = class 1 and see what it  *
             * leads to.  If route throughs are allowed, you may want to increase the  *
             * base cost of OPINs and/or SOURCES so they aren't used excessively.      */
            /* Initialize to unconnected to fix values */
            rr_node[inode].num_edges = 0;
            rr_node[inode].edges = NULL;
            rr_node[inode].switches = NULL;
            rr_node[inode].cost_index = SINK_COST_INDEX;
            rr_node[inode].type = SINK;
        }

        /* Things common to both SOURCEs and SINKs.   */
        rr_node[inode].capacity = class_inf[iclass].num_pins;
        rr_node[inode].occ = 0;
        rr_node[inode].xlow = i;
        rr_node[inode].xhigh = i;
        rr_node[inode].ylow = j;
        rr_node[inode].yhigh = j + type->height - 1;
        rr_node[inode].R = 0;
        rr_node[inode].C = 0;
        rr_node[inode].ptc_num = iclass;
        rr_node[inode].direction = OPEN;
        rr_node[inode].drivers = OPEN;
    }

    /* Connect IPINS to SINKS and dummy for OPINS */
    for (ipin = 0; ipin < num_pins; ipin++) {
        iclass = pin_class[ipin];

        if (class_inf[iclass].type == RECEIVER) {
            inode =
                get_rr_node_index(i, j, IPIN, ipin, rr_node_indices);
            to_node =
                get_rr_node_index(i, j, SINK, iclass,
                                  rr_node_indices);
            rr_node[inode].num_edges = 1;
            rr_node[inode].edges = (int*)my_malloc(sizeof(int));
            rr_node[inode].switches =
                (short*)my_malloc(sizeof(short));
            rr_node[inode].edges[0] = to_node;
            rr_node[inode].switches[0] = delayless_switch;
            ++rr_node[to_node].fan_in;
            rr_node[inode].cost_index = IPIN_COST_INDEX;
            rr_node[inode].type = IPIN;
        } else {
            assert(class_inf[iclass].type == DRIVER);
            inode =
                get_rr_node_index(i, j, OPIN, ipin, rr_node_indices);
            rr_node[inode].num_edges = 0;
            rr_node[inode].edges = NULL;
            rr_node[inode].switches = NULL;
            rr_node[inode].cost_index = OPIN_COST_INDEX;
            rr_node[inode].type = OPIN;
        }

        /* Common to both DRIVERs and RECEIVERs */
        rr_node[inode].capacity = 1;
        rr_node[inode].occ = 0;
        rr_node[inode].xlow = i;
        rr_node[inode].xhigh = i;
        rr_node[inode].ylow = j;
        rr_node[inode].yhigh = j + type->height - 1;
        rr_node[inode].C = 0;
        rr_node[inode].R = 0;
        rr_node[inode].ptc_num = ipin;
        rr_node[inode].direction = OPEN;
        rr_node[inode].drivers = OPEN;
    }
}



static void
build_rr_xchan(IN int i,
               IN int j,
               IN vector_t**** track_to_ipin_lookup,
               IN vector_t** *switch_block_conn,
               IN int cost_index_offset,
               IN int nodes_per_chan,
               IN int* opin_mux_size,
               IN short**** *sblock_pattern,
               IN int Fs_per_side,
               IN segment_details_t* seg_details,
               IN vector_t** * rr_node_indices,
               INOUT boolean* rr_edge_done,
               INOUT rr_node_t* rr_node,
               IN int wire_to_ipin_switch,
               IN directionality_t directionality)
{
    /* Loads up all the routing resource nodes in the x-directed channel      *
     * segments starting at (i,j).                                            */
    int itrack, istart, iend, num_edges, inode, length;
    struct s_linked_edge* edge_list, *next;

    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        istart = get_seg_start(seg_details, itrack, j, i);
        iend = get_seg_end(seg_details, itrack, istart, j, num_grid_columns);

        if (i > istart) {
            continue;    /* Not the start of this segment. */
        }

        edge_list = NULL;
        /* First count number of edges and put the edges in a linked list. */
        num_edges = 0;
        num_edges += get_track_to_ipins(istart, j, itrack,
                                        &edge_list,
                                        rr_node_indices,
                                        track_to_ipin_lookup,
                                        seg_details,
                                        CHANX, num_grid_columns,
                                        wire_to_ipin_switch,
                                        directionality);

        if (j > 0) {
            num_edges += get_track_to_tracks(j, istart, itrack, CHANX,
                                             j, CHANY, num_grid_columns,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        if (j < num_grid_rows) {
            num_edges += get_track_to_tracks(j, istart, itrack, CHANX,
                                             j + 1, CHANY, num_grid_columns,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        if (istart > 1) {
            num_edges += get_track_to_tracks(j, istart, itrack, CHANX,
                                             istart - 1, CHANX, num_grid_columns,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        if (iend < num_grid_columns) {
            num_edges += get_track_to_tracks(j, istart, itrack, CHANX,
                                             iend + 1, CHANX, num_grid_columns,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        inode = get_rr_node_index(i, j, CHANX, itrack, rr_node_indices);
        alloc_and_load_edges_and_switches(rr_node, inode, num_edges,
                                          rr_edge_done, edge_list);

        while (edge_list != NULL) {
            next = edge_list->next;
            free(edge_list);
            edge_list = next;
        }

        /* Edge arrays have now been built up.  Do everything else.  */
        rr_node[inode].cost_index =
            cost_index_offset + seg_details[itrack].index;
        rr_node[inode].occ = 0;
        rr_node[inode].capacity = 1;    /* GLOBAL routing handled elsewhere */
        rr_node[inode].xlow = istart;
        rr_node[inode].xhigh = iend;
        rr_node[inode].ylow = j;
        rr_node[inode].yhigh = j;
        length = iend - istart + 1;
        rr_node[inode].R = length * seg_details[itrack].Rmetal;
        rr_node[inode].C = length * seg_details[itrack].Cmetal;
        rr_node[inode].ptc_num = itrack;
        rr_node[inode].type = CHANX;
        rr_node[inode].direction = seg_details[itrack].direction;
        rr_node[inode].drivers = seg_details[itrack].drivers;
    }
}




static void
build_rr_ychan(IN int i,
               IN int j,
               IN vector_t**** track_to_ipin_lookup,
               IN vector_t** *switch_block_conn,
               IN int cost_index_offset,
               IN int nodes_per_chan,
               IN int* opin_mux_size,
               IN short**** *sblock_pattern,
               IN int Fs_per_side,
               IN segment_details_t* seg_details,
               IN vector_t** * rr_node_indices,
               IN boolean* rr_edge_done,
               INOUT rr_node_t* rr_node,
               IN int wire_to_ipin_switch,
               IN directionality_t directionality)
{
    /* Loads up all the routing resource nodes in the y-directed channel      *
     * segments starting at (i,j).                                            */
    int itrack, istart, iend, num_edges, inode, length;
    struct s_linked_edge* edge_list, *next;

    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        istart = get_seg_start(seg_details, itrack, i, j);
        iend = get_seg_end(seg_details, itrack, istart, i, num_grid_rows);

        if (j > istart) {
            continue;    /* Not the start of this segment. */
        }

        edge_list = NULL;
        /* First count number of edges and put the edges in a linked list. */
        num_edges = 0;
        num_edges += get_track_to_ipins(istart, i, itrack,
                                        &edge_list,
                                        rr_node_indices,
                                        track_to_ipin_lookup,
                                        seg_details,
                                        CHANY, num_grid_rows,
                                        wire_to_ipin_switch,
                                        directionality);

        if (i > 0) {
            num_edges += get_track_to_tracks(i, istart, itrack, CHANY,
                                             i, CHANX, num_grid_rows,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        if (i < num_grid_columns) {
            num_edges += get_track_to_tracks(i, istart, itrack, CHANY,
                                             i + 1, CHANX, num_grid_rows,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        if (istart > 1) {
            num_edges += get_track_to_tracks(i, istart, itrack, CHANY,
                                             istart - 1, CHANY, num_grid_rows,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        if (iend < num_grid_rows) {
            num_edges += get_track_to_tracks(i, istart, itrack, CHANY,
                                             iend + 1, CHANY, num_grid_rows,
                                             nodes_per_chan,
                                             opin_mux_size,
                                             Fs_per_side,
                                             sblock_pattern,
                                             &edge_list,
                                             seg_details,
                                             directionality,
                                             rr_node_indices,
                                             rr_edge_done,
                                             switch_block_conn);
        }

        inode = get_rr_node_index(i, j, CHANY, itrack, rr_node_indices);
        alloc_and_load_edges_and_switches(rr_node, inode, num_edges,
                                          rr_edge_done, edge_list);

        while (edge_list != NULL) {
            next = edge_list->next;
            free(edge_list);
            edge_list = next;
        }

        /* Edge arrays have now been built up.  Do everything else.  */
        rr_node[inode].cost_index =
            cost_index_offset + seg_details[itrack].index;
        rr_node[inode].occ = 0;
        rr_node[inode].capacity = 1;    /* GLOBAL routing handled elsewhere */
        rr_node[inode].xlow = i;
        rr_node[inode].xhigh = i;
        rr_node[inode].ylow = istart;
        rr_node[inode].yhigh = iend;
        length = iend - istart + 1;
        rr_node[inode].R = length * seg_details[itrack].Rmetal;
        rr_node[inode].C = length * seg_details[itrack].Cmetal;
        rr_node[inode].ptc_num = itrack;
        rr_node[inode].type = CHANY;
        rr_node[inode].direction = seg_details[itrack].direction;
        rr_node[inode].drivers = seg_details[itrack].drivers;
    }
}

void
watch_edges(int inode,
            t_linked_edge* edge_list_head)
{
    t_linked_edge* list_ptr;
    int i, to_node;
    list_ptr = edge_list_head;
    i = 0;
    printf("!!! Watching Node %d !!!!\n", inode);
    print_rr_node(stdout, rr_node, inode);
    printf("Currently connects to: \n");

    while (list_ptr != NULL) {
        to_node = list_ptr->tedge;
        print_rr_node(stdout, rr_node, to_node);
        list_ptr = list_ptr->next;
        i++;
    }
}

void
alloc_and_load_edges_and_switches(IN rr_node_t* rr_node,
                                  IN int inode,
                                  IN int num_edges,
                                  INOUT boolean* rr_edge_done,
                                  IN t_linked_edge* edge_list_head)
{
    /* Sets up all the tedge related information for rr_node inode (num_edges,  *
     * the edges array and the switches array).  The edge_list_head points to  *
     * a list of the num_edges edges and switches to put in the arrays.  This  *
     * linked list is freed by this routine. This routine also resets the      *
     * rr_edge_done array for the next rr_node (i.e. set it so that no edges   *
     * are marked as having been seen before).                                 */
    t_linked_edge* list_ptr;
    int i;
    /* Check we aren't overwriting edges */
    assert(rr_node[inode].num_edges < 1);
    assert(NULL == rr_node[inode].edges);
    assert(NULL == rr_node[inode].switches);
    rr_node[inode].num_edges = num_edges;
    rr_node[inode].edges = (int*)my_malloc(num_edges * sizeof(int));
    rr_node[inode].switches = (short*)my_malloc(num_edges * sizeof(short));
    i = 0;
    list_ptr = edge_list_head;

    while (list_ptr && (i < num_edges)) {
        rr_node[inode].edges[i] = list_ptr->tedge;
        rr_node[inode].switches[i] = list_ptr->iswitch;
        ++rr_node[list_ptr->tedge].fan_in;
        /* Unmark the tedge since we are done considering fanout from node. */
        rr_edge_done[list_ptr->tedge] = FALSE;
        list_ptr = list_ptr->next;
        ++i;
    }

    assert(list_ptr == NULL);
    assert(i == num_edges);
}







static int****
alloc_and_load_pin_to_track_map(IN pin_types_t pin_type,
                                IN int nodes_per_chan,
                                IN int Fc,
                                IN block_type_ptr Type,
                                IN boolean perturb_switch_pattern,
                                IN directionality_t directionality)
{
    int** num_dir;      /* [0..height][0..3] Number of *physical* pins on each side.          */
    int** *dir_list;        /* [0..height][0..3][0..num_pins-1] list of pins of correct type  *
                 * * on each side. Max possible space alloced for simplicity */
    int i, j, k, iside, ipin, iclass, num_phys_pins, pindex, ioff;
    int* pin_num_ordering, *side_ordering, *offset_ordering;
    int** num_done_per_dir; /* [0..height][0..3] */
    int**** tracks_connected_to_pin; /* [0..num_pins-1][0..height][0..3][0..Fc-1] */

    /* NB:  This wastes some space.  Could set tracks_..._pin[ipin][ioff][iside] =
     * NULL if there is no pin on that side, or that pin is of the wrong type.
     * Probably not enough memory to worry about, esp. as it's temporary.
     * If pin ipin on side iside does not exist or is of the wrong type,
     * tracks_connected_to_pin[ipin][iside][0] = OPEN.                               */

    if (Type->num_pins < 1) {
        return NULL;
    }

    tracks_connected_to_pin = (int****)
                              alloc_matrix4(0, Type->num_pins - 1,
                                            0, Type->height - 1, 0, 3, 0, Fc - 1, sizeof(int));

    for (ipin = 0; ipin < Type->num_pins; ipin++) {
        for (ioff = 0; ioff < Type->height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                for (i = 0; i < Fc; ++i) {
                    tracks_connected_to_pin[ipin][ioff][iside][i] = OPEN;   /* Unconnected. */
                }
            }
        }
    }

    num_dir = (int**)alloc_matrix(0, Type->height - 1, 0, 3, sizeof(int));
    dir_list =
        (int***)alloc_matrix3(0, Type->height - 1, 0, 3, 0,
                              Type->num_pins - 1, sizeof(int));

    /* Defensive coding.  Try to crash hard if I use an unset entry.  */
    for (i = 0; i < Type->height; i++)
        for (j = 0; j < 4; j++)
            for (k = 0; k < Type->num_pins; k++) {
                dir_list[i][j][k] = (-1);
            }

    for (i = 0; i < Type->height; i++)
        for (j = 0; j < 4; j++) {
            num_dir[i][j] = 0;
        }

    for (ipin = 0; ipin < Type->num_pins; ipin++) {
        iclass = Type->pin_class[ipin];

        if (Type->class_inf[iclass].type != pin_type) { /* Doing either ipins OR opins */
            continue;
        }

        /* Pins connecting only to global resources get no switches->keeps the    *
         * area model accurate.                                                     */

        if (Type->is_global_pin[ipin]) {
            continue;
        }

        for (ioff = 0; ioff < Type->height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                if (Type->pinloc[ioff][iside][ipin] == 1) {
                    dir_list[ioff][iside][num_dir[ioff]
                                          [iside]] = ipin;
                    num_dir[ioff][iside]++;
                }
            }
        }
    }

    num_phys_pins = 0;

    for (ioff = 0; ioff < Type->height; ioff++) {
        for (iside = 0; iside < 4; iside++) {
            num_phys_pins += num_dir[ioff][iside];    /* Num. physical pins per type */
        }
    }

    num_done_per_dir =
        (int**)alloc_matrix(0, Type->height - 1, 0, 3, sizeof(int));

    for (ioff = 0; ioff < Type->height; ioff++) {
        for (iside = 0; iside < 4; iside++) {
            num_done_per_dir[ioff][iside] = 0;
        }
    }

    pin_num_ordering = (int*)my_malloc(num_phys_pins * sizeof(int));
    side_ordering = (int*)my_malloc(num_phys_pins * sizeof(int));
    offset_ordering = (int*)my_malloc(num_phys_pins * sizeof(int));
    /* Connection block I use distributes pins evenly across the tracks      *
     * of ALL sides of the clb at once.  Ensures that each pin connects      *
     * to spaced out tracks in its connection block, and that the other      *
     * pins (potentially in other C blocks) connect to the remaining tracks  *
     * first.  Doesn't matter for large Fc, but should make a fairly         *
     * good low Fc block that leverages the fact that usually lots of pins   *
     * are logically equivalent.                                             */
    iside = LEFT;
    ioff = Type->height - 1;
    ipin = 0;
    pindex = -1;

    while (ipin < num_phys_pins) {
        if (iside == TOP) {
            iside = RIGHT;
        } else if (iside == RIGHT) {
            if (ioff <= 0) {
                iside = BOTTOM;
            } else {
                ioff--;
            }
        } else if (iside == BOTTOM) {
            iside = LEFT;
        } else {
            assert(iside == LEFT);

            if (ioff >= Type->height - 1) {
                pindex++;
                iside = TOP;
            } else {
                ioff++;
            }
        }

        assert(pindex < num_phys_pins); /* Number of physical pins bounds number of logical pins */

        if (num_done_per_dir[ioff][iside] >= num_dir[ioff][iside]) {
            continue;
        }

        pin_num_ordering[ipin] = dir_list[ioff][iside][pindex];
        side_ordering[ipin] = iside;
        offset_ordering[ipin] = ioff;
        assert(Type->pinloc[ioff][iside][dir_list[ioff][iside][pindex]]);
        num_done_per_dir[ioff][iside]++;
        ipin++;
    }

    if (perturb_switch_pattern) {
        load_perturbed_switch_pattern(Type,
                                      tracks_connected_to_pin,
                                      num_phys_pins,
                                      pin_num_ordering,
                                      side_ordering,
                                      offset_ordering,
                                      nodes_per_chan, Fc, directionality);
    } else {
        load_uniform_switch_pattern(Type,
                                    tracks_connected_to_pin,
                                    num_phys_pins,
                                    pin_num_ordering,
                                    side_ordering,
                                    offset_ordering,
                                    nodes_per_chan, Fc, directionality);
    }

    check_all_tracks_reach_pins(Type, tracks_connected_to_pin,
                                nodes_per_chan, Fc, pin_type);
    /* Free all temporary storage. */
    free_matrix(num_dir, 0, Type->height - 1, 0, sizeof(int));
    free_matrix3(dir_list, 0, Type->height - 1, 0, 3, 0, sizeof(int));
    free_matrix(num_done_per_dir, 0, Type->height - 1, 0, sizeof(int));
    free(pin_num_ordering);
    free(side_ordering);
    free(offset_ordering);
    return tracks_connected_to_pin;
}




static void
load_uniform_switch_pattern(IN block_type_ptr type,
                            INOUT int**** tracks_connected_to_pin,
                            IN int num_phys_pins,
                            IN int* pin_num_ordering,
                            IN int* side_ordering,
                            IN int* offset_ordering,
                            IN int nodes_per_chan,
                            IN int Fc,
                            directionality_t directionality)
{
    /* Loads the tracks_connected_to_pin array with an even distribution of     *
     * switches across the tracks for each pin.  For example, each pin connects *
     * to every 4.3rd track in a channel, with exactly which tracks a pin       *
     * connects to staggered from pin to pin.                                   */
    int i, j, ipin, iside, ioff, itrack, k;
    double f_track, fc_step;
    int group_size;
    double step_size;
    /* Uni-directional drive is implemented to ensure no directional bias and this means
     * two important comments noted below                                                */
    /* 1. Spacing should be (W/2)/(Fc/2), and step_size should be spacing/(num_phys_pins),
     *    and lay down 2 switches on an adjacent pair of tracks at a time to ensure
     *    no directional bias. Basically, treat W (even) as W/2 pairs of tracks, and
     *    assign switches to a pair at a time. Can do this because W is guaranteed to
     *    be even-numbered; however same approach cannot be applied to Fc_out pattern
     *    when L > 1 and W <> 2L multiple.
     *
     * 2. This generic pattern should be considered the tileable physical layout,
     *    meaning all track # here are physical #'s,
     *    so later must use vpr_to_phy conversion to find actual logical #'s to connect.
     *    This also means I will not use get_output_block_companion_track to ensure
     *    no bias, since that describes a logical #->that would confuse people.  */
    step_size = (double)nodes_per_chan / (double)(Fc * num_phys_pins);

    if (directionality == BI_DIRECTIONAL) {
        group_size = 1;
    } else {
        assert(directionality == UNI_DIRECTIONAL);
        group_size = 2;
    }

    assert((nodes_per_chan % group_size == 0) && (Fc % group_size == 0));
    fc_step = (double)nodes_per_chan / (double)Fc;

    for (i = 0; i < num_phys_pins; i++) {
        ipin = pin_num_ordering[i];
        iside = side_ordering[i];
        ioff = offset_ordering[i];

        /* Bi-directional treats each track separately, uni-directional works with pairs of tracks */
        for (j = 0; j < (Fc / group_size); j++) {
            f_track = (i * step_size) + (j * fc_step);
            itrack = ((int)f_track) * group_size;
            /* Catch possible floating point round error */
            itrack = min(itrack, nodes_per_chan - group_size);

            /* Assign the group of tracks for the Fc pattern */
            for (k = 0; k < group_size; ++k) {
                tracks_connected_to_pin[ipin][ioff][iside]
                [group_size * j + k] = itrack + k;
            }
        }
    }
}

static void
load_perturbed_switch_pattern(IN block_type_ptr type,
                              INOUT int**** tracks_connected_to_pin,
                              IN int num_phys_pins,
                              IN int* pin_num_ordering,
                              IN int* side_ordering,
                              IN int* offset_ordering,
                              IN int nodes_per_chan,
                              IN int Fc,
                              directionality_t directionality)
{
    /* Loads the tracks_connected_to_pin array with an unevenly distributed     *
     * set of switches across the channel.  This is done for inputs when        *
     * Fc_input = Fc_output to avoid creating "pin domains" -- certain output   *
     * pins being able to talk only to certain input pins because their switch  *
     * patterns exactly line up.  Distribute Fc/2 + 1 switches over half the    *
     * channel and Fc/2 - 1 switches over the other half to make the switch     *
     * pattern different from the uniform one of the outputs.  Also, have half  *
     * the pins put the "dense" part of their connections in the first half of  *
     * the channel and the other half put the "dense" part in the second half,  *
     * to make sure each track can connect to about the same number of ipins.   */
    int i, j, ipin, iside, itrack, ihalf, iconn, ioff;
    int Fc_dense, Fc_sparse, Fc_half[2];
    double f_track, spacing_dense, spacing_sparse, spacing[2];
    double step_size;
    assert(directionality == BI_DIRECTIONAL);
    step_size = (double)nodes_per_chan / (double)(Fc * num_phys_pins);
    Fc_dense = (Fc / 2) + 1;
    Fc_sparse = Fc - Fc_dense;  /* Works for even or odd Fc */
    spacing_dense = (double)nodes_per_chan / (double)(2 * Fc_dense);
    spacing_sparse = (double)nodes_per_chan / (double)(2 * Fc_sparse);

    for (i = 0; i < num_phys_pins; i++) {
        ipin = pin_num_ordering[i];
        iside = side_ordering[i];
        ioff = offset_ordering[i];
        /* Flip every pin to balance switch density */
        spacing[i % 2] = spacing_dense;
        Fc_half[i % 2] = Fc_dense;
        spacing[(i + 1) % 2] = spacing_sparse;
        Fc_half[(i + 1) % 2] = Fc_sparse;
        f_track = i * step_size;    /* Start point.  Staggered from pin to pin */
        iconn = 0;

        for (ihalf = 0; ihalf < 2; ihalf++) {
            /* For both dense and sparse halves. */
            for (j = 0; j < Fc_half[ihalf]; ++j) {
                itrack = (int)f_track;
                /* Can occasionally get wraparound due to floating point rounding.
                   This is okay because the starting position > 0 when this occurs
                   so connection is valid and fine */
                itrack = itrack % nodes_per_chan;
                tracks_connected_to_pin[ipin][ioff][iside][iconn]
                    = itrack;
                f_track += spacing[ihalf];
                iconn++;
            }
        }
    }           /* End for all physical pins. */
}


static void
check_all_tracks_reach_pins(block_type_ptr type,
                            int**** tracks_connected_to_pin,
                            int nodes_per_chan,
                            int Fc,
                            pin_types_t ipin_or_opin)
{
    /* Checks that all tracks can be reached by some pin.   */
    int iconn, iside, itrack, ipin, ioff;
    int* num_conns_to_track;    /* [0..nodes_per_chan-1] */
    assert(nodes_per_chan > 0);
    num_conns_to_track = (int*)my_calloc(nodes_per_chan, sizeof(int));

    for (ipin = 0; ipin < type->num_pins; ipin++) {
        for (ioff = 0; ioff < type->height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                if (tracks_connected_to_pin[ipin][ioff][iside][0]
                        != OPEN) {
                    /* Pin exists */
                    for (iconn = 0; iconn < Fc; iconn++) {
                        itrack =
                            tracks_connected_to_pin[ipin]
                            [ioff][iside][iconn];
                        num_conns_to_track[itrack]++;
                    }
                }
            }
        }
    }

    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        if (num_conns_to_track[itrack] <= 0) {
            printf
            ("Warning (check_all_tracks_reach_pins):  track %d does not \n"
             "\tconnect to any FB ", itrack);

            if (ipin_or_opin == DRIVER) {
                printf("OPINs.\n");
            } else {
                printf("IPINs.\n");
            }
        }
    }

    free(num_conns_to_track);
}

/* Allocates and loads the track to ipin lookup for each physical grid type. This
 * is the same information as the ipin_to_track map but accessed in a different way. */

static vector_t***
alloc_and_load_track_to_pin_lookup(IN int**** pin_to_track_map,
                                   IN int Fc,
                                   IN int height,
                                   IN int num_pins,
                                   IN int nodes_per_chan) {
    int ipin, iside, itrack, iconn, ioff, pin_counter;
    vector_t*** track_to_pin_lookup;

    /* [0..nodes_per_chan-1][0..height][0..3].  For each track number it stores a vector_t
     * for each of the four sides.  x-directed channels will use the TOP and
     * BOTTOM vectors to figure out what clb input pins they connect to above
     * and below them, respectively, while y-directed channels use the LEFT
     * and RIGHT vectors.  Each vector_t contains an nelem field saying how many
     * ipins it connects to.  The list[0..nelem-1] array then gives the pin
     * numbers.                                                                */

    /* Note that a clb pin that connects to a channel on its RIGHT means that  *
     * that channel connects to a clb pin on its LEFT.  The convention used    *
     * here is always in the perspective of the CLB                            */

    if (num_pins < 1) {
        return NULL;
    }

    /* Alloc and zero the the lookup table */
    track_to_pin_lookup =
        (vector_t***)alloc_matrix3(0, nodes_per_chan - 1, 0,
                                        height - 1, 0, 3,
                                        sizeof(vector_t));

    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        for (ioff = 0; ioff < height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                track_to_pin_lookup[itrack][ioff][iside].nelem =
                    0;
                track_to_pin_lookup[itrack][ioff][iside].list =
                    NULL;
            }
        }
    }

    /* Counting pass.  */
    for (ipin = 0; ipin < num_pins; ipin++) {
        for (ioff = 0; ioff < height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                if (pin_to_track_map[ipin][ioff][iside][0] == OPEN) {
                    continue;
                }

                for (iconn = 0; iconn < Fc; iconn++) {
                    itrack =
                        pin_to_track_map[ipin][ioff][iside]
                        [iconn];
                    track_to_pin_lookup[itrack][ioff][iside].
                    nelem++;
                }
            }
        }
    }

    /* Allocate space.  */
    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        for (ioff = 0; ioff < height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                track_to_pin_lookup[itrack][ioff][iside].list = NULL;   /* Defensive code */

                if (track_to_pin_lookup[itrack][ioff][iside].
                        nelem != 0) {
                    track_to_pin_lookup[itrack][ioff][iside].
                    list =
                        (int*)
                        my_malloc(track_to_pin_lookup[itrack]
                                  [ioff][iside].nelem *
                                  sizeof(int));
                    track_to_pin_lookup[itrack][ioff][iside].
                    nelem = 0;
                }
            }
        }
    }

    /* Loading pass. */
    for (ipin = 0; ipin < num_pins; ipin++) {
        for (ioff = 0; ioff < height; ioff++) {
            for (iside = 0; iside < 4; iside++) {
                if (pin_to_track_map[ipin][ioff][iside][0] == OPEN) {
                    continue;
                }

                for (iconn = 0; iconn < Fc; iconn++) {
                    itrack =
                        pin_to_track_map[ipin][ioff][iside]
                        [iconn];
                    pin_counter =
                        track_to_pin_lookup[itrack][ioff]
                        [iside].nelem;
                    track_to_pin_lookup[itrack][ioff][iside].
                    list[pin_counter] = ipin;
                    track_to_pin_lookup[itrack][ioff][iside].
                    nelem++;
                }
            }
        }
    }

    return track_to_pin_lookup;
}

/* A utility routine to dump the contents of the routing resource graph   *
 * (everything -- connectivity, occupancy, cost, etc.) into a file.  Used *
 * only for debugging.                                                    */
void
dump_rr_graph(IN const char* file_name)
{
    int inode;
    FILE* fp;
    fp = my_fopen(file_name, "w");

    for (inode = 0; inode < num_rr_nodes; inode++) {
        print_rr_node(fp, rr_node, inode);
        fprintf(fp, "\n");
    }

#if 0
    fprintf(fp, "\n\n%d rr_indexed_data entries.\n\n", num_rr_indexed_data);

    for (index = 0; index < num_rr_indexed_data; index++) {
        print_rr_indexed_data(fp, index);
        fprintf(fp, "\n");
    }

#endif
    fclose(fp);
}


/* Prints all the data about node inode to file fp.                    */
static void
print_rr_node(FILE* fp,
              rr_node_t* rr_node,
              int inode)
{
    static const char* name_type[] = {
        "SOURCE", "SINK", "IPIN", "OPIN", "CHANX", "CHANY"
    };
    static const char* direction_name[] = {
        "OPEN", "INC_DIRECTION", "DEC_DIRECTION", "BI_DIRECTION"
    };
    static const char* drivers_name[] = {
        "OPEN", "MULTI_BUFFER", "MULTI_MUXED", "MULTI_MERGED", "SINGLE"
    };
    rr_type_t rr_type;
    int iconn;
    rr_type = rr_node[inode].type;
    /* Make sure we don't overrun const arrays */
    assert(rr_type < (sizeof(name_type) / sizeof(char*)));
    assert((rr_node[inode].direction + 1) <
           (sizeof(direction_name) / sizeof(char*)));
    assert((rr_node[inode].drivers + 1) <
           (sizeof(drivers_name) / sizeof(char*)));
    fprintf(fp, "Node: %d %s ", inode, name_type[rr_type]);

    if ((rr_node[inode].xlow == rr_node[inode].xhigh) &&
            (rr_node[inode].ylow == rr_node[inode].yhigh)) {
        fprintf(fp, "(%d, %d) ", rr_node[inode].xlow,
                rr_node[inode].ylow);
    } else {
        fprintf(fp, "(%d, %d) to (%d, %d) ", rr_node[inode].xlow,
                rr_node[inode].ylow, rr_node[inode].xhigh,
                rr_node[inode].yhigh);
    }

    fprintf(fp, "Ptc_num: %d ", rr_node[inode].ptc_num);
    fprintf(fp, "Direction: %s ",
            direction_name[rr_node[inode].direction + 1]);
    fprintf(fp, "Drivers: %s ", drivers_name[rr_node[inode].drivers + 1]);
    fprintf(fp, "\n");
    fprintf(fp, "%d tedge(s):", rr_node[inode].num_edges);

    for (iconn = 0; iconn < rr_node[inode].num_edges; iconn++) {
        fprintf(fp, " %d", rr_node[inode].edges[iconn]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "Switch types:");

    for (iconn = 0; iconn < rr_node[inode].num_edges; iconn++) {
        fprintf(fp, " %d", rr_node[inode].switches[iconn]);
    }

    fprintf(fp, "\n");
    fprintf(fp, "Occ: %d  Capacity: %d\n", rr_node[inode].occ,
            rr_node[inode].capacity);
    fprintf(fp, "R: %g  C: %g\n", rr_node[inode].R, rr_node[inode].C);
    fprintf(fp, "Cost_index: %d\n", rr_node[inode].cost_index);
}


/* Prints all the rr_indexed_data of index to file fp.   */
void
print_rr_indexed_data(FILE* fp,
                      int index)
{
    fprintf(fp, "Index: %d\n", index);
    fprintf(fp, "ortho_cost_index: %d  ",
            rr_indexed_data[index].ortho_cost_index);
    fprintf(fp, "base_cost: %g  ", rr_indexed_data[index].saved_base_cost);
    fprintf(fp, "saved_base_cost: %g\n",
            rr_indexed_data[index].saved_base_cost);
    fprintf(fp, "Seg_index: %d  ", rr_indexed_data[index].seg_index);
    fprintf(fp, "inv_length: %g\n", rr_indexed_data[index].inv_length);
    fprintf(fp, "T_linear: %g  ", rr_indexed_data[index].T_linear);
    fprintf(fp, "T_quadratic: %g  ", rr_indexed_data[index].T_quadratic);
    fprintf(fp, "C_load: %g\n", rr_indexed_data[index].C_load);
}



static void
build_unidir_rr_opins(IN int i,
                      IN int j,
                      IN grid_tile_t** grid,
                      IN int* Fc_out,
                      IN int nodes_per_chan,
                      IN segment_details_t* seg_details,
                      INOUT int** Fc_xofs,
                      INOUT int** Fc_yofs,
                      INOUT rr_node_t* rr_node,
                      INOUT boolean* rr_edge_done,
                      OUT boolean* Fc_clipped,
                      IN vector_t** * rr_node_indices)
{
    /* This routine returns a list of the opins rr_nodes on each
     * side/offset of the block. You must free the result with
     * free_matrix. */
    int ipin, iclass, ofs, chan, seg, max_len, inode;
    side_types_t side;
    rr_type_t chan_type;
    t_linked_edge* edge_list = NULL, *next;
    boolean clipped, vert, pos_dir;
    int num_edges;
    int** Fc_ofs;
    *Fc_clipped = FALSE;

    /* Only the base block of a set should use this function */
    if (grid[i][j].m_offset > 0) {
        return;
    }

    block_type_ptr type = grid[i][j].grid_type;

    /* Go through each pin and find its fanout. */
    for (ipin = 0; ipin < type->num_pins; ++ipin) {
        /* Skip global pins and ipins */
        iclass = type->pin_class[ipin];

        if (type->class_inf[iclass].type != DRIVER) {
            continue;
        }

        if (type->is_global_pin[ipin]) {
            continue;
        }

        num_edges = 0;
        edge_list = NULL;

        for (ofs = 0; ofs < type->height; ++ofs) {
            for (side = 0; side < 4; ++side) {
                /* Can't do anything if pin isn't at this location */
                if (0 == type->pinloc[ofs][side][ipin]) {
                    continue;
                }

                /* Figure out the chan seg at that side.
                 * side is the side of the logic or io block. */
                vert = ((side == TOP) || (side == BOTTOM));
                pos_dir = ((side == TOP) || (side == RIGHT));
                chan_type = (vert ? CHANX : CHANY);
                chan = (vert ? (j + ofs) : i);
                seg = (vert ? i : (j + ofs));
                max_len = (vert ? num_grid_rows : num_grid_columns);
                Fc_ofs = (vert ? Fc_xofs : Fc_yofs);

                if (FALSE == pos_dir) {
                    --chan;
                }

                /* Skip the location if there is no channel. */
                if (chan < 0) {
                    continue;
                }

                if (seg < 1) {
                    continue;
                }

                if (seg > (vert ? num_grid_rows : num_grid_columns)) {
                    continue;
                }

                if (chan > (vert ? num_grid_columns : num_grid_rows)) {
                    continue;
                }

                /* Get the list of opin to mux connections for that chan seg. */
                num_edges +=
                    get_unidir_opin_connections(chan, seg,
                                                Fc_out[type->
                                                       index],
                                                chan_type,
                                                seg_details,
                                                &edge_list,
                                                Fc_ofs,
                                                rr_edge_done,
                                                max_len,
                                                nodes_per_chan,
                                                rr_node_indices,
                                                &clipped);

                if (clipped) {
                    *Fc_clipped = TRUE;
                }
            }
        }

        /* Add the edges */
        inode = get_rr_node_index(i, j, OPIN, ipin, rr_node_indices);
        alloc_and_load_edges_and_switches(rr_node, inode, num_edges,
                                          rr_edge_done, edge_list);

        while (edge_list != NULL) {
            next = edge_list->next;
            free(edge_list);
            edge_list = next;
        }
    }
}


