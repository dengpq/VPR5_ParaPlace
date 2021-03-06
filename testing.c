#include "util.h"
#include "vpr_types_parallel.h"
#include "globals_declare.h"
#include "read_arch.h"
#include "rr_graph.h"
#include "draw.h"
#include "graphics.h"
#include <assert.h>

int
main()
{
    char msg[BUFSIZE] = "This is a test.";
    detail_routing_arch_t det_routing_arch;
    segment_info_t* segment_inf;
    timing_info_t timing_inf;
    subblock_data_t subblock_data;
    chan_width_distr_t chan_width_dist;
    int nodes_per_chan;
    int i, j;
    read_arch("test.arch", DETAILED, &det_routing_arch, &segment_inf,
              &timing_inf, &subblock_data, &chan_width_dist);
    print_arch("test.arch", DETAILED, det_routing_arch, segment_inf,
               timing_inf, subblock_data, chan_width_dist);
    num_clbs = 64;
    init_arch(1., FALSE);
    printf("num_grid_columns = %d; num_grid_rows = %d\n", num_grid_columns, num_grid_rows);
    printf("Setting Nodes Per Channel ...\n");
    /* nodes_per_chan = 32; */
    nodes_per_chan = 8;

    for (i = 0; i < num_grid_columns + 1; i++) {
        chan_width_x[i] = nodes_per_chan;
    }

    for (i = 0; i < num_grid_rows + 1; i++) {
        chan_width_y[i] = nodes_per_chan;
    }

    printf("Building rr_graph ...\n");
    build_rr_graph(DETAILED, det_routing_arch, segment_inf, timing_inf,
                   INTRINSIC_DELAY);
    printf("Dumpping rr_graph ...\n");
    dump_rr_graph("rr_graph.echo");
    printf("Done.\n");
    printf("num_nets = %d\n", num_nets);
    printf("s_net = %d\n", (int)net);
    num_nets = 0;
    num_blocks = 0;
    clb = my_malloc(sizeof(struct s_clb*) * (num_grid_columns + 2));

    for (i = 0; i < num_grid_columns + 2; i++) {
        clb[i] = my_malloc(sizeof(struct s_clb) * (num_grid_rows + 2));
    }

    for (i = 0; i < num_grid_columns + 2; i++) {
        for (j = 0; j < num_grid_rows + 2; j++) {
            clb[i][j].type = CLB;
            clb[i][j].occ = 0;
            clb[i][j].u.block = 0;
        }
    }

    for (i = 0; i < num_grid_columns + 2; i++) {
        clb[i][0].type = IO;
        clb[i][num_grid_rows + 1].type = IO;
        clb[i][0].u.io_blocks = my_malloc(sizeof(int) * io_rat);
        clb[i][num_grid_rows + 1].u.io_blocks = my_malloc(sizeof(int) * io_rat);
    }

    for (j = 0; j < num_grid_rows + 2; j++) {
        clb[0][j].type = IO;
        clb[num_grid_columns + 1][j].type = IO;
        clb[0][j].u.io_blocks = my_malloc(sizeof(int) * io_rat);
        clb[num_grid_columns + 1][j].u.io_blocks = my_malloc(sizeof(int) * io_rat);
    }

    clb[0][0].type = ILLEGAL;
    clb[0][num_grid_rows + 1].type = ILLEGAL;
    clb[num_grid_columns + 1][0].type = ILLEGAL;
    clb[num_grid_columns + 1][num_grid_rows + 1].type = ILLEGAL;
    set_graphics_state(TRUE, 0, DETAILED);
    init_graphics("testing drawing capabilities");
    alloc_draw_structs();
    init_draw_coords(pins_per_clb);
    printf("num_rr_nodes = %d\n", num_rr_nodes);
    update_screen(MAJOR, msg, ROUTING, FALSE);

    while (1);

    close_graphics();
    return 0;
}
