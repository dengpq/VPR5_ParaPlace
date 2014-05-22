#ifndef VPR_UTILS_H
#define VPR_UTILS_H

#include "vpr_types_parallel.h"

boolean is_opin(int ipin,
                block_type_ptr type);

void get_class_range_for_block(IN int iblk,
                               OUT int* class_low,
                               OUT int* class_high);

void load_one_fb_fanout_count(subblock_t* subblock_inf,
                              int num_subblocks,
                              int* num_uses_of_fb_ipin,
                              int** num_uses_of_sblk_opin,
                              int iblk);

void sync_nets_to_blocks(IN int num_blocks,
                         IN const block_t block_list[],
                         IN int num_nets,
                         INOUT net_t net_list[]);

void sync_grid_to_blocks(IN int num_blocks,
                         IN const block_t block_list[],
                         IN int num_grid_columns,
                         IN int num_grid_rows,
                         INOUT grid_tile_t** grid);
#endif

