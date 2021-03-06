#include <assert.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "vpr_utils.h"

/* This module contains subroutines that are used in several unrelated parts *
 * of VPR.  They are VPR-specific utility routines.                          */


/******************** Subroutine definitions ********************************/
/* Points the grid structure back to the blocks list */
void sync_grid_to_blocks(IN int num_blocks,
                         IN const block_t block_list[],
                         IN int num_grid_columns,
                         IN int num_grid_rows,
                         INOUT grid_tile_t** grid)
{
    int i, j, k;

    /* Reset usage and allocate blocks list if needed */
    for (j = 0; j <= (num_grid_rows + 1); ++j) {
        for (i = 0; i <= (num_grid_columns + 1); ++i) {
            grid[i][j].m_usage = 0;

            if (grid[i][j].grid_type != NULL) {
                /* If already allocated, leave it since size doesn't change */
                if (NULL == grid[i][j].in_blocks) {
                    grid[i][j].in_blocks = (int*)my_malloc(sizeof(int) *
                                                        grid[i][j].grid_type->capacity);
                    /* Set them as unconnected */
                    for (k = 0; k < grid[i][j].grid_type->capacity; ++k) {
                        grid[i][j].in_blocks[k] = OPEN;
                    }
                }
            }
        }
    }

    /* Go through each block */
    for (i = 0; i < num_blocks; ++i) {
        /* Check range of block coords */
        if (blocks[i].x < 0 || blocks[i].x > (num_grid_columns + 1)
              || blocks[i].y < 0
              || (blocks[i].y + blocks[i].block_type->height - 1) > (num_grid_rows + 1)
              || blocks[i].z < 0 || blocks[i].z > (blocks[i].block_type->capacity)) {
            printf(ERRTAG
                   "Block %d is at invalid location (%d, %d, %d)\n",
                   i, blocks[i].x, blocks[i].y, blocks[i].z);
            exit(1);
        }

        /* Check types match */
        if (blocks[i].block_type != grid[blocks[i].x][blocks[i].y].grid_type) {
            printf(ERRTAG "A block is in a grid location "
                   "(%d x %d) with a conflicting type.\n", blocks[i].x,
                   blocks[i].y);
            exit(1);
        }

        /* Check already in use */
        if (OPEN != grid[blocks[i].x][blocks[i].y].in_blocks[blocks[i].z]) {
            printf(ERRTAG
                   "Location (%d, %d, %d) is used more than once\n",
                   blocks[i].x, blocks[i].y, blocks[i].z);
            exit(1);
        }

        if (grid[blocks[i].x][blocks[i].y].m_offset != 0) {
            printf(ERRTAG
                   "Large block not aligned in placment for block %d at (%d, %d, %d)",
                   i, blocks[i].x, blocks[i].y, blocks[i].z);
            exit(1);
        }

        /* Set the block */
        for (j = 0; j < blocks[i].block_type->height; j++) {
            grid[blocks[i].x][blocks[i].y + j].in_blocks[blocks[i].z] = i;
            ++(grid[blocks[i].x][blocks[i].y + j].m_usage);
            assert(grid[blocks[i].x][blocks[i].y + j].m_offset == j);
        }
    }
}

/* This function updates the nets list to point back to blocks list */
void sync_nets_to_blocks(IN int num_blocks,
                    IN const block_t block_list[],
                    IN int num_nets,
                    INOUT net_t net_list[])
{
    int i, j, k, l;
    block_type_ptr cur_type;

    /* Count the number of sinks for each net */
    for (j = 0; j < num_blocks; ++j) {
        cur_type = block_list[j].block_type;

        for (k = 0; k < cur_type->num_type_pins; ++k) {
            i = block_list[j].nets[k];

            if ( i >= 0 && RECEIVER == cur_type->class_inf[cur_type->pin_class[k]].type) {
                ++net_list[i].num_net_pins;
            }
        }
    }

    /* Alloc and load block lists of nets */
    int* next_sink = (int*)my_malloc( num_nets * sizeof(int) );

    for (i = 0; i < num_nets; ++i) {
        /* The list should be num_sinks + 1 driver. Re-alloc if already allocated. */
        if (net_list[i].node_blocks != NULL) {
            free(net_list[i].node_blocks);
        }

        net_list[i].node_blocks =
            (int*)my_malloc(sizeof(int) * (net_list[i].num_net_pins + 1));

        if (net_list[i].node_block_pins != NULL) {
            free(net_list[i].node_block_pins);
        }

        net_list[i].node_block_pins =
            (int*)my_malloc(sizeof(int) * (net_list[i].num_net_pins + 1));

        next_sink[i] = 1;       /* First sink goes at position 1, since 0 is for driver */
    }

    for (j = 0; j < num_blocks; ++j) {
        cur_type = block_list[j].block_type;

        for (k = 0; k < cur_type->num_type_pins; ++k) {
            i = block_list[j].nets[k];

            if (i >= 0 && RECEIVER == cur_type->class_inf[cur_type->pin_class[k]].type) {
                l = next_sink[i];
                net_list[i].node_blocks[l] = j;
                net_list[i].node_block_pins[l] = k;
                ++next_sink[i];
            } else if (i >= 0) {
                assert(DRIVER == cur_type->class_inf[cur_type->pin_class[k]].type);
                net_list[i].node_blocks[0] = j;
                net_list[i].node_block_pins[0] = k;
            }
        }
    }

    free( next_sink );

}

boolean
is_opin(int ipin,
        block_type_ptr type)
{

    /* Returns TRUE if this clb pin is an output, FALSE otherwise. */
    int iclass = type->pin_class[ipin];
    if (type->class_inf[iclass].type == DRIVER) {
        return (TRUE);
    } else {
        return (FALSE);
    }
}

void
get_class_range_for_block(IN int iblk,
                          OUT int* class_low,
                          OUT int* class_high)
{
    /* Assumes that the placement has been done so each block has a set of pins allocated to it */
    block_type_ptr type = blocks[iblk].block_type;
    assert(type->num_class % type->capacity == 0);
    *class_low = blocks[iblk].z * (type->num_class / type->capacity);
    *class_high =
        (blocks[iblk].z + 1) * (type->num_class / type->capacity) - 1;
}


void
load_one_fb_fanout_count(subblock_t* subblock_inf,
                         int num_subblocks,
                         int* num_uses_of_fb_ipin,
                         int** num_uses_of_sblk_opin,
                         int iblk)
{

    /* Loads the fanout counts for one block (iblk).  */
    block_type_ptr type = blocks[iblk].block_type;
    int isub, ipin, conn_pin, opin;
    int internal_sub, internal_pin;

    /* Reset ipin counts */
    for (ipin = 0; ipin < type->num_type_pins; ipin++) {
        num_uses_of_fb_ipin[ipin] = 0;
    }

    /* First pass, reset fanout counts */
    for (isub = 0; isub < num_subblocks; isub++) {
        for (opin = 0; opin < type->max_subblock_outputs; opin++) {
            num_uses_of_sblk_opin[isub][opin] = 0;
        }
    }

    for (isub = 0; isub < num_subblocks; isub++) {
        /* Is the subblock output connected to a FB opin that actually goes *
         * somewhere?  Necessary to check that the FB opin connects to      *
         * something because some logic blocks result in netlists where      *
         * subblock outputs being automatically hooked to a FB opin under   *
         * all conditions.                                                   */
        for (opin = 0; opin < type->max_subblock_outputs; opin++) {
            conn_pin = subblock_inf[isub].outputs[opin];

            if (conn_pin != OPEN) {
                if (blocks[iblk].nets[conn_pin] != OPEN) {
                    /* FB output is used */
                    num_uses_of_sblk_opin[isub][opin]++;
                }
            }
        }

        for (ipin = 0; ipin < type->max_subblock_inputs; ipin++) {
            conn_pin = subblock_inf[isub].inputs[ipin];
            if (conn_pin != OPEN) {
                if (conn_pin < type->num_type_pins) {
                    /* Driven by FB ipin */
                    num_uses_of_fb_ipin[conn_pin]++;
                } else {
                    /* Driven by sblk output in same fb */
                    internal_sub = (conn_pin - type->num_type_pins) /
                                     type->max_subblock_outputs;
                    internal_pin = (conn_pin - type->num_type_pins) %
                                     type->max_subblock_outputs;
                    num_uses_of_sblk_opin[internal_sub][internal_pin]++;
                }
            }
        }       /* End for each sblk ipin */

        conn_pin = subblock_inf[isub].clock;    /* Now do clock pin */
        if (conn_pin != OPEN) {
            if (conn_pin < type->num_type_pins) {
                /* Driven by FB ipin */
                num_uses_of_fb_ipin[conn_pin]++;
            } else {
                /* Driven by sblk output in same clb */
                internal_sub =
                  (conn_pin - type->num_type_pins) / type->max_subblock_outputs;
                internal_pin =
                  (conn_pin - type->num_type_pins) % type->max_subblock_outputs;
                num_uses_of_sblk_opin[internal_sub][internal_pin]++;
            }
        }
    }  /* End for each subblock */
}

