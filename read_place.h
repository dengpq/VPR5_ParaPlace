#ifndef READ_PLACE_H
#define READ_PLACE_H

#include "util.h"

void read_place(IN const char* place_file,
                IN const char* arch_file,
                IN const char* net_file,
                IN int num_grid_columns,
                IN int num_grid_rows,
                IN int num_blocks,
                INOUT block_t block_list[]);

void print_place(IN char* place_file,
                 IN char* net_file,
                 IN char* arch_file);

void read_user_pad_loc(IN char* pad_loc_file);

#endif

