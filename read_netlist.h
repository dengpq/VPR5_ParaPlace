#ifndef READ_NETLIST_H
#define READ_NETLIST_H

void read_netlist(IN const char* net_file,
                  IN int num_types,
                  IN const struct s_type_descriptor block_types[],
                  IN block_type_ptr IO_type,
                  IN int io_ipin,
                  IN int io_opin,
                  OUT subblock_data_t* subblock_data_ptr,
                  OUT int* num_blocks,
                  OUT block_t* block_list[],
                  OUT int* num_nets,
                  OUT net_t* net_list[]);
#endif

