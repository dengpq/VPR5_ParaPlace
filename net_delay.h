#ifndef NET_DELAY_H
#define NET_DELAY_H

#include "util.h"

double** alloc_net_delay(linked_vptr_t** chunk_list_head_ptr);

void free_net_delay(double** net_delay,
                    linked_vptr_t** chunk_list_head_ptr);

void load_net_delay_from_routing(double** net_delay);

void load_constant_net_delay(double** net_delay,
                             double delay_value);

void print_net_delay(double** net_delay,
                     char* fname);

#endif
