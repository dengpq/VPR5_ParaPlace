#ifndef RR_GRAPH_AREA_H
#define RR_GRAPH_AREA_H

#include "vpr_types_parallel.h"

void count_routing_transistors(directionality_t directionality,
                               int num_switch,
                               segment_info_t* segment_inf,
                               double R_minW_nmos,
                               double R_minW_pmos);
#endif

