#ifndef ROUTE_BREADTH_FIRST_H
#define ROUTE_BREADTH_FIRST_H

#include "vpr_types_parallel.h"

boolean try_breadth_first_route(router_opts_t router_opts,
                                vector_t** fb_opins_used_locally,
                                int width_fac);
#endif

