#ifndef ROUTE_DIRECTED_SEARCH_H
#define ROUTE_DIRECTED_SEARCH_H

#include "vpr_types_parallel.h"

boolean try_directed_search_route(router_opts_t router_opts,
                                  vector_t** fb_opins_used_locally,
                                  int width_fac,
                                  t_mst_edge** mst);

#endif

