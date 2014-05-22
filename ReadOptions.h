#ifndef READOPTIONS_H
#define READOPTIONS_H

#include "vpr_types_parallel.h"
#include "OptionTokens.h"

typedef struct s_options {
    char* ArchFile;
    char* NetFile;
    char* PlaceFile;
    char* RouteFile;

    /* General options */
    int GraphPause;
    double constant_net_delay;
    boolean TimingAnalysis;
    char* OutFilePrefix;

    /* Placement options */
    place_algorithm_t PlaceAlgorithm;
    double PlaceInitT;
    double PlaceExitT;
    double PlaceAlphaT;
    double PlaceInnerNum;
    int Seed;
    double place_cost_exp;
    place_cong_types_t PlaceCostType;
    int PlaceChanWidth;
    int PlaceNonlinearRegions;
    char* PinFile;
    boolean ShowPlaceTiming;
    int block_dist;

    /* Timing-driven placement options only */
    double PlaceTimingTradeoff;
    int RecomputeCritIter;
    int inner_loop_recompute_divider;
    double place_exp_first;
    double place_exp_last;

    /* Router Options */
    int max_router_iterations;
    int bb_factor;
    double initial_pres_fac;
    double pres_fac_mult;
    double acc_fac;
    double first_iter_pres_fac;
    double bend_cost;
    router_types_t RouteType;
    int RouteChanWidth;
    router_algorithm_t RouterAlgorithm;
    router_base_cost_t base_cost_type;

    /* Timing-driven router options only */
    double astar_fac;
    double criticality_exp;
    double max_criticality;

    int Count[OT_BASE_UNKNOWN];
} t_options;

void ReadOptions(IN int argc,
                 IN char** argv,
                 OUT t_options* Options);

#endif

