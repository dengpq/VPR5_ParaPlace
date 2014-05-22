#ifndef SETUPVPR_H
#define SETUPVPR_H

#include "util.h"

boolean IsTimingEnabled(IN t_options Options);

void SetupVPR(IN t_options Options,
              IN boolean TimingEnabled,
              OUT t_arch* Arch,
              OUT operation_types_t* Operation,
              OUT placer_opts_t* PlacerOpts,
              OUT annealing_sched_t* AnnealSched,
              OUT router_opts_t* RouterOpts,
              OUT detail_routing_arch_t* RoutingArch,
              OUT segment_info_t** Segments,
              OUT timing_info_t* Timing,
              OUT subblock_data_t* Subblocks,
              OUT boolean* ShowGraphics,
              OUT int* GraphPause);

void CheckSetup(IN operation_types_t Operation,
                IN placer_opts_t PlacerOpts,
                IN annealing_sched_t AnnealSched,
                IN router_opts_t RouterOpts,
                IN detail_routing_arch_t RoutingArch,
                IN segment_info_t* Segments,
                IN timing_info_t Timing,
                IN subblock_data_t Subblocks,
                IN chan_width_distr_t Chans);

void CheckArch(IN t_arch Arch,
               IN boolean TimingEnabled);

void CheckOptions(IN t_options Options,
                  IN boolean TimingEnabled);

void ShowSetup(IN t_options Options,
               IN t_arch Arch,
               IN boolean TimingEnabled,
               IN operation_types_t Operation,
               IN placer_opts_t PlacerOpts,
               IN annealing_sched_t AnnealSched,
               IN router_opts_t RouterOpts,
               IN detail_routing_arch_t RoutingArch,
               IN segment_info_t* Segments,
               IN timing_info_t Timing,
               IN subblock_data_t Subblocks);

#endif

