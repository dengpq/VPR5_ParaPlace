#ifndef XML_ARCH_H
#define XML_ARCH_H

#include "vpr_types_parallel.h"

/* type definitions */
struct s_clb_grid {
    boolean IsAuto;
    double Aspect;
    int W;
    int H;
};

typedef struct s_arch t_arch;

struct s_arch {
    chan_width_distr_t Chans;
    int N;          /* Cluster size */
    int K;         /* LUT size */
    switch_block_t SBType;
    double R_minW_nmos;
    double R_minW_pmos;
    int Fs;
    double C_ipin_cblock;
    double T_ipin_cblock;
    double grid_logic_tile_area;
    double ipin_mux_trans_size;
    struct s_clb_grid clb_grid;
    segment_info_t* Segments;
    int num_segments;
    switch_info_t* Switches;
    int num_switches;
};


/* function declarations */
void XmlReadArch(IN const char* ArchFile,
                 IN boolean timing_enabled,
                 OUT struct s_arch* arch,
                 OUT type_descriptor_t** Types,
                 OUT int* NumTypes);
void EchoArch(IN const char* EchoFile,
              IN const type_descriptor_t* Types,
              IN int NumTypes);

#endif

