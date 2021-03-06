void check_rr_graph(IN t_graph_type graph_type,
                    IN int num_types,
                    IN block_type_ptr types,
                    IN int num_grid_columns,
                    IN int num_grid_rows,
                    IN grid_tile_t** grid,
                    IN int nodes_per_chan,
                    IN int Fs,
                    IN int num_seg_types,
                    IN int num_switches,
                    IN segment_info_t* segment_inf,
                    IN int global_route_switch,
                    IN int delayless_switch,
                    IN int wire_to_ipin_switch,
                    segment_details_t* seg_details,
                    int* Fc_in,
                    int* Fc_out,
                    vector_t** * rr_node_indices,
                    int**** *opin_to_track_map,
                    int**** *ipin_to_track_map,
                    vector_t**** track_to_ipin_lookup,
                    vector_t** * switch_block_conn,
                    boolean* perturb_ipins);

void check_node(int inode,
                router_types_t route_type);
