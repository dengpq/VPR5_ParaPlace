#ifndef DRAW_H
#define DRAW_H

#include "vpr_types_parallel.h"

void update_screen(int priority,
                   char* msg,
                   pic_type_t pic_on_screen_val,
                   boolean crit_path_button_enabled);

void alloc_draw_structs(void);

void init_draw_coords(double clb_width);

void set_graphics_state(boolean show_graphics_val,
                        int gr_automode_val,
                        router_types_t route_type);
#endif

