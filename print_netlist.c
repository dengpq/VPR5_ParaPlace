#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "util.h"
#include "vpr_types_parallel.h"
#include "globals.h"
#include "print_netlist.h"


/******************** Subroutines local to this module ***********************/

static void print_pinnum(FILE* fp,
                         int pinnum);


/********************* Subroutine definitions ********************************/
void print_netlist(char* foutput,
                   char* net_file,
                   subblock_data_t subblock_data)
{
    /* Prints out the netlist related data structures into the file    *
     * fname.                                                          */
    int i, j, ipin, max_pin;
    int num_global_nets;
    int num_p_inputs, num_p_outputs;
    FILE* fp;
    subblock_t** subblock_inf;
    int* num_subblocks_per_block;
    num_global_nets = 0;
    num_p_inputs = 0;
    num_p_outputs = 0;

    /* Count number of global nets */
    for (i = 0; i < num_nets; i++) {
        if (!net[i].is_global) {
            num_global_nets++;
        }
    }

    /* Count I/O input and output pads */
    for (i = 0; i < num_blocks; i++) {
        if (blocks[i].block_type == IO_TYPE) {
            for (j = 0; j < IO_TYPE->num_type_pins; j++) {
                if (blocks[i].nets[j] != OPEN) {
                    if (IO_TYPE->class_inf[IO_TYPE->pin_class[j]].type == DRIVER) {
                        num_p_inputs++;
                    } else {
                        assert(IO_TYPE->class_inf[IO_TYPE->pin_class[j]].type ==
                                RECEIVER);
                        num_p_outputs++;
                    }
                }
            }
        }
    }

    fp = my_fopen(foutput, "w");
    fprintf(fp, "Input netlist file: %s\n", net_file);
    fprintf(fp, "num_p_inputs: %d, num_p_outputs: %d, num_clbs: %d\n",
            num_p_inputs, num_p_outputs, num_blocks);
    fprintf(fp, "num_blocks: %d, num_nets: %d, num_globals: %d\n",
            num_blocks, num_nets, num_global_nets);
    fprintf(fp, "\nNet\tName\t\t#Pins\tDriver\t\tRecvs. (block, pin)\n");

    for (i = 0; i < num_nets; ++i) {
        fprintf(fp, "\n%d\t%s\t", i, net[i].name);

        if (strlen(net[i].name) < 8) {
            fprintf(fp, "\t");    /* Name field is 16 chars wide */
        }

        const int knum_net_pins = net[i].num_net_pins;
        fprintf(fp, "%d", knum_net_pins + 1);

        for (j = 0; j <= knum_net_pins; ++j)
            fprintf(fp, "\t(%4d,%4d)", net[i].node_blocks[j],
                    net[i].node_block_pins[j]);
    }

    fprintf(fp, "\nBlock\tName\t\tType\tPin Connections\n\n");

    for (i = 0; i < num_blocks; i++) {
        fprintf(fp, "\n%d\t%s\t", i, blocks[i].name);

        if (strlen(blocks[i].name) < 8) {
            fprintf(fp, "\t");    /* Name field is 16 chars wide */
        }

        fprintf(fp, "%s", blocks[i].block_type->name);
        max_pin = blocks[i].block_type->num_type_pins;

        for (j = 0; j < max_pin; j++) {
            print_pinnum(fp, blocks[i].nets[j]);
        }
    }

    fprintf(fp, "\n");
    /* Now print out subblock info. */
    subblock_inf = subblock_data.subblock_inf;
    num_subblocks_per_block = subblock_data.num_subblocks_per_block;
    fprintf(fp, "\n\nSubblock List:\n\n");

    for (i = 0; i < num_blocks; i++) {
        fprintf(fp, "\nBlock: %d (%s)\tNum_subblocks: %d\n", i,
                blocks[i].name, num_subblocks_per_block[i]);
        /* Print header. */
        fprintf(fp, "Index\tName\t\tInputs");

        for (j = 0; j < blocks[i].block_type->max_subblock_inputs; j++) {
            fprintf(fp, "\t");
        }

        fprintf(fp, "Outputs");

        for (j = 0; j < blocks[i].block_type->max_subblock_outputs; j++) {
            fprintf(fp, "\t");
        }

        fprintf(fp, "Clock\n");

        /* Print subblock info for block i. */

        for (j = 0; j < num_subblocks_per_block[i]; j++) {
            fprintf(fp, "%d\t%s", j, subblock_inf[i][j].name);

            if (strlen(subblock_inf[i][j].name) < 8) {
                fprintf(fp, "\t");    /* Name field is 16 characters */
            }

            for (ipin = 0; ipin < blocks[i].block_type->max_subblock_inputs;
                    ipin++) {
                print_pinnum(fp, subblock_inf[i][j].inputs[ipin]);
            }

            for (ipin = 0; ipin < blocks[i].block_type->max_subblock_outputs;
                    ipin++) {
                print_pinnum(fp, subblock_inf[i][j].outputs[ipin]);
            }

            print_pinnum(fp, subblock_inf[i][j].clock);
            fprintf(fp, "\n");
        }
    }

    fclose(fp);
}


static void
print_pinnum(FILE* fp,
             int pinnum)
{
    /* This routine prints out either OPEN or the pin number, to file fp. */
    if (pinnum == OPEN) {
        fprintf(fp, "\tOPEN");
    } else {
        fprintf(fp, "\t%d", pinnum);
    }
}
