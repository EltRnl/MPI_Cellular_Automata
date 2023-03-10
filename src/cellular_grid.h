#include <stdlib.h>
#include <stdio.h>
#include "grid.h"

enum side{North,East,South,West};

struct _cellular_grid{
    grid grid;
    bit (* convolution) (bit *);
    int width;
    int height;    
    int inner_width;
    int inner_height;
};

struct _cell_point{
    int gen;
    int x;
    int y;
};


typedef struct _cellular_grid * cellular_grid;
typedef struct _cell_point cell_point;

cellular_grid create_cell_grid(uint width, uint height, bit (* convolution) (bit *));

void delete_cell_grid(cellular_grid CG);

int get_cell(cellular_grid CG, int x, int y);

int set_cell(cellular_grid CG, int x, int y, int new_value);

void get_wall(cellular_grid CG, enum side s, int* values);

int set_wall(cellular_grid CG, enum side s, int* values);

void next_generation(cellular_grid CG);

void print_cell_grid(cellular_grid CG);