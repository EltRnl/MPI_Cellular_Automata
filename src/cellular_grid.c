#include <stdlib.h>
#include "cellular_grid.h"


_Bool valid_coordinates_cell(cellular_grid CG, int x, int y){
    return x>=-1 && x<=CG->inner_width+1 && y>=-1 && y<=CG->inner_height+1;
}

bit* get_neighbors(cellular_grid CG, int x, int y){
    bit* neighbors = malloc(9);
    for(int j=0; j<3; j++){
        for(int i=0; i<3; i++){
            neighbors[j*3+i] = get_cell(CG,x-1+i,y-1+j)>0;
        }
    }
    return neighbors;
}


cellular_grid create_cell_grid(uint width, uint height, bit (* convolution) (bit *)){
    cellular_grid CG = malloc(sizeof(struct _cellular_grid));
    CG->grid = create_grid(width+2,height+2);
    CG->convolution = convolution;
    CG->width = width + 2;
    CG->height = height + 2;    
    CG->inner_width = width;
    CG->inner_height = height;
}

void delete_cell_grid(cellular_grid CG){
    delete_grid(CG->grid);
    free(CG);
}


int get_cell(cellular_grid CG, int x, int y){
    if (!valid_coordinates_cell(CG,x,y)) return -1;
    return get_bit(CG->grid,x+1,y+1);
}

int set_cell(cellular_grid CG, int x, int y, int new_value){
    if (!valid_coordinates_cell(CG,x,y)) return -1;
    return set_bit(CG->grid,x+1,y+1,new_value>0);
}

bit* get_wall(cellular_grid CG, enum side s){
    bit* values;
    switch (s){
    case North:
        values = malloc(sizeof(bit)*CG->inner_width);
        for(int x=-1; x<=CG->inner_width; x++) values[x] = get_cell(CG,x,0);
        break;
    case South:
        values = malloc(sizeof(bit)*CG->inner_width);
        for(int x=-1; x<=CG->inner_width; x++) values[x] = get_cell(CG,x,CG->inner_height-1);
        break;
    
    case West:
        values = malloc(sizeof(bit)*CG->inner_height);
        for(int y=-1; y<=CG->inner_width; y++) values[y] = get_cell(CG,0,y);
        break;
    case East:
        values = malloc(sizeof(bit)*CG->inner_height);
        for(int y=-1; y<=CG->inner_width; y++) values[y] = get_cell(CG,CG->inner_width-1,y);
        break;

    default:
        return NULL;
        break;
    }

    return values;
}

int set_wall(cellular_grid CG, enum side s, bit* values){
    switch (s){
    case North:
        for(int x=-1; x<=CG->inner_width; x++) set_cell(CG,x,-1,values[x]?1:0);
        break;
    case South:
        for(int x=-1; x<=CG->inner_width; x++) set_cell(CG,x,CG->inner_height,values[x]?1:0);
        break;
    
    case West:
        for(int y=-1; y<=CG->inner_width; y++) set_cell(CG,-1,y,values[y]?1:0);
        break;
    case East:
        for(int y=-1; y<=CG->inner_width; y++) set_cell(CG,CG->inner_width,y,values[y]?1:0);
        break;

    default:
        return -1;
        break;
    }
    return 1;
}

void next_generation(cellular_grid CG){
    grid new_generation = copy(CG->grid);
    for(int y=0; y<CG->inner_height; y++){
        for(int x=0; x<CG->inner_width; x++){
            bit* n = get_neighbors(CG,x,y);
            set_bit(new_generation,x+1,y+1,CG->convolution(n));
            free(n);
        }
    }
    delete_grid(CG->grid);
    CG->grid = new_generation;
}

void print_cell_grid(cellular_grid CG){
    printf("\e[1;1H\e[2J");
    for(int y=0; y<CG->inner_height; y++){
        for(int x=0; x<CG->inner_width; x++){
            if(get_cell(CG,x,y)>0)
                printf("\u2B1B");
            else 
                printf("\u2B1C");
        }
        printf("\n");
    }	
}

int generate_svg_points(cellular_grid CG, int generation, cell_point** points){
    int nb_points = 0;
    for (int y=0; y<CG->inner_height; y++) for (int x=0; x<CG->inner_width; x++) if (get_cell(CG,x,y) == 1) nb_points++;

    *points = (cell_point*)malloc(nb_points*sizeof(cell_point));
    cell_point* list = *points;
    int index = 0;
    for (int y=0; y<CG->inner_height; y++)
        for (int x=0; x<CG->inner_width; x++){
            if(index==nb_points) break;

            if (get_cell(CG,x,y) == 1){
                list[index].gen = generation;
                list[index].x = x;
                list[index].y = y;
                index++;
            }
        }

    return nb_points;

}