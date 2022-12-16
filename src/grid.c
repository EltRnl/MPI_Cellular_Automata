#include <stdlib.h>
#include "grid.h"


_Bool valid_coordinates(grid G, uint x, uint y){
    return x<G -> width && y<G -> height;
}

grid create_grid(uint width, uint height){
    grid G = malloc(sizeof(struct _grid));
    G -> width = width;
    G -> height = height;
    G -> size = width*height;
    G -> value = (bit *) malloc(height*width*sizeof(bit));
    return G;
}

void delete_grid(grid G){
    free(G->value);
    free(G);
}

int get_bit(grid G, uint x, uint y){
    if (!valid_coordinates(G,x,y)) return -1;
    return G -> value[y*G->width + x];
}

int set_bit(grid G, uint x, uint y, bit new_bit){
    if (!valid_coordinates(G,x,y)) return -1;
    G -> value[y*G->width + x] = new_bit;
    return 1;
}

int set_bits(grid G, bit * new_values){
    if(sizeof(new_values)<G -> size) return -1;
    for(uint i=0; i<G -> size; i++){
        G->value[i] = new_values[i];
    }
    return 1;
}

int set_bits_points(grid G, point * p, uint nb_points){
    uint status = 1;
    for(uint i=0; i<nb_points; i++){
        if(valid_coordinates(G,p[i]->x,p[i]->y))
            set_bit(G,p[i]->x,p[i]->y,1);
        else 
            status--;
    }
    return status;
}

grid copy(grid G){
    grid g = create_grid(G->width,G->height);
    for(uint i=0; i<G -> size; i++){
        g->value[i] = G->value[i];
    }
    return g;
}