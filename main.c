#include "src/cellular_grid.h"
#include "settings.h"
#include <time.h>
#include <unistd.h>

bit conway (bit* neighbors){
    int count = 0;
    bit self = neighbors[4];
    for(int y=0; y<3; y++){
        for(int x=0; x<3; x++){
            if(y==1 && x==1) continue;
            if(neighbors[y*3+x]) count++;
        }
    }
    return self?(count==2 || count==3):(count==3);
}

void main(){
    cellular_grid CG = create_cell_grid(WIDTH,HEIGHT,conway);

    time_t t;
    srand((unsigned) time(&t));

    for(int i=0; i<WIDTH*HEIGHT/2; i++){
        set_cell(CG,rand()%WIDTH,rand()%HEIGHT,1);
    }

    
    for(int i=0; i<ITERATIONS; i++){
        printf("\e[1;1H\e[2J");
        print_cell_grid(CG);
        next_generation(CG);
        usleep(TIME_INTERVAL_U);
    }

}