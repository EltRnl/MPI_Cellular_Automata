#include "communication_utils.h"
#include "cellular_grid.h"
#include "settings.h"

#include <time.h>
#include <unistd.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>

/***************************** Math and Coordinate functions *****************************/

int find_factor(int n){
    int max_factor = (int)sqrt((double)n);
    int k;
    for(k=max_factor; k>1; k--)
        if(n%k==0) 
            return k;
    return 1;
}

int rounded_division(int dividend, int divisor){
    int q = dividend/divisor;
    int r = (dividend%divisor > divisor/2)?1:0;
    return q+r;
}

int position_to_rank(int width, int height, int x, int y){
    y = (y+height)%height;
    x = (x+width)%width;

    return x+y*width;
}

/***************************** Convolution functions *****************************/

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

/***************************** Communication functions *****************************/

/**
 * @brief Sends the walls of our local grid to our neighbors, while receiving theirs (might be ourself if one of the dimensions is 1)
 * 
 * @param CG Our local Cellular Grid
 * @param comm The communication schema
 */
void transmit_walls(cellular_grid CG, struct comm_schema comm){

}

/**
 * @brief Gather all the data to 1 node for later display
 * 
 * @param CG Our local Cellular Grid
 * @param comm The communication schema
 * @param master The rank of the node choosen to gather everything
 */
void gather_to_one_and_print(cellular_grid CG, struct comm_schema comm, int master){
    if(comm.rank==0) print_cell_grid(CG);
}

/***************************** Main loop function *****************************/

int automata_loop(int argc, char** argv){
    // MPI Initialization 
    MPI_Init(&argc,&argv);
    struct comm_schema comm;

    MPI_Comm_size(MPI_COMM_WORLD,&comm.size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm.rank);

    // Communication schema creation (virtual grid of automata cells)
    comm.height = find_factor(comm.size);
    comm.width = comm.size/comm.height;

    comm.y = comm.rank/comm.width;
    comm.x = comm.rank%comm.width;

    // Automata grid creation

    /* We divide the whole cellular_grid into equal pieces, rounding up to the closest integer (with 0.5 being 0) */
    int local_width = rounded_division(WIDTH,comm.width);
    int local_height = rounded_division(HEIGHT,comm.height);
    /* If it's a border on the right or bottom, we give it the rest of the width of height, to assure that we comply with the grid size we want. */
    if(comm.x==comm.width-1) local_width = WIDTH - local_width * (comm.width - 1);
    if(comm.y==comm.height-1) local_height = HEIGHT - local_height * (comm.height - 1);

    cellular_grid CG = create_cell_grid(local_width,local_height,conway);

    // Automata grid values initialization at random
    time_t t;
    srand((unsigned) time(&t));

    for(int i=0; i<local_width*local_height/2; i++){
        set_cell(CG,rand()%local_width,rand()%local_height,1);
    }

    // Main loop
    for(int i=0; i<ITERATIONS; i++){
        // Generation computation
        transmit_walls(CG,comm);
        next_generation(CG);

        // Display/Saving
        gather_to_one_and_print(CG,comm,0); // TODO : Maybe have a master node be elected (or determine who can display/save the result on user's machine)
        if(comm.rank==0){
            printf("\nComm : %d x %d\nNode : %d x %d",comm.width,comm.height,local_width,local_height); fflush(stdout);
        }
        
        // Time skip
        usleep(TIME_INTERVAL_U);
    }

    MPI_Finalize();
    return 0;
}