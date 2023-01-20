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
 * @brief Sends the walls of our local grid to our neighbors, while receiving theirs (might be ourself if one of the dimensions is 1).
 * Uses blocking communications without any optimisation, and without taking into account the case where one dimension is odd.
 * 
 * @param CG Our local Cellular Grid
 * @param comm The communication schema
 */
void transmit_walls_V1_0(cellular_grid CG, struct comm_schema comm){
    // Bit arrays used to send and receiving the walls
    bit HWall_recv[CG->height * 2];    
    bit* HWall_send;    
    bit VWall_recv[CG->width * 2];    
    bit* VWall_send;    

    MPI_Status status;
    MPI_Status * used_status = &status; // MPI_STATUS_IGNORE

    // Phase 1 : Horizontal Transfer
    switch (comm.x%2){
    case 0:
        // Sub-phase A : sending and receiving East
        HWall_send = get_wall(CG,East);
        MPI_Send( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 0 , MPI_COMM_WORLD);
        MPI_Recv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 0 , MPI_COMM_WORLD , used_status);
        set_wall(CG,East,HWall_recv);

        // printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        // Sub-phase B : sending and receiving West
        MPI_Recv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 1 , MPI_COMM_WORLD , used_status);
        set_wall(CG,West,HWall_recv);
        HWall_send = get_wall(CG,West);
        MPI_Send( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 1 , MPI_COMM_WORLD);

        // printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        break;

    case 1:
        // Sub-phase A : sending and receiving West
        MPI_Recv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 0 , MPI_COMM_WORLD , used_status);
        set_wall(CG,West,HWall_recv);
        HWall_send = get_wall(CG,West);
        MPI_Send( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 0 , MPI_COMM_WORLD);

        // printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        // Sub-phase B : sending and receiving East
        HWall_send = get_wall(CG,East);
        MPI_Send( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 1 , MPI_COMM_WORLD);
        MPI_Recv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 1 , MPI_COMM_WORLD , used_status);
        set_wall(CG,East,HWall_recv);

        // printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        break;

    default:
        break;
    }

    // Phase 2 : Vertical Transfer
    switch (comm.y%2){
    case 0:
        // Sub-phase A : sending and receiving North
        VWall_send = get_wall(CG,North);
        MPI_Send( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 0 , MPI_COMM_WORLD);
        MPI_Recv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 0 , MPI_COMM_WORLD , used_status);
        set_wall(CG,North,VWall_recv);

        // printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        // Sub-phase B : sending and receiving South
        MPI_Recv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 1 , MPI_COMM_WORLD , used_status);
        set_wall(CG,South,VWall_recv);
        VWall_send = get_wall(CG,South);
        MPI_Send( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 1 , MPI_COMM_WORLD);

        // printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        break;

    case 1:
        // Sub-phase A : sending and receiving South
        MPI_Recv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 0 , MPI_COMM_WORLD , used_status);
        set_wall(CG,South,VWall_recv);
        VWall_send = get_wall(CG,South);
        MPI_Send( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 0 , MPI_COMM_WORLD);

        // printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        // Sub-phase B : sending and receiving North
        VWall_send = get_wall(CG,North);
        MPI_Send( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 1 , MPI_COMM_WORLD);
        MPI_Recv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 1 , MPI_COMM_WORLD , used_status);
        set_wall(CG,North,VWall_recv);

        // printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);

        break;

    default:
        break;
    }
}

/**
 * @brief Sends the walls of our local grid to our neighbors, while receiving theirs (might be ourself if one of the dimensions is 1).
 * Uses non-blocking communications without any optimisation, and without taking into account the case where one dimension is odd.
 * 
 * @param CG Our local Cellular Grid
 * @param comm The communication schema
 */
void transmit_walls_V1_1(cellular_grid CG, struct comm_schema comm){
    // Bit arrays used to send and receiving the walls
    bit HWall_recv[CG->height * 2];    
    bit* HWall_send;    
    bit VWall_recv[CG->width * 2];    
    bit* VWall_send;    

    MPI_Status status;
    MPI_Status * used_status = &status; // MPI_STATUS_IGNORE
    MPI_Request req_send;
    MPI_Request req_recv;

    // Phase 1 : Horizontal Transfer
    switch (comm.x%2){
    case 0:
        // Sub-phase A : sending and receiving East
        HWall_send = get_wall(CG,East);
        MPI_Isend( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 0 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 0 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,East,HWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        // Sub-phase B : sending and receiving West
        HWall_send = get_wall(CG,West);
        MPI_Isend( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 1 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 1 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,West,HWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        break;

    case 1:
        // Sub-phase A : sending and receiving West
        HWall_send = get_wall(CG,West);
        MPI_Isend( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 0 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x-1,comm.y) , 0 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,West,HWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        // Sub-phase B : sending and receiving East
        HWall_send = get_wall(CG,East);
        MPI_Isend( HWall_send , CG->height , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 1 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( HWall_recv , CG->height * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x+1,comm.y) , 1 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,East,HWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Horizontal : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        break;

    default:
        break;
    }

    // Phase 2 : Vertical Transfer
    switch (comm.y%2){
    case 0:
        // Sub-phase A : sending and receiving North
        VWall_send = get_wall(CG,North);
        MPI_Isend( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 0 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 0 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,North,VWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        // Sub-phase B : sending and receiving South
        VWall_send = get_wall(CG,South);
        MPI_Isend( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 1 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 1 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,South,VWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        break;

    case 1:
        // Sub-phase A : sending and receiving South
        VWall_send = get_wall(CG,South);
        MPI_Isend( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 0 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y+1) , 0 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,South,VWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        // Sub-phase B : sending and receiving North
        VWall_send = get_wall(CG,North);
        MPI_Isend( VWall_send , CG->width , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 1 , MPI_COMM_WORLD , &req_send);
        MPI_Irecv( VWall_recv , CG->width * 2 , MPI_C_BOOL , position_to_rank(comm.width,comm.height,comm.x,comm.y-1) , 1 , MPI_COMM_WORLD , &req_recv);
        MPI_Wait( &req_recv , used_status);
        set_wall(CG,North,VWall_recv);
        MPI_Wait( &req_send , NULL);

        #ifdef V2
        printf("Vertical : Node %d received something from node %d.\n",comm.rank,status.MPI_SOURCE);fflush(stdout);
        #endif

        break;

    default:
        break;
    }
}


/**
 * @brief Gather all the data to 1 node for later display
 * 
 * @param CG Our local Cellular Grid
 * @param comm The communication schema
 * @param master The rank of the node choosen to gather everything
 */
void gather_to_one_and_print(cellular_grid CG, struct comm_schema comm){
    if(comm.rank==comm.master){
        //print_cell_grid(CG);
    }
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

    comm.master = 0; // TODO : Maybe have a master node be elected (or determine who can display/save the result on user's machine)

    // Automata grid creation

    /* We divide the whole cellular_grid into equal pieces, rounding up to the closest integer (with 0.5 being 0) */
    int local_width = rounded_division(WIDTH,comm.width);
    int local_height = rounded_division(HEIGHT,comm.height);
    /* If it's a border on the right or bottom, we give it the rest of the width of height, to assure that we comply with the grid size we want. */
    if(comm.x==comm.width-1) local_width = WIDTH - local_width * (comm.width - 1);
    if(comm.y==comm.height-1) local_height = HEIGHT - local_height * (comm.height - 1);

    cellular_grid CG = create_cell_grid(local_width,local_height,conway);

    #ifdef V1
    if(comm.rank==comm.master){
        printf("\nComm : %d x %d\nNode : %d x %d\n",comm.width,comm.height,local_width,local_height); fflush(stdout);
    }
    #endif

    // Automata grid values initialization at random
    time_t t;
    srand((unsigned) time(&t));

    for(int i=0; i<local_width*local_height/2; i++){
        set_cell(CG,rand()%local_width,rand()%local_height,1);
    }

    clock_t start, end;        

    // Main loop
    for(int i=0; i<ITERATIONS; i++){
        #ifdef V1
        if(comm.rank==comm.master) start = clock();
        #endif

        // Generation computation
        transmit_walls_V1_0(CG,comm);
        next_generation(CG);

        // Display/Saving
        gather_to_one_and_print(CG,comm); 
        
        #ifdef V1
        if(comm.rank==comm.master) {
            end = clock();
            printf("Generation %d done in %lfs.\n",i,(double)(end-start)/CLOCKS_PER_SEC);
        }
        #endif

        // // Time skip
        // usleep(TIME_INTERVAL_U);
    }

    MPI_Finalize();
    return 0;
}