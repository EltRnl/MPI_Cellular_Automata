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

/***************************** SVG saving function *****************************/

int generate_svg_points(cellular_grid CG, int generation, cell_point** points, struct comm_schema comm){
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
                list[index].x = x + comm.x * rounded_division(WIDTH,comm.width);
                list[index].y = y + comm.y * rounded_division(HEIGHT,comm.height);
                index++;
            }
        }

    return nb_points;

}

void save_generation_svg(cell_point* points, int nb_points, FILE* svg){
    for(int i=0; i<nb_points; i++){
        if(points[i].gen == 0)
            fprintf(svg,"<rect width='0' height='1' x='%d' y='%d' fill='black'><animate id='gen%d' attributeName='width' values='1' begin='0s;gen%d.end' dur='%s'/></rect>\n",points[i].x,points[i].y,points[i].gen,ITERATIONS-1,SVG_GEN_DURATION);
        else
            fprintf(svg,"<rect width='0' height='1' x='%d' y='%d' fill='black'><animate id='gen%d' attributeName='width' values='1' begin='gen%d.end' dur='%s'/></rect>\n",points[i].x,points[i].y,points[i].gen,points[i].gen - 1,SVG_GEN_DURATION);
    }
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
void gather_to_one(cellular_grid CG, struct comm_schema comm, int generation, FILE* svg, MPI_Datatype mpi_point){
    // Retrieving points to send
    cell_point ** points = malloc(sizeof(cell_point*));
    int nb_points = generate_svg_points(CG,generation,points,comm);

    // Gathering the number of points each of the processes have
    int *incoming_sizes = malloc(sizeof(int)*comm.size);
    MPI_Gather( &nb_points , 1 , MPI_INT , incoming_sizes , 1 , MPI_INT , comm.master , MPI_COMM_WORLD);

    int max_points = 0;
    int total_points = 0;
    for(int i=0; comm.rank==comm.master && i<comm.size; i++){
        if (max_points<incoming_sizes[i]) max_points = incoming_sizes[i];
        total_points += incoming_sizes[i];
    }

    // Gathering the points
    cell_point * gather_buff;
    int displacements[comm.size]; // Array containing the displacements of each incoming arrays in the buffer from the start of the address.
    if(comm.rank==comm.master){
        gather_buff = malloc(total_points*sizeof(cell_point));
        displacements[0] = 0;
        for (int i=1; i<comm.size; ++i) { 
           displacements[i] = displacements[i-1]+incoming_sizes[i-1]; 
        }
    }
    MPI_Gatherv( *points , nb_points , mpi_point , gather_buff , incoming_sizes , displacements , mpi_point , comm.master , MPI_COMM_WORLD);

    // Saving result in svg
    if(comm.rank==comm.master){
        save_generation_svg(gather_buff,total_points,svg);
        free(gather_buff);
    }
    free(incoming_sizes);
    free(*points);
    free(points);
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
    srand((unsigned) time(&t) + comm.rank);

    for(int i=0; i<local_width*local_height/2; i++){
        set_cell(CG,rand()%local_width,rand()%local_height,1);
    }

    // Creating SVG file
    FILE* fp;
    if(comm.rank == comm.master){
        fp = fopen(SVG_FILE_PATH,"w");
        fprintf(fp,"<svg version='1.1' viewBox='0 0 %d %d' xmlns='http://www.w3.org/2000/svg'>\n",WIDTH,HEIGHT);
        fprintf(fp,"<rect width='%d' height='%d' x='0' y='0' fill='white'/>\n",WIDTH,HEIGHT);
    }

    // Creating MPI structure later used for gathering

    // Displacement of attributes inside the structure
    cell_point dummy_point; // Variable used only to get displacement
    MPI_Aint base_address;
    MPI_Aint displacements[3];
    MPI_Get_address(&dummy_point, &base_address);
    MPI_Get_address(&dummy_point.gen, &displacements[0]);
    MPI_Get_address(&dummy_point.x, &displacements[1]);
    MPI_Get_address(&dummy_point.y, &displacements[2]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);

    // MPI types inside the structure
    MPI_Datatype types[3] = { MPI_INT, MPI_INT, MPI_INT };
    // Length of each block
    int block_lens[3] = { 1, 1, 1 };

    // Finalization of the type
    MPI_Datatype cell_point_type;
    MPI_Type_create_struct(3, block_lens, displacements, types, &cell_point_type);
    MPI_Type_commit(&cell_point_type);
    



    // Main loop
    clock_t start, end;        

    #ifdef V1
    if(comm.rank==comm.master)
        printf("Main loop is starting.\n");
    #endif
    for(int i=0; i<ITERATIONS; i++){
        #ifdef V1
        if(comm.rank==comm.master) start = clock();
        #endif
        // Gather generations points to one process so it can be saved in svg
        gather_to_one(CG,comm,i,fp,cell_point_type); 

        // Next Generation computation
        transmit_walls_V1_0(CG,comm);
        next_generation(CG);

        #ifdef V1
        if(comm.rank==comm.master) {
            end = clock();
            printf("Generation %d done in %lfs.\n",i,(double)(end-start)/CLOCKS_PER_SEC);
        }
        #endif

        // // Time skip
        // usleep(TIME_INTERVAL_U);
    }

    if(comm.rank==comm.master){
        fprintf(fp,"</svg>");
        fclose(fp);
    }

    MPI_Finalize();
    return 0;
}