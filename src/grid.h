#include <stdio.h>

#define bit _Bool

struct _grid{
    bit * value;    // Array containing the bit values of the grid
    uint width;     // Width of the grid
    uint height;    // Height of the grid
    uint size;      // Size of the value array = width*height
};

struct _point{
    uint x;
    uint y;
};


typedef struct _grid* grid;
typedef struct _point* point;

/**
 * @brief Create a grid structure.
 * 
 * @param width Grid's width
 * @param height Grid's height
 * @return grid The created grid
 */
grid create_grid(uint width, uint height);

/**
 * @brief Deletes grid structure.
 * 
 * @param G The grid to delete
 */
void delete_grid(grid G);

/**
 * @brief Gets the value of a bit at a given position of a grid.
 * 
 * @param G The referenced grid
 * @param x Position x of the point to get
 * @param y Position y of the point to get
 * @return int The value read (0 or 1, -1 if invalid position)
 */
int get_bit(grid G, uint x, uint y);

/**
 * @brief Set the value of a bit at a given position of a grid.
 * 
 * @param G The referenced grid
 * @param x Position x of the point to set
 * @param y Position y of the point to set
 * @param new_bit 
 * @return int Status = 1 for no error | -1 invalid position
 */
int set_bit(grid G, uint x, uint y, bit new_bit);

/**
 * @brief Set values of a grid from a set of values.
 * 
 * @param G The referenced grid
 * @param new_values 
 * @return int Status = 1 for no error | -1 invalid position
 */
int set_bits(grid G, bit * new_values);

/**
 * @brief Set values of a grid to 1 from a set of points.
 * 
 * @param G The referenced grid
 * @param p 
 * @param nb_points 
 * @return int Status = 1 for no error | <=0 invalid positions
 */
int set_bits_points(grid G, point * p, uint nb_points);

/**
 * @brief Copy a grid structure.
 * 
 * @param G The referenced grid
 * @return grid The copied grid
 */
grid copy(grid G);