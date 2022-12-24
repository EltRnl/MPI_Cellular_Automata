CC = mpicc
CFLAGS = -g -lm

OBJECTS = grid.o cellular_grid.o automata.o

all: main

main: main.c $(OBJECTS)
	$(CC) $(OBJECTS) main.c -o main $(CFLAGS)

%.o: src/%.c
	$(CC) -c $^ $(CFLAGS)

run : main
	mpirun -np 4 main

clean:
	rm main $(OBJECTS)