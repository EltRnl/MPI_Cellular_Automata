CC = mpicc
CFLAGS = -g -lm

OBJECTS = grid.o cellular_grid.o automata.o

# Variables
## 0 : No prints
## 1 : Print per generations (time spent)
## 2 : All prints (all processes)
VERBOSE = 1

ifeq ($(VERBOSE),0)
	CFLAGS += -DV0
endif
ifeq ($(VERBOSE),1)
	CFLAGS += -DV1
endif
ifeq ($(VERBOSE),2)
	CFLAGS += -DV2
endif

# Compilation commands

all: main

main: main.c $(OBJECTS)
	$(CC) $(OBJECTS) main.c -o main $(CFLAGS)

%.o: src/%.c
	$(CC) -c $^ $(CFLAGS)

run : main
	mpirun -np 8 --use-hwthread-cpus main

clean:
	rm main $(OBJECTS)