CC = mpicc
CFLAGS = -g -Wall
LDFLAGS = -lm
VARFLAGS = 

OBJECTS = grid.o cellular_grid.o rendering.o automata.o

# Variables
## How verbose the application is :
## 		0 : No prints
## 		1 : Print per generations (time spent)
## 		2 : All prints (all processes)
VERBOSE = 0

ifeq ($(VERBOSE),0)
	VARFLAGS += -DV0
endif
ifeq ($(VERBOSE),1)
	VARFLAGS += -DV1 -DV0
endif
ifeq ($(VERBOSE),2)
	VARFLAGS += -DV2 -DV1 -DV0
endif

## Choice of display :
## 		x11 : X11 display on screen
## 		svg : Saving as SVG file
## 		default : No display or file save (used for performance mesurement)
DISPLAY_MODE = x11

ifeq ($(DISPLAY_MODE), x11)
	VARFLAGS += -DX11
	LDFLAGS += -L/usr/X11/lib -lX11
else ifeq ($(DISPLAY_MODE), svg)
	VARFLAGS += -DSVG
else
	VARFLAGS += -DNORENDER
endif

# Compilation commands

all: main

main: main.c $(OBJECTS)
	$(CC) $(VARFLAGS) $(OBJECTS) main.c -o main $(CFLAGS) $(LDFLAGS)

%.o: src/%.c
	$(CC) $(VARFLAGS) -c $^ $(CFLAGS) 

run : main
	mpirun -np 8 --use-hwthread-cpus -x DISPLAY=:0 main

clean:
	rm main $(OBJECTS)