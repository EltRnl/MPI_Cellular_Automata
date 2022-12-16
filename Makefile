CC = gcc


all: main

main: main.c grid.o cellular_grid.o
	$(CC) grid.o cellular_grid.o main.c -o main -g

%.o: src/%.c
	$(CC) -c $^ -g

run : $(TARGET)
	./$(TARGET)

clean:
	rm main *.o