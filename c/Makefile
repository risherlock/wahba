
SRC = main.c esoq2.c

CC = gcc
CFLAGS = -Wall -lm
OBJ = $(SRC:%.c=build/%.o)
TARGET = build/a.out

all: $(TARGET)

$(TARGET): $(OBJ) | build
	$(CC) $(OBJ) $(CFLAGS) -o $(TARGET)

build/%.o: %.c | build
	$(CC) $(CFLAGS) -c $< -o $@

build:
	mkdir -p build

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf build

.PHONY: all run clean
