CC = gcc
CFLAGS = -O3 -Wall -fPIC -I./include -fopenmp -march=native -ffast-math -ftree-vectorize
LDFLAGS = -shared -fopenmp
LIBNAME = libsuperdiff.a
SONAME = libsuperdiff.so

SRC_DIR = src
OBJ_DIR = build
INC_DIR = include

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

all: $(LIBNAME) $(SONAME)

$(LIBNAME): $(OBJS)
	ar rcs $@ $^

$(SONAME): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(LIBNAME) $(SONAME)

.PHONY: all clean
