# Compiler and flags
CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -I/usr/local/include
LDFLAGS = -L/usr/local/lib
LDLIBS  = -ltrexio

# Output binary and source file
BIN = mp2_energy
SRC = src/main.c

# Default target
all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) $(LDLIBS) -o $@

# Remove compiled files
clean:
	rm -f $(BIN)

.PHONY: all clean

# Simple test target
test: mp2_energy
	./tests/run_tests.sh

