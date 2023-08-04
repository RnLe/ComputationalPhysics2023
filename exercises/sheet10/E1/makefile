# Makefile for pybind11
ifeq ($(OS),Windows_NT)
	current_os := Windows
else
	current_os := $(shell uname -s)
endif

ifeq ($(current_os),Darwin)
DARWIN_LINKER_FLAGS = -undefined dynamic_lookup
endif

CC = g++
CFLAGS = -std=c++11 -fPIC -Ofast
LINKER_FLAGS = -shared $(DARWIN_LINKER_FLAGS)
INCLUDE = $(shell python3 -m pybind11 --includes)
SUFFIX_PYBIND = $(shell python3-config --extension-suffix)

Solver$(SUFFIX_PYBIND): pybind_wrapper.cpp
	$(CC) $(CFLAGS) $(LINKER_FLAGS) $(INCLUDE) $^ -o Solver$(SUFFIX_PYBIND)

clean:
	rm Solver$(SUFFIX_PYBIND)
