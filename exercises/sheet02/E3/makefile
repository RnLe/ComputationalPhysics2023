CXX = gcc
CXXFLAGS = -Wall -O0 -std=c++20

SRC_DIR = .
OBJ_DIR = ./build/Debug
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(SRC_FILES:%.cpp=$(OBJ_DIR)/%.o)

EXECUTABLE = outDebug

all: $(EXECUTABLE)

$(OBJ_FILES): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(EXECUTABLE): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lstdc++ -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: clean run

clean:
	rm -f $(OBJ_DIR)/*.o $(EXECUTABLE)

args = 
run: $(EXECUTABLE)
	./$(EXECUTABLE) $(args)