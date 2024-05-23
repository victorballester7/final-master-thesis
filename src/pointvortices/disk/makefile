# Description: Makefile for the project

# Compiler and flags
CXX := g++
CC := gcc
CXXFLAGS := -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -g3 -Ofast -std=c++11
CFLAGS := -Wall -Wextra -pedantic -Ofast -std=c11

# Libraries and paths
LIBS := -lm
INCLUDES := 

# Folders
SRC := .
SRC2 := ..
INCLUDE := ../../../include
BIN := ../../../bin/pointvortices/disk

# Executable name
TARGET := main

# Sources and objects
SOURCES := $(wildcard $(SRC)/*.c) $(wildcard $(SRC)/*.cpp) $(wildcard $(SRC2)/*.c) $(wildcard $(SRC2)/*.cpp)
# OBJECTS := $(patsubst $(SRC)/%.cpp, $(BIN)/%.o, $(filter %.cpp, $(SOURCES))) $(patsubst $(SRC)/%.c, $(BIN)/%.o, $(filter %.c, $(SOURCES))) $(patsubst $(SRC2)/%.cpp, $(BIN)/%.o, $(filter %.cpp, $(SOURCES))) $(patsubst $(SRC2)/%.c, $(BIN)/%.o, $(filter %.c, $(SOURCES)))
OBJECTS := $(patsubst $(SRC)/%.cpp, $(BIN)/%.o, $(filter $(SRC)/%.cpp, $(SOURCES))) $(patsubst $(SRC)/%.c, $(BIN)/%.o, $(filter $(SRC)/%.c, $(SOURCES))) $(patsubst $(SRC2)/%.cpp, $(BIN)/%.o, $(filter $(SRC2)/%.cpp, $(SOURCES))) $(patsubst $(SRC2)/%.c, $(BIN)/%.o, $(filter $(SRC2)/%.c, $(SOURCES)))

# Headers
HEADERS := $(wildcard $(INCLUDE)/*.h) $(wildcard $(INCLUDE)/*.hpp)

.PHONY: all clean

all: $(BIN)/$(TARGET)

$(BIN)/$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(BIN)/%.o: $(SRC)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN)/%.o: $(SRC2)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BIN)/%.o: $(SRC)/%.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(BIN)/%.o: $(SRC2)/%.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -f $(BIN)/*.o $(BIN)/$(TARGET)
