CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -fPIC
BUILD_DIR := build
SRC := cpp/polynomial.cpp
HDR := cpp/polynomial.h
PYBIND_SRC := cpp/bindings.cpp

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
  LIB_EXT := so
  EXE_EXT :=
endif

ifeq ($(OS),Windows_NT)
  LIB_EXT := dll
  EXE_EXT := .exe
endif

OBJ_DIR := $(BUILD_DIR)/obj
LIB := $(BUILD_DIR)/libpoly.$(LIB_EXT)
TEST_EXE := $(BUILD_DIR)/test_polynomial$(EXE_EXT)
PYEXT := $(BUILD_DIR)/polylib_ext$(EXE_EXT)

.PHONY: all lib test test-cpp test-py clean

all: lib

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(OBJ_DIR): $(BUILD_DIR)
	mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/polynomial.o: $(SRC) $(HDR) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -DBUILDING_POLY -c $(SRC) -o $@

lib: $(LIB)

$(LIB): $(OBJ_DIR)/polynomial.o | $(BUILD_DIR)
	$(CXX) -shared $(OBJ_DIR)/polynomial.o -o $(LIB)
	mkdir -p src/polylib/lib
	cp $(LIB) src/polylib/lib/

pyext: $(OBJ_DIR)/polynomial.o $(OBJ_DIR)/bindings.o | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -shared $(OBJ_DIR)/polynomial.o $(OBJ_DIR)/bindings.o -o $(BUILD_DIR)/polylib_ext$(if $(filter $(LIB_EXT),dylib),.dylib,.so) `python3 -m pybind11 --includes`

test: test-cpp test-py

test-cpp: $(OBJ_DIR)/polynomial.o | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) cpp/tests/test_polynomial.cpp $(OBJ_DIR)/polynomial.o -o $(TEST_EXE)
	$(TEST_EXE)

test-py: lib
	python -m unittest discover -s tests -p "test*.py" -v

clean:
	rm -rf $(BUILD_DIR)
	rm -rf __pycache__ src/__pycache__

$(OBJ_DIR)/bindings.o: $(PYBIND_SRC) $(HDR) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $(PYBIND_SRC) -o $@ `python3 -m pybind11 --includes`
