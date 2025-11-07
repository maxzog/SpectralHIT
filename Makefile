# Compiler and flags
FC = gfortran
FFLAGS = -O3 -march=native -J./obj -I./obj -I$(FFTW_DIR)/include -L$(FFTW_DIR)/lib -lfftw3 -lgomp -lm -Wall -ffast-math -fopenmp

# Directories
SRC_DIR = ./src
CASE_DIR = ./cases
OBJ_DIR = ./obj

# Source files in src
SRC_FILES = $(wildcard $(SRC_DIR)/*.f90)

# Object files for src
OBJ_FILES = $(patsubst $(SRC_DIR)/%.f90,$(OBJ_DIR)/%.o,$(SRC_FILES))

# Case files
CASE_FILES = q1 q2 q3 q4 q5 q6

# Default target: show help message
all:
	@echo "Specify one of the following targets: $(CASE_FILES)"

# Pattern rule for compiling each case file
$(CASE_DIR)/%: $(OBJ_DIR)/%.o $(OBJ_FILES)
	$(FC) $(FFLAGS) -o $@ $^

# Compile object file for each case file in ./cases
$(OBJ_DIR)/%.o: $(CASE_DIR)/%.f90 $(OBJ_FILES)
	$(FC) $(FFLAGS) -c $< -o $@

# Compile source files in src
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Targets for each case file
q1: $(CASE_DIR)/q1
q2: $(CASE_DIR)/q2
q3: $(CASE_DIR)/q3
q4: $(CASE_DIR)/q4
q5: $(CASE_DIR)/q5
q6: $(CASE_DIR)/q6

# Clean the build
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*.mod $(CASE_DIR)/q1 $(CASE_DIR)/q2 $(CASE_DIR)/q3 $(CASE_DIR)/q4 $(CASE_DIR)/q5 $(CASE_DIR)/q6

# Phony targets
.PHONY: all clean q1 q2 q3 q4 q5 q6
