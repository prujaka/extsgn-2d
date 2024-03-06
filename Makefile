FC = gfortran
FCFLAGS = -g -O3
# FLFLAGS = -g -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation \
# -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file \
# -fcheck=all -std=f2008 -pedantic -fbacktrace
EXEC = exe
BLD_DIR = build
SRC_DIR = src

MOD_SRCS = $(wildcard $(SRC_DIR)/m_*.f90)
SRCS = $(wildcard $(SRC_DIR)/*.f90)
OBJS = $(SRCS:$(SRC_DIR)/%.f90=$(BLD_DIR)/%.o)
MODS = $(MOD_SRCS:$(SRC_DIR)/%.f90=$(BLD_DIR)/%.mod)
MOD_OBJS = $(MOD_SRCS:$(SRC_DIR)/%.f90=$(BLD_DIR)/%.o)

# Compile
$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(BLD_DIR)
	$(FC) $(FCFLAGS) -Jbuild -c $< -o $@

# Link
$(BLD_DIR)/$(EXEC): $(OBJS)
	@mkdir -p $(BLD_DIR)
	$(FC) $(FCFLAGS) -o $@ $(OBJS)

# Dependencies
# $(BLD_DIR)/main.o : $(MOD_OBJS)
$(BLD_DIR)/main.o: $(BLD_DIR)/aux.o $(BLD_DIR)/methods.o $(BLD_DIR)/model.o $(BLD_DIR)/parameters.o
$(BLD_DIR)/aux.o $(BLD_DIR)/methods.o $(BLD_DIR)/model.o: $(BLD_DIR)/parameters.o
$(BLD_DIR)/methods.o: $(BLD_DIR)/model.o $(BLD_DIR)/aux.o


run: $(BLD_DIR)/$(EXEC)
	$(BLD_DIR)/$(EXEC)

.PHONY: clean debug plot
clean:
	$(RM) -r $(BLD_DIR)
	$(RM) $(SRC_DIR)/*.mod

debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "EXEC = $(EXEC)"

plot:
	@python postprocessing/plot.py
