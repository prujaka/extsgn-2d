MOD1=parameters.f90
MOD2=aux.f90
MOD3=model.f90
MOD4=methods.f90
MAIN=main.f90
SRC_DIR=src

FCOMP = gfortran
OBJ = $(SRC_DIR)/$(MOD1) $(SRC_DIR)/$(MOD2) $(SRC_DIR)/$(MOD3) $(SRC_DIR)/$(MOD4) $(SRC_DIR)/$(MAIN)
EXEC = exe

FFLAGS_DEBUG = -g -g0 -fimplicit-none -fcheck=all -fbacktrace
FFLAGS_OPT = -O3

.PHONY: compile
compile:
	@$(FCOMP) $(FFLAGS_OPT) -o $(EXEC) $(OBJ)

.PHONY: debug
debug:
	@$(FCOMP) $(FFLAGS_DEBUG) -o $(EXEC) $(OBJ)

.PHONY: run
run:
	@./$(EXEC)

.PHONY: makerun
makerun: compile run

.PHONY: clean
clean:
	@rm exe
	@rm *mod*
