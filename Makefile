TARGET_NAME = extsgn
BUILD_DIR   = build
SRC_DIR     = src

MKDIR_P ?= mkdir -p

SRCS := $(wildcard $(SRC_DIR)/*.f90)
OBJS := $(SRCS:$(SRC_DIR)/%.f90=$(BUILD_DIR)/%.o)

FC      = gfortran
FCFLAGS = -g -O2 -J $(BUILD_DIR)
LDFLAGS =

$(BUILD_DIR)/$(TARGET_NAME): $(OBJS)
	$(FC) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 $(BUILD_DIR)/%.mod
	@$(MKDIR_P) $(dir $@)
	$(FC) $(FCFLAGS) -c $< -o $@

$(BUILD_DIR)/%.mod: $(SRC_DIR)/%.f90
	@$(MKDIR_P) $(dir $@)
	$(FC) $(FCFLAGS) -c $< -o $(BUILD_DIR)/$*.o

$(BUILD_DIR)/aux.mod: $(BUILD_DIR)/parameters.mod
$(BUILD_DIR)/main.mod: $(BUILD_DIR)/aux.mod $(BUILD_DIR)/parameters.mod \
	$(BUILD_DIR)/methods.mod
$(BUILD_DIR)/model.mod: $(BUILD_DIR)/parameters.mod
$(BUILD_DIR)/methods.mod: $(BUILD_DIR)/parameters.mod $(BUILD_DIR)/model.mod \
	$(BUILD_DIR)/aux.mod

$(BUILD_DIR)/aux.o: $(BUILD_DIR)/parameters.mod
$(BUILD_DIR)/main.o: $(BUILD_DIR)/aux.mod $(BUILD_DIR)/parameters.mod \
	$(BUILD_DIR)/methods.mod
$(BUILD_DIR)/model.o: $(BUILD_DIR)/parameters.mod
$(BUILD_DIR)/methods.o: $(BUILD_DIR)/parameters.mod $(BUILD_DIR)/model.mod \
	$(BUILD_DIR)/aux.mod

.PHONY: all
all: $(BUILD_DIR)/$(TARGET_NAME)

.PHONY: clean
clean:
	$(RM) -r $(BUILD_DIR)
	$(RM) $(SRC_DIR)/*.mod

.PHONY: run
run: $(BUILD_DIR)/$(TARGET_NAME)
	cd $(BUILD_DIR); ./$(TARGET_NAME)

.PHONY: debug
debug:
	@echo SRCS = $(SRCS)
	@echo OBJS = $(OBJS)

.PHONY: contour_plot
contour_plot:
	@python postprocessing/contourplot.py

.PHONY: plot
plot:
	@python postprocessing/plot.py
