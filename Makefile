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
$(BLD_DIR)/main.o : $(MOD_OBJS)
$(BLD_DIR)/m_aux.o $(BLD_DIR)/m_methods.o $(BLD_DIR)/m_model.o: $(BLD_DIR)/m_parameters.o
$(BLD_DIR)/m_methods.o: $(BLD_DIR)/m_model.o $(BLD_DIR)/m_aux.o


run: $(BLD_DIR)/$(EXEC)
	$(BLD_DIR)/$(EXEC)

.PHONY: clean debug plot generate_images rm_seq rm_img_seq generate_video \
	generate_matrix
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
	@mkdir -p img
	@python plot/xsgnplot.py

generate_images:
	@python plot/plot-sequence.py

generate_video:
	@mkdir -p vid
	@python plot/make-video.py

generate_matrix:
	@python plot/pic-data.py

rm_seq:
	@rm -f out/*_t=*.dat

rm_img_seq:
	@rm -f img/*_t=*.png

clean_sequences:
	@make rm_seq
	@make rm_img_seq

complete_video_circuit_artsy:
	@make clean_sequences
	@make generate_matrix
	@make run
	@make generate_images
	@make generate_video
