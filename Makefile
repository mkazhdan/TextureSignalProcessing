
##export CFLAGS += -DGLM_FORCE_RADIANS=1  # avoid some warnings in include/glm
## For some reason, enabling GLM_FORCE_RADIANS seems to invert the camera

programs = Geodesics LineIntegralConvolution TextureFiltering ReactionDiffusion

# Allow "make -j" to operate in parallel over the programs.
all: $(programs)
$(programs):
	$(MAKE) -C $@

programs_debug = $(foreach n,$(programs),debug_$(n))  # pseudo-dependency to allow "make -j" parallelism
debug: $(programs_debug)
$(programs_debug):
	$(MAKE) -C $(@:debug_%:%) debug

programs_clean = $(foreach n,$(programs),clean_$(n))  # pseudo-dependency to allow "make -j" parallelism
clean: $(programs_clean)
$(programs_clean):
	$(MAKE) -C $(@:clean_%=%) clean

.PHONY: $(programs) $(programs_debug) $(programs_clean)
