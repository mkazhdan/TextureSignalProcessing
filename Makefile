COMPILER ?= gcc
#COMPILER ?= clang
NO_VISUAL ?= false

ifeq ($(NO_VISUAL),true)
	programs = TextureFiltering TextureStitching TextureMasking TextureDilation
else
	programs = Geodesics LineIntegralConvolution TextureFiltering TextureStitching ReactionDiffusion  TextureMasking TextureDilation
endif

# Allow "make -j" to operate in parallel over the programs.
all: $(programs)
$(programs):
	$(MAKE) -C $@ COMPILER=$(COMPILER) NO_VISUAL=$(NO_VISUAL)

programs_debug = $(foreach n,$(programs),debug_$(n))  # pseudo-dependency to allow "make -j" parallelism
debug: $(programs_debug)
$(programs_debug):
	$(MAKE) -C $(@:debug_%:%) debug COMPILER=$(COMPILER)

programs_clean = $(foreach n,$(programs),clean_$(n))  # pseudo-dependency to allow "make -j" parallelism
clean: $(programs_clean)
$(programs_clean):
	$(MAKE) -C $(@:clean_%=%) clean

.PHONY: $(programs) $(programs_debug) $(programs_clean)
