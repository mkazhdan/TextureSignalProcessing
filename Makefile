COMPILER ?= gcc
#COMPILER ?= clang
NO_VISUAL ?= false
USE_PARDISO ?= false

ifeq ($(NO_VISUAL),true)
	programs = TextureFiltering TextureStitching TextureMasking TextureDilation GradientDomain.example
else
	programs = Geodesics LineIntegralConvolution TextureFiltering TextureStitching ReactionDiffusion  TextureMasking TextureDilation GradientDomain.example
endif

# Allow "make -j" to operate in parallel over the programs.
all: $(programs)
$(programs):
	$(MAKE) -C $@ COMPILER=$(COMPILER) NO_VISUAL=$(NO_VISUAL) USE_PARDISO=$(USE_PARDISO)

programs_debug = $(foreach n,$(programs),debug_$(n))  # pseudo-dependency to allow "make -j" parallelism
debug: $(programs_debug)
$(programs_debug):
	$(MAKE) -C $(@:debug_%:%) debug COMPILER=$(COMPILER)

programs_clean = $(foreach n,$(programs),clean_$(n))  # pseudo-dependency to allow "make -j" parallelism
clean: $(programs_clean)
$(programs_clean):
	$(MAKE) -C $(@:clean_%=%) clean

.PHONY: $(programs) $(programs_debug) $(programs_clean)
