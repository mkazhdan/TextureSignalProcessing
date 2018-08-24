all:
	cd Geodesics               && make
	cd LineIntegralConvolution && make
	cd TextureFiltering        && make
	cd ReactionDiffusion       && make

debug:
	cd Geodesics               && make debug
	cd LineIntegralConvolution && make debug
	cd TextureFiltering        && make debug
	cd ReactionDiffusion       && make debug
	
clean:
	cd Geodesics               && make clean
	cd LineIntegralConvolution && make clean
	cd TextureFiltering        && make clean
	cd ReactionDiffusion       && make clean
