# Makefile for Program TJ

all:
	(cd Utility;     make)
	(cd Equilibrium; make)
	(cd Layer;	 make)
	(cd WriteEFIT;   make)
	(cd Flux;        make)
	(cd Tear;        make)
	(cd Island;      make)
	(cd ECE;      	 make)
	(cd StartUp;     make)
	(cd Vertical;    make)
	(cd TJ;        	 make)
	(cd TearX;       make)

clean:
	(cd Utility;     make clean)
	(cd Equilibrium; make clean)
	(cd Layer;	 make clean)
	(cd WriteEFIT;   make clean)
	(cd Flux;        make clean)
	(cd Tear;        make clean)
	(cd Island;      make clean)
	(cd ECE;      	 make clean)
	(cd StartUp;     make clean)
	(cd Vertical;    make clean)
	(cd TJ;          make clean)
	(cd TearX;       make clean)
