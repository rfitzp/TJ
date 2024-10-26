# Makefile for Program TJ

all:
	(cd Equilibrium; make)
	(cd Layer;	 make)
	(cd TJ;        	 make)
	(cd WriteEFIT;   make)
	(cd Flux;        make)
	(cd Tear;        make)

clean:

	(cd Equilibrium; make clean)
	(cd Layer;	 make clean)
	(cd TJ;          make clean)
	(cd WriteEFIT;   make clean)
	(cd Flux;        make clean)
	(cd Tear;        make clean)
