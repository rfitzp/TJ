# Makefile for Program TJ

all:
	(cd TJ;        make)
	(cd WriteEFIT; make)
	(cd Flux;      make)
	(cd Tear;      make)

clean:

	(cd TJ;        make clean)
	(cd WriteEFIT; make clean)
	(cd Flux;      make clean)
	(cd Tear;      make clean)
