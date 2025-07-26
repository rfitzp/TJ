# Makefile for Class Utility

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

all:
	$(C_COMPILER) $(C_FLAGS) -o Utility.o -c $(C_INCLUDE) Utility.cpp
	$(C_COMPILER) $(C_FLAGS) -o Faddeeva.o -c $(C_INCLUDE) Faddeeva.cpp

clean: 
	rm -f Utility.o Faddeeva.o

