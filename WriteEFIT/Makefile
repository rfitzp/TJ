# Makefile for Program WriteEFIT

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

BUILD_DIR = ../bin

F_SOURCE = WriteEFIT.f90
F_OBJECT = $(F_SOURCE:%.f90=%.o)

EXECUTABLE = ../bin/write_efit

$(EXECUTABLE):	 	$(F_OBJECT) $(BUILD_DIR)
			$(F_COMPILER) -o $(EXECUTABLE) $(F_OBJECT) $(F_LIBS)

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

clean: 
			rm -f $(F_OBJECT) $(EXECUTABLE)

