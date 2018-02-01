DATE=`date +%Y%m%d-%H%M`
MYDIR=$(shell pwd)
MYCASE=$(notdir $(CURDIR))

EXEDIR=${MYDIR}/bin
EX=aikef_mpi
EXE=bin/${EX}

NODES=4
PPN=24
RUNTIME=12:00:00
MAIL=j.riousset@tu-bs.de

PROCS= $(shell echo ${NODES}*${PPN} | bc)

SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:src/%.cpp=build/%.o)

#Numbers of cores to use for compilation
MAKEFLAGS+="-j 4"

#Number of cores to use / in case of cluster this is the number of cores per node in case of home this ist the number of cores/cpus to use
CORES=4


INCLUDE  := \
	-Iinclude/ \
	-Iinclude/vectorclass \
	-I/usr/local/include \
	-I/usr/local/Cellar/cspice/66/include \
	-I/usr/local/Cellar/gsl/2.4/include

LIB :=  \
	-L/usr/local/lib \
	-L/usr/local/Cellar/gsl/2.4/lib/ \
	-L/usr/local/Cellar/cspice/66/lib/ \
	-lsiloh5 -lhdf5 -lhdf5_hl -lgsl -lgslcblas -lz -lm -ldl -lc #-static-libgcc

OMPI_CC := clang
OMPI_CXX := clang++
ARCH := -march=native -mtune=native -m64 -msse4 -msse4.1 -msse4.2
OPTIMIZE := $(ARCH) -O3 -fno-builtin -fno-rtti -ffunction-sections -fno-signed-zeros
STANDARD := -Wall -Wextra -pedantic -std=c++11
COMPILER := /usr/local/Cellar/open-mpi/3.0.0_2/bin/mpicxx

export OMPI_CC
export OMPI_CXX

all: $(OBJS) $(EXE)

$(EXE): $(OBJS)
	$(COMPILER) $(INCLUDE) $(OBJS) -o $(EXE) $(OPTIMIZE) $(STANDARD) $(LIB)

build/%.o: src/%.cpp
	$(COMPILER) $(INCLUDE) -c $< -o $@ $(OPTIMIZE) $(STANDARD)


clean:
	find . \( -name "*.o" -o -name "${EX}" -o -name "process_*.log" \) -exec rm -rfv {} \;
	rm -rfv bin/*

run:
		#@echo Deleting old files...
		#rm -f bin/State/*
		#rm -f bin/Last_State/*
		#rm -f $(wildcard bin/process_*.log)
		@echo Starting AIKEF...
		#mpirun -np $(CORES) -wdir ${EXEDIR} nice -19 ./${EX} &
		/usr/local/Cellar/open-mpi/3.0.0_2/bin/mpirun -np 2 bin/./aikef_mpi > log/aikef.log &
		#msub -l walltime=${RUNTIME},nodes=${NODES}:ppn=${PPN} -m abe -M ${MAIL} -j oe -o ${EXEDIR}/log_${MYCASE}.out -v EXEDIR=${EXEDIR},EX=${EX},MYCASE=${MYCASE} start.job > ${EXEDIR}/JOB.ID

stop:
	ifeq ($(TYPE),$(filter $(TYPE),home_gcc home_icc home_gcc_O3 home_icc_O3 home_clang home_clang_O3))
		skill ${EX}
	else ifeq ($(TYPE),$(filter $(TYPE),hlrn_gcc hlrn_icc juropa_gcc juropa_icc))
		@echo Aborting AIKEF...
		canceljob `cat ${EXEDIR}/JOB.ID`
	endif

rmstate:
	rm -r bin/State
	rm -r bin/Last_State
	rm $(wildcard bin/process_*.log)

rmdata:
	find .    -name ".*.sw*" -exec rm -rfv {} \;
	find output -name "*.*"    -exec rm -rfv {} \;
	find bin  \( -name "*log*" -o -name "JOB.ID" -o -name "State*" -o -name "Last_State" \) -exec rm -rfv {} \;

clear: clean rmdata

rmdir:
	rm -r bin output build

make_doc:
	@echo Doxygen is needed for automatic documentation generation
	doxygen Doxyfile

info:
	@echo
	@echo COMPILER = $(COMPILER)
	@echo OMPI_CC  = $(OMPI_CC)
	@echo OMPI_CXX = $(OMPI_CXX)
	@echo
	@echo INCLUDE  = $(INCLUDE)
	@echo
	@echo OPTIMIZE = $(OPTIMIZE)
	@echo
	@echo STANDARD = $(STANDARD)
	@echo
	@echo LIB  = $(LIB)
	@echo
	@echo SRCS = $(SRCS)
	@echo
	@echo OBJS = $(OBJS)
	@echo
	@echo

mkdir:
	-mkdir bin
	-mkdir output
	-mkdir build
	-mkdir output/lineout
	-mkdir output/particle_detector
	-mkdir output/particle_tracks
	-mkdir output/uniform_output
	-mkdir output/silo
	-mkdir output/silo_3D
	-mkdir output/trajectories

