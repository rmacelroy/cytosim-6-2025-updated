# Linux installation

On Linux, you need to install the GNU compiler collection, BLAS and LApack. This is sufficient to compile sim. Recent Linux distributions provide precompiled BLAS/LAPACK as optional installation. Please, check your distribution.

# Prerequisite

Update your linux and packages:

	sudo apt update

To install the required tools and libraries with `Ubuntu`:

	sudo apt-get install make g++
	sudo apt-get install libncurses-dev
	sudo apt-get install libblas.dev liblapack.dev
	sudo apt-get install freeglut3-dev libxi-dev libxmu-dev libglew-dev

# Docker

FROM: ubuntu:18.04
MAINTAINER: SergeDmi www.biophysics.fr

	apt-get update
	apt-get install -y libncurses-dev
	apt-get install -y libblas.dev liblapack.dev 
	apt-get install -y freeglut3-dev libxi-dev libxmu-dev libglew-dev
	apt-get install -y git make
	apt-get install -y g++

	git clone https://gitlab.com/f-nedelec/cytosim.git cytosim
	
 RUN in WORKDIR cytosim
	
    make -j4


# Contact

Maintained by Serge Dmitrieff

