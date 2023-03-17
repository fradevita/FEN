#!/bin/bash

# Define environmental variables
# (NOTE: considering to add this line to your bashrc to avoid exporting every time)
export FEN_DIR=??
export _2DECOMP_DIR=??
export FFTW3_DIR=??

# Install FFTW3 library
if [ ! -d "$FFTW3_DIR" ];
then
	echo "Please provide a correct path to the variable FFTW3_DIR"
	exit 1
fi
cd $FFTW3_DIR
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar -zxvf fftw-3.3.10.tar.gz
cd fftw-3.3.10
./configure --prefix=$FFTW3_DIR
make
make install

# Install 2decomp&FFT library
if [ ! -d "${_2DECOMP_DIR}" ];
then
        echo "Please provide a correct path to the variable _2DECOMP_DIR"
        exit 1
fi
git clone https://github.com/numericalalgorithmsgroup/2decomp_fft.git ${_2DECOMP_DIR}
cd ${_2DECOMP_DIR}
sed -i '13s/OPTIONS=/OPTIONS= -DDOUBLE_PREC/' src/Makefile.inc
sed -i '62s/OPTIM=-O3/OPTIM=-O2/' src/Makefile.inc
make

echo 'Installation completed.'
