#!/bin/bash

export SITE=`pwd`

 ../src/pylith/pylith_installer/configure \
     --enable-python=yes \
     --enable-cython=yes \
     --enable-six=yes \
     --enable-numpy=yes \
     --enable-cftime=yes \
     --with-petsc-options='--with-debugging=0 --download-f2cblaslapack=1 --with-COPTFLAGS="-g -march=native"' \
     --enable-swig=yes \
     --with-debugging=yes \
     --enable-mpi=mpich \
     --with-pylith-repo=https://github.com/rwalkerlewis/pylith.git \
     --with-pylith-git=multi-fault \
     --with-make-threads=2 \
    --prefix=$SITE/pylith \
    CFLAGS="-g -march=native" CXXFLAGS="-g -march=native -std=c++11" \
    --enable-force-install


##     --with-pylith-git=baagaard/fix-poroelasticity-fullscale-tests \    
##     --with-pylith-git=knepley/feature-petsc-fe \
## CPPFLAGS="-std=c++11"
# Pylith
#     --with-pylith-repo=https://github.com/rwalkerlewis/pylith.git \
#     --with-pylith-git=rlwalker/cryer-test \

#     --with-pylith-repo=https://github.com/rwalkerlewis/pylith.git \
#     --with-pylith-git=rlwalker/poroelastic-mms \
#     --with-pylith-repo=https://github.com/rwalkerlewis/pylith.git \
#     --with-pylith-git=rlwalker/mms \
