#!/bin/bash

cd /opt/pylith
mkdir src
mkdir -p ${TOPBUILD_DIR} && pushd ${TOPBUILD_DIR} && mkdir pythia spatialdata pylith && popd
mkdir -p ${INSTALL_DIR}
cd /opt/pylith/src
git clone --recursive https://github.com/geodynamics/pythia.git
git clone --recursive https://github.com/geodynamics/spatialdata.git
git clone --recursive --branch multi-fault https://github.com/rwalkerlewis/pylith.git
git clone --branch knepley/pylith-dev https://gitlab.com/petsc/petsc.git
# PyLith repository (repeat for other repositories you forked)
cd /opt/pylith/src/pylith
git remote add upstream https://github.com/geodynamics/pylith.git
# Set URLs for submodules in `.git/config` to geodynamics repository (PyLith repository).
cd /opt/pylith/src/pylith
git config submodule.m4.url https://github.com/geodynamics/autoconf_cig.git
git config submodule.templates/friction/m4.url https://github.com/geodynamics/autoconf_cig.git
git config submodule.templates/materials/m4.url https://github.com/geodynamics/autoconf_cig.git

# Update submodules
git submodule update

# Pythia
cd ${TOPBUILD_DIR}/pythia
pushd ${TOPSRC_DIR}/pythia && autoreconf -if && popd
${TOPSRC_DIR}/pythia/configure --prefix=${INSTALL_DIR} --enable-testing \
    CC=mpicc CXX=mpicxx CFLAGS="-g -Wall" CXXFLAGS="-std=c++11 -g -Wall"
make install
make check

# Spatialdata
cd ${TOPBUILD_DIR}/spatialdata
pushd ${TOPSRC_DIR}/spatialdata && autoreconf -if && popd
${TOPSRC_DIR}/spatialdata/configure --prefix=${INSTALL_DIR} \
    --enable-swig --enable-testing --enable-test-coverage \
    --with-python-coverage=python3-coverage \
	CPPFLAGS="-I${DEPS_DIR}/include -I${INSTALL_DIR}/include" \
	LDFLAGS="-L${DEPS_DIR}/lib -L${INSTALL_DIR}/lib --coverage" \
	CXX=mpicxx CXXFLAGS="-std=c++11 -g -Wall --coverage"
make install -j$(nproc)
make check -j$(nproc)

# PETSc
cd ${TOPSRC_DIR}/petsc
python3 ./configure --with-c2html=0 --with-lgrind=0 --with-fc=0 \
    --with-x=0 --with-clanguage=C --with-mpicompilers=1 \
    --with-shared-libraries=1 --with-64-bit-points=1 --with-large-file-io=1 \
    --with-hdf5=1 --download-chaco=1 --download-ml=1 \
    --download-f2cblaslapack=1 --with-debugging=1 CFLAGS="-g -O -Wall" \
    CPPFLAGS="-I${HDF5_INCDIR} -I${DEPS_DIR}/include -I${INSTALL_DIR}/include" \
    LDFLAGS="-L${HDF5_LIBDIR} -L${DEPS_DIR}/lib -L${INSTALL_DIR}/lib"
make 
make check

# Pylith
cd ${TOPBUILD_DIR}/pylith
pushd ${TOPSRC_DIR}/pylith && autoreconf -if && popd
${TOPSRC_DIR}/pylith/configure --prefix=${INSTALL_DIR} \
    --enable-cubit --enable-hdf5 --enable-swig --enable-testing \
    --enable-test-coverage --with-python-coverage=python3-coverage \
    CPPFLAGS="-I${HDF5_INCDIR} -I${DEPS_DIR}/include -I${INSTALL_DIR}/include" \
    LDFLAGS="-L${HDF5_LIBDIR} -L${DEPS_DIR}/lib -L${INSTALL_DIR}/lib --coverage" \
    CC=mpicc CFLAGS="-g -Wall" CXX=mpicxx CXXFLAGS="-std=c++11 -g -Wall --coverage"
make install -j$(nproc)
make check -j$(nproc)

