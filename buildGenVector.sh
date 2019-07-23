#!/bin/bash

CURDIR="$(pwd)"

SOURCE="${CURDIR}/source"
BUILD="${CURDIR}/build/"
INSTALL="${CURDIR}/install"
mkdir -p ${SOURCE} ${BUILD} ${INSTALL}
#BUILDTYPES=(Default Release Debug RelWithDebInfo MinSizeRel)
BUILDTYPES=( Debug Release)

cd ${SOURCE}
### in SOURCE genvector libarry source should be cloned.
git clone https://github.com/AndStorm/genvector.git

    export CC=pgcc
    export CXX=pgc++
    export LD=ld.gold
    export AR=ar
    export CXXFLAGS="-mcmodel=medium -Minfo=all -acc -ta=host"
    CXXFLAGSRELEASE="-O3 -DNDEBUG"
    CXXFLAGSCOMMON="-std=c++17"

  for BUILDTYPE in ${BUILDTYPES[@]}
  do
    PROJECT=genvector
    BUILDDIR="${BUILD}/${PROJECT}/${COMP}/${BUILDTYPE}/"
    INSTALLDIR="${INSTALL}/${PROJECT}/${COMP}/${BUILDTYPE}/"
    mkdir -p ${BUILDDIR} && cd ${BUILDDIR}
    cmake ${SOURCE}/${PROJECT} -DCMAKE_BUILD_TYPE=${BUILDTYPE} -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} -DCMAKE_CXX_FLAGS="${CXXFLAGS}" -DCMAKE_CXX_FLAGS_RELEASE="${CXXFLAGSRELEASE}" -DCMAKE_LINKER=${LD} -DCMAKE_AR="$(which ${AR})"  -DBOOST_ROOT="${BOOSTROOT}" -DBOOSTROOT="${BOOSTROOT}" -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_CXX_FLAGS_DEBUG="-g" -DCMAKE_CXX_STANDARD=17
    cmake --build . -- install
    GENVECTOR_INCLUDEDIR="${INSTALLDIR}/include/"
    GENVECTOR_LIBRARY="${INSTALLDIR}/lib/libgenvector.a"
  done
  


