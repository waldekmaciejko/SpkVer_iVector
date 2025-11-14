cmake_minimum_required(VERSION 3.19)
project(SpkVer_iVector LANGUAGES CXX)

find_package(Qt6 6.5 REQUIRED COMPONENTS Core)

qt_standard_project_setup()

qt_add_executable(SpkVer_iVector
    #__main.cpp
    doc/doxyconfig.cfg
    src/enrollFRRFAR.cpp
    src/filelists.cpp
    src/funhelpers.cpp
    src/gpldamodel.cpp
    src/Log4AVA.cpp
    src/mfccefeatures.cpp
    src/mfccextractor.cpp
    src/trainTVPMiVector.cpp
    src/unittests.cpp
    src/config-fileparser.cpp
    WORKSPACE/calc/placeholder.stayempty
    WORKSPACE/cfg/placeholder.stayempty
    WORKSPACE/feat/enroll/placeholder.stayempty
    WORKSPACE/feat/train/placeholder.stayempty
    WORKSPACE/feat/placeholder.stayempty
    WORKSPACE/list/placeholder.stayempty
    WORKSPACE/logs/placeholder.stayempty
    README.md
    include/enrollFRRFAR.h
    include/filelists.h
    include/funhelpers.h
    include/gpldamodel.h
    include/Log4AVA.h
    include/mfccefeatures.h
    include/mfccextractor.h
    include/trainTVPMiVector.h
    include/unittests.h
    include/config-fileparser.h
    setup
    include/progressbar.h
    #_main.cpp
    examples/config.cfg
    src/config-fileparser.cpp
    main.cpp
    examples/workspace.cfg
    examples/data.cfg
    include/timeheleprs.h
    src/timehelpers.cpp

)

## ===========
#
# MANUAL SEARCH lapack and blas ON UBUNTU SYSTEM
# FIRST INSTALL COMMAND IN SETUP FILE
#
# dpkg --get-selections | grep -v deinstall | grep blas
# dpkg -L libblas-dev
#
# dpkg --get-selections | grep -v deinstall | grep lapack
# dpkg -L liblapack-dev
#
## ===========


target_link_libraries(SpkVer_iVector
    PRIVATE
        Qt::Core
        /usr/lib/x86_64-linux-gnu/lapack/liblapack.so
        /usr/lib/x86_64-linux-gnu/blas/libblas.so
)

target_include_directories(SpkVer_iVector PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(SpkVer_iVector PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

include(GNUInstallDirs)

# INSTALATION EIGEN

#find_package (Eigen3 3.4.0 REQUIRED NO_MODULE)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/libs/eigen-3.4.0 )

# INSTALATION ARMADILLO
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/libs/armadillo-15.2.1/include )

# INSTALATION QCUSTOMPLOT

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/libs/qcustomplot )

install(TARGETS SpkVer_iVector
    BUNDLE  DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# ldd <binary> to see dependencies
# before use linuxdeployqt-continuous-x86_64 add to system variable path to Qt:
# export PATH=/home/ai/Qt/6.6.3/gcc_64/bin/:$PATH



