## This file is part of mcannotate.
##
## mcannotate is an annotation tools for PDB files
## that produce information about nucleotides
## conformation and their spatials interfactions
## Copyright (C) 2008,2009,2010, 2011 Universite de Montreal
##

cmake_minimum_required (VERSION 2.6)

# list the source code
file(GLOB MCANNOTATE_SOURCES_CC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cc)

# add the lib
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable (mcannotate ${MCANNOTATE_SOURCES_CC} )
set_target_properties(mcannotate PROPERTIES VERSION ${MCANNOTATE_VERSION_STRING})
set_target_properties(mcannotate PROPERTIES SKIP_BUILD_RPATH TRUE)
target_link_libraries(mcannotate ${EXT_LIBS})

# add install target
install (TARGETS mcannotate DESTINATION bin)



