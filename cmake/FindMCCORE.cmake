## This file is part of mcsym.
##
## mcsym is a computer program that predicts all-atoms RNA tertiary structure from 
## sequence and secondary structure
## Copyright (C) 2008,2009,2010,2011 Université de Montréal
##


##
## This module finder will define the following variables
##
## MCCORE_FOUND 
## MCCORE_INCLUDE_DIRS
## MCCORE_LIBRARIES
##

include(LibFindMacros)

if (HANDLE_GCC_VAR)
  find_path(MCCORE_INCLUDE_DIR  mccore/Version.h
    PATHS
    /usr/include
    /usr/local/include
    ENV CPLUS_INCLUDE_PATH)

  find_library(MCCORE_LIBRARY mccore
    PATHS
    /usr/lib
    /usr/local/lib
    ENV LIBRARY_PATH )
else()
  find_path(MCCORE_INCLUDE_DIR  mccore/Version.h
    PATHS
    /usr/include
    /usr/local/include)

  find_library(MCCORE_LIBRARY mccore
    PATHS
    /usr/lib
    /usr/local/lib)
endif()

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(MCCORE_PROCESS_INCLUDES MCCORE_INCLUDE_DIR)
set(MCCORE_PROCESS_LIBS MCCORE_LIBRARY)
libfind_process(MCCORE)

# advanced view only
mark_as_advanced(MCCORE_INCLUDE_DIR MCCORE_LIBRARY)



