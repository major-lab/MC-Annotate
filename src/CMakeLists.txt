## This file is part of mcannotate.
##
## mcannotate is an annotation tools for PDB files that produce information about nucleotides 
## conformation and their spatials interfactions
## Copyright (C) 2008,2009,2010, 2011 Institut de recherche en immunologie et en cancérologie
##
## mcannotate is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## mcannotate is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with mcannotate  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.6)

# liste de tous les fichiers source
file(GLOB MCANNOTATE_SOURCES_CC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cc)

# ajoute la librarie
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable (mcannotate ${MCANNOTATE_SOURCES_CC} )
set_target_properties(mcannotate PROPERTIES VERSION ${MCANNOTATE_VERSION_STRING})
set_target_properties(mcannotate PROPERTIES SKIP_BUILD_RPATH TRUE)
target_link_libraries(mcannotate ${EXT_LIBS})

# ajoute le target d'installation
install (TARGETS mcannotate DESTINATION bin)



