/* This file is part of mcannotate.
 *
 * mcsym is a computer program that predicts all-atoms RNA tertiary structure from 
 * sequence and secondary structure
 * Copyright (C) 2008,2009,2010,2011 Institut de recherche en immunologie et en cancérologie
 *
 * mcannotate is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * mcannotate is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with mcannotate.  If not, see <http://www.gnu.org/licenses/>.
 */  

#ifndef CONFIG_H
#define CONFIG_H

// package information
#define PACKAGE_NAME "${PACKAGE_NAME}"
#define PACKAGE_VERSION_STRING "${MCANNOTATE_VERSION_STRING}"

// checks for functions

// needed for actual version handling of Version.cc
#define VERSION_CPU "${CMAKE_SYSTEM_PROCESSOR}"
#define VERSION_OS "${CMAKE_SYSTEM_NAME}"

#endif /*CONFIG_H*/