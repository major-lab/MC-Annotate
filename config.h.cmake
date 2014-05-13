// This file is part of mcannotate.

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