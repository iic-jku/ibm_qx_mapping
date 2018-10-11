# Try to find the MPFR C++ interface
# See http://www.holoborodko.com/pavel/mpfr/
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPFR++ 2.3.0)
# to require version 2.3.0 to newer of MPFR C++.
#
# Once done this will define
#
#  MPREAL_FOUND - system has MPFR C++ interface with correct version
#  MPREAL_INCLUDES - the MPFR C++ include directory
#  MPREAL_VERSION - MPFR C++ version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2010 Jitse Niesen, <jitse@maths.leeds.ac.uk>
# Copyright (c) 2015 Jack Poulson, <jack.poulson@gmail.com>
# Redistribution and use is allowed according to the terms of the BSD license.

find_path(MPREAL_INCLUDES NAMES mpreal.h PATHS ${INCLUDE_INSTALL_DIR})

# Set MPREAL_FIND_VERSION to 1.0.0 if no minimum version is specified
if(NOT MPREAL_FIND_VERSION)
  if(NOT MPREAL_FIND_VERSION_MAJOR)
    set(MPREAL_FIND_VERSION_MAJOR 1)
  endif()
  if(NOT MPREAL_FIND_VERSION_MINOR)
    set(MPREAL_FIND_VERSION_MINOR 0)
  endif()
  if(NOT MPREAL_FIND_VERSION_PATCH)
    set(MPREAL_FIND_VERSION_PATCH 0)
  endif()
  set(MPREAL_FIND_VERSION
    "${MPREAL_FIND_VERSION_MAJOR}.${MPREAL_FIND_VERSION_MINOR}.${MPREAL_FIND_VERSION_PATCH}")
endif()

if(MPREAL_INCLUDES)
	# Query MPREAL_VERSION
  file(READ "${MPREAL_INCLUDES}/mpreal.h" _mpreal_version_header)

  string(REGEX MATCH "define[ \t]+MPREAL_VERSION_MAJOR[ \t]+([0-9]+)"
    _mpreal_major_version_match "${_mpreal_version_header}")
  set(MPREAL_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPREAL_VERSION_MINOR[ \t]+([0-9]+)"
    _mpreal_minor_version_match "${_mpreal_version_header}")
  set(MPREAL_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPREAL_VERSION_PATCHLEVEL[ \t]+([0-9]+)"
    _mpreal_patchlevel_version_match "${_mpreal_version_header}")
  set(MPREAL_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(MPREAL_VERSION
    ${MPREAL_MAJOR_VERSION}.${MPREAL_MINOR_VERSION}.${MPREAL_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum required
  if(${MPREAL_VERSION} VERSION_LESS ${MPREAL_FIND_VERSION})
    set(MPREAL_VERSION_OK FALSE)
	message(STATUS "MPFR C++ version ${MPREAL_VERSION} found in ${MPREAL_INCLUDES}, "
                   "but at least version ${MPREAL_FIND_VERSION} is required")
  else()
    set(MPREAL_VERSION_OK TRUE)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPREAL DEFAULT_MSG MPREAL_INCLUDES MPREAL_VERSION_OK)
if (MPREAL_FOUND)
	message(STATUS "MPFR C++ version: ${MPREAL_VERSION}")
endif (MPREAL_FOUND)
mark_as_advanced(MPREAL_INCLUDES)

