include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)

pkg_search_module(_LIBXC libxc>=${Libxc_FIND_VERSION})

find_library(LIBXC_LIBRARIES NAMES xc
  PATH_SUFFIXES lib
  HINTS
  ENV EBROOTLIBXC
  ENV LIBXCROOT
  ${_LIBXC_LIBRARY_DIRS}
  DOC "libxc libraries list")

find_library(LIBXC_LIBRARIES_F03 NAMES xcf03
  PATH_SUFFIXES lib
  HINTS
  ENV EBROOTLIBXC
  ENV LIBXCROOT
  ${_LIBXC_LIBRARY_DIRS}
  DOC "libxc libraries list")

find_path(LIBXC_INCLUDE_DIR NAMES xc.h xc_f90_types_m.mod xc_f03_lib_m.mod
  PATH_SUFFIXES inc include
  HINTS
  ${_LIBXC_INCLUDE_DIRS}
  ENV EBROOTLIBXC
  ENV LIBXCROOT
  )

find_package_handle_standard_args(Libxc DEFAULT_MSG LIBXC_LIBRARIES LIBXC_INCLUDE_DIR)

if (Libxc_FOUND)
  set(Libxc_VERSION ${Libxc_FIND_VERSION})
  set(Libxc_INCLUDE_DIR ${LIBXC_INCLUDE_DIR})
  add_library(Libxc::xcf03 INTERFACE IMPORTED)
  set_target_properties(Libxc::xcf03 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${LIBXC_LIBRARIES};${LIBXC_LIBRARIES_F03}")
endif()

