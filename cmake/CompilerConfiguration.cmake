if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS          "-ffree-form -ffree-line-length-none -std=f2008ts")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-O3 -funroll-loops")
  set(CMAKE_Fortran_FLAGS_DEBUG    "-O0 -ggdb")
  set(CMAKE_Fortran_FLAGS_COVERAGE "-O0 -fprofile-arcs -ftest-coverage")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS          "-free -stand f08 -fpp")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-O3 -diag-disable=5268")  # Disable the line-length-extension warning #5268
  set(CMAKE_Fortran_FLAGS_DEBUG    "-O0 -debug")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_Fortran_FLAGS          "-Mfreeform -Mextend -Mallocatable=03")  # -Mallocatable=03: enable F2003+ assignment semantics
  set(CMAKE_Fortran_FLAGS_RELEASE  "-fast")
  set(CMAKE_Fortran_FLAGS_DEBUG    "-g")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set(CMAKE_Fortran_FLAGS          "-f2008 -free -Warn=reallocation -Warn=subnormal")
  set(CMAKE_Fortran_FLAGS_RELEASE  "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG    "-g -C")
  if (NOT OpenMP_FOUND)
    set(CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_Fortran_FLAGS_RELEASE} -gline")  # -gline is only supported without OpenMP
    set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_Fortran_FLAGS_DEBUG} -C=all")  # some checks are not available with OpenMP
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_Fortran_FLAGS          "-f free -M 3105")  # -M 3105: hide a false-positive warning about modified loop variables due to loop fusing
  set(CMAKE_Fortran_FLAGS_RELEASE  "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG    "-g -h error_on_warning")
  set(CMAKE_Fortran_MODOUT_FLAG    "-ef")  # override to get lower-case module file names
else ()
  message(WARNING "\
Unknown Fortran compiler, trying without any additional (optimization) flags.\n\
You will most likely have to specify extra options for free-form Fortran 2008 for your compiler!\n\
Please open an issue at https://github.com/cp2k/dbcsr/issues with the reported compiler name and the required flags.")
  message("-- CMAKE_Fortran_COMPILER_ID: " ${CMAKE_Fortran_COMPILER_ID})
  message("-- CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
endif ()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS              "-std=c++11")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -funroll-loops")
  set(CMAKE_CXX_FLAGS_COVERAGE     "-O0 -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG        "-O0 -ggdb")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS              "-std=c++11")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -funroll-loops")
  set(CMAKE_CXX_FLAGS_DEBUG        "-O0 -g")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(CMAKE_CXX_FLAGS              "-std=c++11")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3 -funroll-loops")
  set(CMAKE_CXX_FLAGS_DEBUG        "-O0 -g")

  if (CMAKE_BUILD_TYPE MATCHES Coverage)
    # when building with coverage suppport, shared libs/executables must be explicitly linked to avoid undefined symbols
    SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lgcov --coverage")
  endif ()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3")
  set(CMAKE_CXX_FLAGS_DEBUG        "-O0 -debug")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_CXX_FLAGS_RELEASE      "-fast")
  set(CMAKE_CXX_FLAGS_DEBUG        "-g")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Cray")
  set(CMAKE_CXX_FLAGS_RELEASE      "-O3")
  set(CMAKE_CXX_FLAGS_DEBUG        "-g")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    # prevent deallocation failures due to tcmalloc's free with glibc's aligned_alloc, see https://bugzilla.redhat.com/show_bug.cgi?id=1569391
    set(CMAKE_C_FLAGS              "${CMAKE_C_FLAGS} -h system_alloc")
    set(CMAKE_CXX_FLAGS            "${CMAKE_CXX_FLAGS} -h system_alloc")
    set(CMAKE_Fortran_FLAGS        "${CMAKE_Fortran_FLAGS} -h system_alloc")
    # since the detection of the implicitly linked libraries occurs bevore we can intervene, filter them out again
    list(FILTER CMAKE_C_IMPLICIT_LINK_LIBRARIES EXCLUDE REGEX "tcmalloc")
    list(FILTER CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES EXCLUDE REGEX "tcmalloc")
  endif ()
else ()
  message(WARNING "\
Unknown C++ compiler, trying without any additional (optimization) flags.\n\
You may have to specify flags for C++11/14 support manually.\n\
Please open an issue at https://github.com/cp2k/dbcsr/issues with the reported compiler name and the required flags.")
  message("-- CMAKE_CXX_COMPILER_ID: " ${CMAKE_CXX_COMPILER_ID})
  message("-- CMAKE_CXX_COMPILER full path: " ${CMAKE_CXX_COMPILER})
endif ()
