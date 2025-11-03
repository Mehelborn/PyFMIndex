cmake_minimum_required(VERSION 3.21)
project(
        AWFMIndex
        HOMEPAGE_URL https://github.com/TravisWheelerLab/AvxWindowFmIndex
        LANGUAGES C
)

set(CMAKE_C_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#add support for OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY src/DNAFMIndex/lib)

# RPATH settings for relocatable libraries
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_INSTALL_RPATH "@loader_path")
    set(DETECTED_RPATH "@loader_path")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
    set(DETECTED_RPATH "$ORIGIN")
endif()

#PUBLIC_H_FILES is missing divsufsort64.h, because it doesn't exist until after divsufsort is built
set(
    PUBLIC_H_FILES
    lib/AvxWindowFmIndex/src/AwFmIndex.h
    lib/AvxWindowFmIndex/lib/FastaVector/src/FastaVector.h
    lib/AvxWindowFmIndex/lib/FastaVector/src/FastaVectorMetadataVector.h
    lib/AvxWindowFmIndex/lib/FastaVector/src/FastaVectorString.h
)

set(
    C_FILES
    lib/AvxWindowFmIndex/src/AwFmCreate.c
    lib/AvxWindowFmIndex/src/AwFmFile.c
    lib/AvxWindowFmIndex/src/AwFmIndexStruct.c
    lib/AvxWindowFmIndex/src/AwFmKmerTable.c
    lib/AvxWindowFmIndex/src/AwFmLetter.c
    lib/AvxWindowFmIndex/src/AwFmOccurrence.c
    lib/AvxWindowFmIndex/src/AwFmParallelSearch.c
    lib/AvxWindowFmIndex/src/AwFmSearch.c
    lib/AvxWindowFmIndex/src/AwFmSimdConfig.c
    lib/AvxWindowFmIndex/src/AwFmSuffixArray.c
)

add_library(awfmindex SHARED ${C_FILES})

set_target_properties(
    awfmindex
    PROPERTIES
    PUBLIC_HEADER "${PUBLIC_H_FILES}"
    INSTALL_RPATH "${DETECTED_RPATH}"
    BUILD_WITH_INSTALL_RPATH TRUE
)

target_compile_options(
    awfmindex
    PRIVATE
    -mtune=native
    -march=native
    -Wall
    -Wextra
    -O3
    -fPIC
    -lgomp
    ${OpenMP_C_FLAGS}
)

# Check if the target architecture is x86_64 (Intel) or aarch64 (ARM)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    # Check if the compiler supports AVX2
    include(CheckCCompilerFlag)
    CHECK_C_COMPILER_FLAG("-mavx2" COMPILER_SUPPORTS_AVX2)
    if(COMPILER_SUPPORTS_AVX2)
        # Set AVX2 flags for x86_64 (Intel)
        target_compile_options(awfmindex PRIVATE -mavx2)
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    # Set NEON flags for ARM
    target_compile_options(awfmindex PRIVATE -march=armv8-a+simd)
endif()

install(TARGETS awfmindex)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries" FORCE)
set(CMAKE_BUILD_TYPE "Release")
add_subdirectory(lib/AvxWindowFmIndex/lib)

target_include_directories(
    awfmindex PRIVATE
    lib/AvxWindowFmIndex/lib/FastaVector/src
    lib/AvxWindowFmIndex/lib/libdivsufsort/include
)
target_link_libraries(awfmindex PRIVATE fastavector divsufsort64)

# Copy dependencies to build directory for convenience
add_custom_command(
    TARGET awfmindex POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:fastavector>
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:divsufsort64>
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE_NAME:divsufsort64>
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libdivsufsort64.so.3
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE_NAME:divsufsort64>
        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libdivsufsort64.so
    COMMENT "Copying dependencies to build directory"
)
