cmake_minimum_required (VERSION 3.0.0)

project (mappability CXX)
set (CXX11_FOUND "1")

#set (CMAKE_CXX_COMPILER "g++-6")

if (NOT DEFINED ${CMAKE_BUILD_TYPE})
    set (CMAKE_BUILD_TYPE Release)
    message("*** CMAKE_BUILD_TYPE was not passed to CMake, therefore it was set to Release. ***")
endif()

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan/include")

# Load the SeqAn module and fail if not found.
# Search for zlib as a dependency for SeqAn.
find_package (SeqAn REQUIRED)
find_package (OpenMP)

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "-std=c++14 -DSEQAN_DISABLE_VERSION_CHECK=YES -I ${CMAKE_CURRENT_SOURCE_DIR}/sdsl-lite/include ${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.
add_executable (mappability mappability.cpp)
target_link_libraries (mappability ${SEQAN_LIBRARIES})

add_executable (tests tests.cpp)
target_link_libraries (tests ${SEQAN_LIBRARIES})

add_executable (create_index create_index.cpp)
target_link_libraries (create_index ${SEQAN_LIBRARIES})

add_executable (mappability_dump mappability_dump.cpp)
target_link_libraries (mappability_dump ${SEQAN_LIBRARIES})

add_executable (heatmap heatmap.cpp)
target_link_libraries (heatmap ${SEQAN_LIBRARIES})

add_executable (cut_n cut_n.cpp)
target_link_libraries (cut_n ${SEQAN_LIBRARIES})

add_executable (insert_n insert_n.cpp)
target_link_libraries (insert_n ${SEQAN_LIBRARIES})
