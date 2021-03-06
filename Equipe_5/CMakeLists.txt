cmake_minimum_required(VERSION 2.8)
project(mc-pricer CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


# Detect Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})

# Add math lib
set(LIBS ${LIBS} m)

# Testing the parser
add_executable(parser-test src/parser.cpp src/test_parser.cpp)
add_executable(opt-test src/basket.cpp src/barrier_l.cpp src/barrier_u.cpp src/asian.cpp src/barrier.cpp src/performance.cpp src/parser.cpp src/test_options.cpp) 
add_executable(bs-test src/basket.cpp src/barrier_l.cpp src/barrier_u.cpp src/asian.cpp src/barrier.cpp src/performance.cpp src/bs.cpp src/parser.cpp src/test_bs.cpp)
add_executable(pricer src/basket.cpp src/barrier_l.cpp src/barrier_u.cpp src/asian.cpp src/barrier.cpp src/performance.cpp src/bs.cpp src/mc.cpp src/parser.cpp src/pricer.cpp)
add_executable(generate-past src/basket.cpp src/barrier_l.cpp src/barrier_u.cpp src/asian.cpp src/barrier.cpp src/performance.cpp src/bs.cpp src/parser.cpp src/past.cpp)
target_link_libraries(parser-test ${LIBS})
target_link_libraries(opt-test ${LIBS})
target_link_libraries(bs-test ${LIBS})
target_link_libraries(pricer ${LIBS})
target_link_libraries(generate-past ${LIBS})
install(TARGETS parser-test DESTINATION bin)
install(TARGETS opt-test DESTINATION bin)
install(TARGETS bs-test DESTINATION bin)
install(TARGETS pricer DESTINATION bin)
install(TARGETS generate-past DESTINATION bin)
# add_executable(mc-pricer src/parser.cpp src/pricer.cpp)
# target_link_libraries(mc-pricer ${LIBS})
# install(TARGETS mc-pricer DESTINATION bin)


