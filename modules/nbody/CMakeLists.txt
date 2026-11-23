add_library(nbody STATIC
  src/ForceCalculator.cpp
  src/Integrator.cpp
)

target_include_directories(nbody PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(nbody PRIVATE -Wall -Wextra -Wpedantic)
endif()