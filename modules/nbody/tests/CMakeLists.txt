include(FetchContent)

# GoogleTest is not assumed to be installed on the host. Fetch a pinned release
# so builds are reproducible. Remove this block and link an installed GTest
# instead if a system installation is preferred.
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

add_executable(nbody_tests
  vector3_tests.cpp
  particle_system_tests.cpp
  gravity_tests.cpp
  integrator_tests.cpp
  physical_tests.cpp
)

target_link_libraries(nbody_tests PRIVATE nbody GTest::gtest_main)

gtest_discover_tests(nbody_tests)