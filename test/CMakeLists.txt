add_executable(libkepler_test
    numtest.c
    testlist.c
    kepler_test.c
    )
target_link_libraries(libkepler_test kepler m)
