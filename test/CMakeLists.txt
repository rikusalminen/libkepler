link_libraries(kepler)
link_libraries(m)

add_executable(libkepler_test
    numtest.c
    testlist.c
    kepler_test.c
    universal_var_test.c
    )

add_executable(cs_series cs_series.c)
