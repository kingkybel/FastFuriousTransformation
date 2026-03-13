add_executable(run_tests
    run_tests.cc
    fft_tests.cc
    hough_transform_tests.cc
    spectral_filter_tests.cc
)

target_link_libraries(run_tests
    FastFuriousTransformation
    gtest
    gtest_main
)
add_test(NAME run_tests COMMAND run_tests)
