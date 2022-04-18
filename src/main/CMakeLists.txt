add_executable(main main.cpp)
target_include_directories(main PRIVATE ../lib)
target_link_libraries(main PRIVATE ${lib})
target_compile_features(main PRIVATE cxx_std_17)
