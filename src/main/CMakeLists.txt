add_executable(main main.cpp)
target_include_directories(main PRIVATE ../lib)
target_link_libraries(main PRIVATE ${lib})
target_compile_features(main PRIVATE cxx_std_17)


find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        f23347fca25219d1c42cbb91608b5556814bf572
)

FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()

include_directories(SYSTEM ${matplotlib_SOURCE_DIR})
target_link_libraries(main
    PRIVATE
        ${PYTHON_LIBRARIES}
        Python3::NumPy
)