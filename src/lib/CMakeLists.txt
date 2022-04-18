file(GLOB hdrs "*.h*") 
file(GLOB srcs "*.cpp")

add_library(${lib} STATIC ${srcs} ${hdrs})
target_compile_features(${lib} PRIVATE cxx_std_17)