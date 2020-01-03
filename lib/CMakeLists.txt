set(OBJS utils/initGpu.cpp utils/parse.cpp utils/memAlloc.cpp fields/fields.cpp gen/generator.cpp propagator/propagator.cpp propagator/propKernel.cu io/Hdf5IO.cpp)

include_directories(.)
cuda_add_library(Axitons STATIC ${OBJS})
