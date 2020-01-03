set(OBJS utils/initGpu.cpp utils/parse.cpp utils/memAlloc.cpp fields/fields.cpp propagator/propagator.cpp propagator/propKernel.cu)

include_directories(.)
cuda_add_library(Axitons STATIC ${OBJS})
