include_directories(.)
set(LIBS Axiton ${CMAKE_THREAD_LIBS_INIT} )
LIST(APPEND LIBS ${HDF5_LIBRARY})
LIST(APPEND LIBS ${CUDA_cufft_LIBRARY} ${CUDA_cudart_LIBRARY})

add_executable(axiton axiton.cpp)

target_link_libraries(axiton ${LIBS} )
