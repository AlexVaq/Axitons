include_directories(.)
set(LIBS Axitons ${CMAKE_THREAD_LIBS_INIT} )
LIST(APPEND LIBS ${FFTW_LIBRARY})
LIST(APPEND LIBS ${FFTW_OMP_LIBRARY})
LIST(APPEND LIBS ${HDF5_LIBRARIES})
LIST(APPEND LIBS ${CUDA_cufft_LIBRARY} ${CUDA_cudart_LIBRARY})

add_executable(axiton axiton.cpp)

target_link_libraries(axiton ${LIBS} )
