#include <cstdio>

#include <cuda.h>
#include <cuda_runtime.h>

#include <omp.h>

// TODO Compute ranks per processor with hwloc (total cache size)

static int nThreads      =  1;

static int tPerBlock	 = 0;
static int maxThreads[3] = { 0, 0, 0 };
static int maxGrid[3]    = { 0, 0, 0 };

static size_t gpuMem     = 0;

int	maxThreadsPerBlock() {
    return	tPerBlock;
}

int	maxThreadsPerDim(const int dim) {
    return	(dim < 3) ? maxThreads[dim] : 0;
}

int	maxGridSize(const int dim) {
    return	(dim < 3) ? maxGrid[dim] : 0;
}

size_t	gpuMemAvail() {
    return	gpuMem;
}

int	initGpu (int accId = 0) {
    int nAccs = 0;

    cudaError_t cErr = cudaGetDeviceCount(&nAccs);

    if (cErr != cudaSuccess) {
        printf ("CUDA error: %s", cudaGetErrorString(cErr));
        return	-1;
    }

    printf ("Got %d accelerators, selecting accelerator %d\n", nAccs, accId);

    cudaSetDevice(accId);
    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

    cudaDeviceProp gpuProp;
    cudaGetDeviceProperties(&gpuProp, accId);

    printf ("  Peak Memory Bandwidth of Gpu %d (GB/s): %f\n", accId, 2.0*gpuProp.memoryClockRate*(gpuProp.memoryBusWidth/8)/1.0e6);
    gpuMem	      = gpuProp.totalGlobalMem;
    tPerBlock     = gpuProp.maxThreadsPerBlock;
    maxThreads[0] = gpuProp.maxThreadsDim[0];
    maxThreads[1] = gpuProp.maxThreadsDim[1];
    maxThreads[2] = gpuProp.maxThreadsDim[2];
    maxGrid[0]    = gpuProp.maxGridSize[0];
    maxGrid[1]    = gpuProp.maxGridSize[1];
    maxGrid[2]    = gpuProp.maxGridSize[2];

    int nProcs, mThreads;

    #pragma omp parallel
    {
        nProcs   = omp_get_num_procs();
        nThreads = omp_get_num_threads();
        mThreads = omp_get_max_threads();
    }

    printf ("Cpu will use %d threads for %d processors (max %d)\n", nThreads, nProcs, mThreads);

    return	0;
}
