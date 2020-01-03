#ifndef	InitializeGpuGuard
	#define	InitializeGpuGuard

	int	initGpu(int accId = 0);
	size_t	gpuMemAvail();
	int	maxThreadsPerBlock();
	int	maxThreadsPerDim(const int dim);
	int	maxGridSize(const int dim);
#endif
