#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <errno.h>
#include "enum-vars.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include "cudaErrors.h"

static std::unordered_map<void *, std::pair<AllocType, size_t>> allocTable;
static size_t trackDeviceMem = 0;
static size_t trackAllocMem  = 0;

void	trackFree (void *ptr)
{
	AllocType aType = allocTable[ptr].first;
	size_t    bytes = allocTable[ptr].second;

	switch (aType) {
		case TrackHostAlloc:
		free (ptr);
		trackAllocMem  -= bytes;
		break;

		case TrackDeviceAlloc:
		cudaFree (ptr);
		trackDeviceMem -= bytes;
		break;

		default:
		printf ("Wrong alloc type\n");
		break;
	}

	allocTable.erase(ptr);
	ptr = nullptr;
}

void	trackAlloc (void **ptr, size_t size)
{
	if (((*ptr) = malloc(size)) == nullptr)
	{
		printf ("Error allocating %lu bytes of host memory\n", size);
		exit (1);
	}

	allocTable.insert(std::make_pair(*ptr, std::make_pair(TrackHostAlloc, size)));
	trackAllocMem += size;
}

void	deviceAlloc (void **ptr, size_t size)
{
	if ((cudaMalloc(ptr, size)) != cudaSuccess)
	{
		printf ("Error allocating %lu bytes of device memory\n", size);
		exit (1);
	}

	allocTable.insert(std::make_pair(*ptr, std::make_pair(TrackDeviceAlloc, size)));
	trackDeviceMem += size;
}

void	printMemStats	()
{
	printf ("Total allocated host memory   %lu\n", trackAllocMem);
	printf ("Total allocated device memory %lu\n", trackDeviceMem);

	printf ("\nCurrent pointers in memory:\n");
	printf ("\tHost\n");

	for (auto &data : allocTable) {
		void *ptr   = data.first;
		size_t size = data.second.second;
		if (data.second.first == TrackHostAlloc)
			printf ("Pointer %p\tSize %lu\n", ptr, size);
	}

	printf ("\n\tDevice\n");

	for (auto &data : allocTable) {
		void *ptr   = data.first;
		size_t size = data.second.second;
		if (data.second.first == TrackDeviceAlloc)
			printf ("Pointer %p\tSize %lu\n", ptr, size);
	}
	printf ("\n");
}
